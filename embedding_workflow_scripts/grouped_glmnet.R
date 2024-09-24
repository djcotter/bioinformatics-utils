# grab_top_embeddings_by_variance.R
# Daniel Cotter
# 2024-09-18

# This script takes in an ordering file and an embeddings file and decomposes it into clusters. 
# it then runs a glmnet on each cluster for a given metadata category to determine a smaller
# set of nonzero coeficients and then runs a larger glmnet on the entire set of nonzero embeddings
# from each cluster to determine the most important embeddings overall for that metadata category


## import packages --------
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(resample))
suppressPackageStartupMessages(library(furrr))


## parse arguments --------
# define command line arguments
# define command line arguments
option_list <- list(
  make_option(c("-m", "--metadata"), help="Metadata file", type="character"),
  make_option(c("-e", "--embeddings"), help="Embeddings file", type="character"),
  make_option(c("-o", "--ordering"), help="Ordering file", type="character"),
  make_option(c("-p", "--output_prefix"), help="Output prefix.", type="character"),
  make_option(c("--temp_dir"), help="Temporary directory to store intermediate files", 
              type="character"),
  make_option(c("--num_threads"), help="Number of threads to use for parallel operations",
              type="integer", default=1),
  make_option(c("--metadata_categories"), help="Metadata categories to run glmnet on. Comma separated", 
              type="character"),
  make_option(c("--even_classes"), help="Sample equal numbers of each class", action="store_true", default=FALSE),
  make_option(c("-s", "--min_samples_per_category"), help="Minimum number of samples per metadata category", 
              type="integer", default = 30),
)

# parse command line arguments
opt <- parse_args(OptionParser(option_list = option_list))

# check that user specified all files
if (!file.exists(opt$metadata) | !file.exists(opt$embeddings) | !file.exists(opt$ordering) | is.null(opt$output_prefix)) {
  stop("Must provide embeddings, ordering, and output prefix")
}

# create a temporary directory to store intermediate files
if (!is.null(opt$temp_dir)) {
  temp_dir <- ifelse(grepl("/$", opt$temp_dir),
                     opt$temp_dir, 
                     paste0(opt$temp_dir, "/"))
  system(paste("mkdir -p", temp_dir))
} else {
  temp_dir <- file.path(dirname(opt$output_prefix), "tmp/")
  system(paste("mkdir -p", temp_dir))
}

# define future plans
setDTthreads(opt$num_threads)
plan(multicore, workers=opt$num_threads)
options(future.globals.maxSize=4000*1024^2)
min_num_per_category <- opt$min_samples_per_category

## print a summary of the arguments
cat("\n####################\n")
cat("Running grab_top_embeddings_by_variance.R with the following arguments:\n")
cat("Metadata file: ", opt$metadata, "\n")
cat("Ordering file: ", opt$ordering, "\n")
cat("Embeddings file: ", opt$embeddings, "\n")
cat("Output file: ", opt$output_prefix, "\n")
cat("Num embeddings to keep per anchor: ", opt$num_to_keep, "\n")
cat("Temporary directory: ", temp_dir, "\n")
cat("####################\n\n")

## load data --------
# load embeddings
cat("\nLoading embeddings...\n")
# copy the embeddings file to the temp directory to speed up I/O
embeddings_temp <- file.path(temp_dir, "raw_embeddings_temp.tsv")
if (!file.exists(embeddings_temp)) {
  system(paste("cp", opt$embeddings, embeddings_temp))
}
embeddings <- fread(embeddings_temp, header = F)
colnames(embeddings) <- c("kmer", paste0("embedding_", 1:(ncol(embeddings)-1)))

# load the ordering file
cat("Loading the ordering file...\n")
# copy the ordering file to the temp directory to speed up I/O
ordering_temp <- file.path(temp_dir, "ordering_temp.tsv")
if (!file.exists(ordering_temp)) {
  system(paste("cp", opt$ordering, ordering_temp))
}
ordering <- fread(ordering_temp, header=F, sep="\t", 
                  col.names = c("sample_name", "seq", "kmer", "start", "end")) 
ordering <- ordering %>% select(sample_name, kmer)

ordering <- ordering %>% group_by(sample_name) %>% mutate(cluster=row_number()) %>%
  mutate(cluster=cluster-1) %>% ungroup()

clusters <- ordering %>% group_by(cluster) %>% group_split()

join_and_write_clusters <- function(cluster_df, filename, all_embeddings) {
  cluster_df %>% select(-cluster) %>% left_join(all_embeddings, by="kmer") %>%
    select(-kmer) %>% fwrite(file=filename,nThread = 1,col.names = T)
}

temp_embeddings_dir <- file.path(temp_dir, "embeddings_per_cluster/")
system(paste("mkdir -p", temp_embeddings_dir))
cluster_files <- paste0(temp_embeddings_dir, "embeddings_cluster_", 0:(length(clusters)-1), ".csv")

cat("Writing all clusters and their embeddings out to file in: ", temp_embeddings_dir)

if (!sum(file.exists(cluster_files)) == length(cluster_files)) {
  future_walk2(clusters, cluster_files, \(x,y) join_and_write_clusters(x, y, all_embeddings = embeddings), .progress = T)
}

# read in the metadata file
# load metadata
cat("Loading metadata...\n")
my_metadata <- fread(opt$metadata, header =T)
# if sample_name is not in the metadata rename the first column to sample name
if (!"sample_name" %in% colnames(my_metadata)) {
  cat("Cannot find sample_name in metadata. Renaming first column to sample_name...\n")
  setnames(my_metadata, names(my_metadata)[1], "sample_name")
}

# filter through metadata columns and decide which to use 
# change this to fit your specific metadata
if (!is.null(opt$metadata_categories)) {
  metadata_labels <- strsplit(opt$metadata_categories, ",") %>% unlist()
} else {
  metadata_labels <- NULL
  for (i in colnames(my_metadata %>% select(-sample_name))) {
    num_cats <- my_metadata %>% group_by(get(i)) %>% filter(n() > min_num_per_category) %>% 
      filter(!is.na(get(i))) %>% pull(get(i)) %>% n_distinct()
    if (num_cats >= 2) {
      metadata_labels = c(metadata_labels, i)
    }
  }
}

# print the metadata labels that are being used
cat("Using the following metadata labels:\n\t", paste(metadata_labels, collapse = "\n \t"), "\n")


# define a function to grab
grab_glmnet_significant_features <- function(in_file, metadata_label, metadata) {
  cluster_num = str_extract(in_file, "cluster_(\\d+).csv", group=1) %>% as.integer()
  temp_dt <- fread(in_file, header=T, nThread = 1) %>% select(sample_name, starts_with("embedding")) # first filter for only one cluster
  colnames(temp_dt) <- ifelse(grepl("embedding", colnames(temp_dt)),
                      yes=paste0("cluster_", cluster_num, "_", colnames(temp_dt)),
                      no = colnames(temp_dt))
  glm_dt <- temp_dt %>% left_join(metadata, by="sample_name") %>% select(-sample_name)
  glm_dt <- glm_dt %>% rename(class=!!metadata_label) %>% filter(!is.na(class))
  glm_dt <- glm_dt %>% group_by(class) %>% filter(n() > min_num_per_category) %>% ungroup()
  X <- as.matrix(glm_dt %>% select(starts_with("cluster_")))
  y <- as.factor(glm_dt$class)
  cv_fit <- tryCatch(cv.glmnet(X, y, family="multinomial", type.measure="class"), error=function(e) NULL)
  if (is.null(cv_fit)) {
    cv_fit <- tryCatch(cv.glmnet(X, y, family="multinomial", type.measure="class"), error=function(e) NULL)
  }
  if (is.null(cv_fit)) {
    return(NULL)
  }
  coef <- coef(cv_fit, s="lambda.min")
  coefs <- map2_dfc(coef,names(coef), \(x,y) as.matrix(x) %>% as.data.frame() %>% setnames(y)) %>% 
    rownames_to_column("feature") %>% filter(!grepl("Intercept", feature)) %>% 
    filter(if_any(2:ncol(.), ~.x != 0)) %>% pull(feature)
  
  temp_dt <- temp_dt %>% select(sample_name, all_of(coefs)) %>% arrange(sample_name)
  if (cluster_num!=0) {
    temp_dt <- temp_dt %>% select(-sample_name)
  }
  return(temp_dt)
}

for (i in metadata_labels) {
  metadata_label <- i
  cat("Fitting glmnet model for metadata label: ", i, "\n")
  
  # grabbing the top glmnet features for each cluster
  cat("Calculating top glmnet features for each cluster...\n")
  top_glm_dt <- future_map_dfc(cluster_files,
                               \(x) grab_glmnet_significant_features(x, metadata_label = i, metadata=my_metadata),
                               .progress = T)
  top_glm_dt <- top_glm_dt %>% relocate(sample_name)
  
  # add the metadata back onto the embeddings
  curr_metadata <- my_metadata %>% select(sample_name, all_of(i))
  top_glm_dt <- top_glm_dt %>% 
    left_join(curr_metadata, by="sample_name") %>% 
    rename(class=!!i)

  dt <- top_glm_dt %>% filter(!is.na(class))
  dt <- my_dt %>% group_by(class) %>% filter(n() > min_num_per_category) %>% ungroup()
  
  if (nrow(dt) < 90) {
    cat(paste("Only", nrow(dt), "observations in data. Skipping", metadata_label, "\n\n"))
    next
  }
  
  if (dt %>% pull(class) %>% unique() %>% length() <2) {
    cat("Not enough metadata categories after filtering. Skipping...\n\n")
    next
  }
  
  dt <- as.data.table(dt)
  print(metadata_label)
  cat("Data Breakdown...\n")
  print(dt[, .N, by=class])
  cat("\n\n")
  # decide on the number of samples to keep in the training set
  # such that the number of samples is approximately equal for each class
  n_samples <- min(dt[, .N, by=class]$N)
  n_train <- floor(n_samples * 0.5)
  cat(paste("Training on", n_train, "samples per category.\n"))
  
  # separate into training and testing
  # keep equal numbers of each class in training and testing
  set.seed(123)
  if (opt$even_classes) {
    cat("Sampling equal numbers of each class...\n")
    # sample n_sample from each class before splitting into train and test
    sample_idx <- dt[, sample(.I, n_samples), by=class]$V1
    dt <- dt[sample_idx]
    train_idx <- dt[, sample(.I, n_train), by=class]$V1
    train <- dt[train_idx]
    test <- dt[-train_idx]
  } else {
    cat("Using remaining samples in each class for testing...\n")
    # sample n_train from each class
    train_idx <- dt[, sample(.I, n_train), by=class]$V1
    train <- dt[train_idx]
    test <- dt[-train_idx]
  }
  
  # prepare data for glmnet
  X <- as.matrix(train[, -c("sample_name", "class")])
  y <- train$class
  
  # fit glmnet
  cat("Fitting cv glmnet model...\n")
  time_fit <- Sys.time()
  fit <- tryCatch(cv.glmnet(X, y, family="multinomial", type.measure="class"), error=function(e) NULL)
  if (is.null(fit)) {
    fit <- tryCatch(cv.glmnet(X, y, family="multinomial", type.measure="class"), error=function(e) NULL)
  }
  if (is.null(fit)) {
    cat("Error in glmnet after trying twice. Skipping...\n\n")
    next
  }
  cat("Model fit in ", round(Sys.time() - time_fit, 2), " seconds.\n")
  cat("Using the model for lambda.min\n")
  
  # predict on test set
  X_test <- as.matrix(test[, -c("sample_name", "class")])
  y_test <- test$class
  y_pred <- predict(fit, X_test, s="lambda.min", type="class")
  
  # confusion matrix
  confusion_matrix <- table(y_test, y_pred)
  cat("Confusion matrix:\n")
  print(confusion_matrix)
  cat("\n\n")
  
  # get the accuracy
  acc <- sum(diag(confusion_matrix)) / sum(confusion_matrix)
  cat(paste("accuracy:", acc))
  cat("\n")
  
  # get the sensitivity (use a trycatch and return NA if there's an error)
  sens <- tryCatch(caret::sensitivity(confusion_matrix), error=function(e) NA)
  
  # get the specificity (use a trycatch and return NA if there's an error)
  spec <- tryCatch(caret::specificity(confusion_matrix), error=function(e) NA)
  
  # extract non-zero coefficients
  coef <- coef(fit, s="lambda.min")
  
  # transform multinomial coefficients into data frame
  coef_df <- map2_dfc(coef,names(coef), \(x,y) as.matrix(x) %>% as.data.frame() %>% setnames(y)) %>% 
    rownames_to_column("feature") %>% 
    filter(if_any(2:ncol(.), ~.x != 0)) %>%
    unite("coefficients", -feature, sep=",") %>%
    mutate(coefficients=paste0("[", coefficients, "]")) %>%
    mutate(classes=paste0("[", paste(names(coef), collapse=","), "]")) %>%
    mutate(metadata_category=metadata_label, accuracy = acc,
           sensitivity=sens, specificity=spec) %>%
    select(metadata_category, feature, accuracy,
           sensitivity, specificity, classes, coefficients)
  
  # plot confusion matrix
  confusion_matrix <- table(y_test, y_pred)
  confusion_matrix <- as.data.frame.matrix(confusion_matrix)
  confusion_matrix <- rownames_to_column(confusion_matrix, "true_class")
  confusion_matrix <- gather(confusion_matrix, "predicted_class", "feature_count", -true_class)
  confusion_matrix <- confusion_matrix %>% mutate(true_class = factor(true_class, levels=unique(true_class)),
                                                  predicted_class = factor(predicted_class, levels=unique(predicted_class)))
  
  p <- ggplot(confusion_matrix, aes(x=true_class, y=predicted_class, fill=feature_count)) +
    geom_tile(color="white") +
    scale_fill_gradient(low="white", high="steelblue") +
    theme_light() +
    theme(axis.text.x = element_text(angle=45, hjust=1)) +
    labs(x="True class", y="Predicted class", fill="Count") 
  
  # add accuracy and all counts to plot
  p <- p + geom_text(aes(label=feature_count), vjust=1, size=3) +
    ggtitle(paste0(metadata_label, " | Accuracy: ", round(acc, 2),
                   "\nSensitivity: ", round(sens, 2),
                   " | Specificity: ", round(spec, 2)))
  print(p)
  all_coef = rbind(all_coef, coef_df)
  cat("\n\n")
}

for (i in metadata_labels) {
  top_glm_dt <- future_map_dfc(cluster_files,
                               \(x) grab_glmnet_significant_features(x, metadata_label = i),
                               .progress = T)
  top_glm_dt <- top_glm_dt %>% relocate(sample_name)
  
  # perform the larger glmnet
}

cat("Calculating top variance components per cluster...\n")
top_glm_dt <- future_map_dfc(cluster_files,
                             \(x) grab_glmnet_significant_features(x, metadata_label = metadata_labels[1]),
                             .progress = T)
top_glm_dt <- top_glm_dt %>% relocate(sample_name)

# write out the embeddings matrix to a temp file
embeddings_metadta_file <- paste0(opt$output_prefix, "_top_glm_embeddings.tsv")
cat("Writing top variance embeddings to ", embeddings_topVar_file, "\n")
write_tsv(top_var_dt, embeddings_topVar_file, col_names = T, quote="needed")
embeddings_feather <- file.path(temp_dir, "top_variance_embeddings.feather")
feather::write_feather(top_var_dt, embeddings_feather)
