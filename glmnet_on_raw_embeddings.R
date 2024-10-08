# glmnet_on_embeddings_PCs.R
# Daniel Cotter
# 2024-09-15

# This script takes in three files. 1) a metadata file, 2) a tsv of kmer embeddings, 
# and 3) a tsv ordering file mapping each kmer to it's position in a sample. It then calculates 
# PCA on the embeddings and maps N PCs to the ordering file. The script then fits a glmnet 
# model to the samples using the embedding PCs as features and each metadata category for a
# different response. The output is a file containing the nonzero coefficients for each model
# and a pdf of the confusion matrix for each model. 


## import packages --------
suppressPackageStartupMessages(library(glmnet))
suppressPackageStartupMessages(library(caret))
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(optparse))


## parse arguments --------
# define command line arguments
option_list <- list(
  make_option(c("-m", "--metadata"), help="Metadata file", type="character"),
  make_option(c("-e", "--embeddings"), help="Embeddings file", type="character"),
  make_option(c("-o", "--ordering"), help="Ordering file", type="character"),
  make_option(c("-n", "--num_embeddings"), help="Number of embeddings per cluster to use", type="integer", default = 10),
  make_option(c("-p", "--output_prefix"), help="Output prefix.", type="character"),
  make_option(c("-s", "--min_samples_per_category"), help="Minimum number of samples per metadata category", 
              type="integer", default = 30),
  make_option(c("--even_classes"), help="Sample equal numbers of each class", action="store_true", default=FALSE),
  make_option(c("--temp_dir"), help="Temporary directory to store intermediate files", 
              type="character")
)


# parse command line arguments
opt <- parse_args(OptionParser(option_list = option_list))

# check that user specified all files
if (!file.exists(opt$metadata) | !file.exists(opt$embeddings) | !file.exists(opt$ordering) | is.null(opt$output_prefix)) {
  stop("Must provide metadata, embeddings, ordering, and output prefix")
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

## specify the output files
coefficients_out = paste0(opt$output_prefix, "_nonzero_coefficients.tsv")
confusion_matrix_out = paste0(opt$output_prefix, "_confusion_matrices.pdf")
min_num_per_category = opt$min_samples_per_category

## print a summary of the arguments
cat("\n####################\n")
cat("Running glmnet_on_embeddings_PCs.R with the following arguments:\n")
cat("Metadata file: ", opt$metadata, "\n")
cat("Embeddings file: ", opt$embeddings, "\n")
cat("Ordering file: ", opt$ordering, "\n")
cat("Number of PCs: ", opt$num_pcs, "\n")
cat("Output coefficients: ", coefficients_out, "\n")
cat("Output confusion matrices: ", confusion_matrix_out, "\n")
cat("Min number of samples per metadata label: ", min_num_per_category, "\n")
cat("Temporary directory: ", temp_dir, "\n")
cat("####################\n\n")

## load data --------
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
metadata_labels <- NULL
for (i in colnames(my_metadata %>% select(-sample_name))) {
  num_cats <- my_metadata %>% group_by(get(i)) %>% filter(n() > min_num_per_category) %>% 
    filter(!is.na(get(i))) %>% pull(get(i)) %>% n_distinct()
  if (num_cats >= 2) {
    metadata_labels = c(metadata_labels, i)
  }
}
# print the metadata labels that are being used
cat("Using the following metadata labels:\n\t", paste(metadata_labels, collapse = "\n \t"), "\n")

# load embeddings
cat("\nLoading embeddings...\n")
# copy the embeddings file to the temp directory to speed up I/O
embeddings_temp <- file.path(temp_dir, "embeddings_temp.tsv")
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

# merge the ordering file with the embeddings PCs
cat("Merging the ordering file with the embeddings...\n")

main_dt <- ordering %>% left_join(embeddings, by="kmer") %>%
  group_by(sample_name) %>% mutate(cluster = row_number()) %>% 
  mutate(cluster=cluster-1) %>% ungroup()

cluster_to_kmer_mapping <- main_dt %>% select(cluster, kmer) %>% distinct()
cluster_to_kmer_mapping <- cluster_to_kmer_mapping %>% 
  arrange(cluster) %>% group_by(cluster) %>% 
  summarise(kmers = str_c(kmer, collapse=",")) %>% ungroup()

cat("Pivoting the embeddings to wide format...\n")
main_dt <- main_dt %>% select(sample_name, cluster, starts_with("embedding"))
main_dt <- data.table::dcast(setDT(main_dt), sample_name ~ cluster, 
                             value.var = paste("embedding", 
                                               1:(ncol(embeddings)-1),
                                               sep="_"), sep="-")

cat("Grabbing the top", opt$num_embeddings, "embeddings by variance per cluster...\n")
col_variances <- resample::colVars(main_dt %>% select(starts_with("embedding"))) %>% 
  enframe() %>% separate(name, into=c("embedding", "cluster"), sep="-")

top_var_cols <- col_variances %>% 
  group_by(cluster) %>% 
  slice_max(value, n=opt$num_embeddings, with_ties=F) %>% ungroup() %>% 
  unite(name, embedding, cluster, sep="-") %>% 
  pull(name)

main_dt <- main_dt %>% select(sample_name, all_of(top_var_cols))

new_names <- colnames(main_dt %>% select(-sample_name)) %>% 
  str_extract("(\\d+)-(\\d+)", group=c(1,2)) %>%
  as.data.frame()
new_names <- paste("cluster_", new_names$V2, "_embedding_", new_names$V1, sep="") %>% 
  c("sample_name", .)

colnames(main_dt) <- new_names

# write out the embeddings matrix to a temp file
embeddings_topVar_file <- file.path(temp_dir, "embeddings_topVar.tsv")
cat("Writing top variance embeddings to ", embeddings_topVar_file, "\n")
write_tsv(main_dt, embeddings_topVar_file)

rm(ordering, embeddings_pcs, embeddings)
gc()

## Fit all glmnet models --------
all_coef <- NULL

pdf(confusion_matrix_out)

for (i in metadata_labels) {
  cat("Fitting glmnet model for metadata label: ", i, "\n")
  # filter metadata down to only samples that have a value for the current metadata label
  all_metadata <- my_metadata %>% filter(!is.na(get(i)))
  
  # get the samples that are present in the embeddings and metadata
  curr_samples <- all_metadata %>% pull(sample_name) %>% unique()
  my_dt <- main_dt %>% filter(sample_name %in% curr_samples)
  
  # filter to make sure the metadata calculations are done on the samples that are actually present
  all_metadata <- all_metadata %>% filter(sample_name %in% my_dt$sample_name)
  my_dt <- my_dt %>% filter(sample_name %in% all_metadata$sample_name)
  
  # pivot main dt wider to have one row per sample
  wide_dt <- my_dt
  
  # rename the metadata column to class for further operations
  metadata_label <- i
  metadata <- all_metadata %>%
    select(sample_name, all_of(metadata_label)) %>%
    rename(class=all_of(metadata_label)) %>%
    mutate(class=ifelse(class=="", NA, class))
  
  # merge metadata with embeddings
  dt <- merge(wide_dt, metadata, by.x="sample_name", by.y="sample_name")
  
  # filter out samples with missing serotype
  dt <- dt %>% filter(!is.na(class))
  
  # filter out serotypes with less than 100 samples
  dt <- dt %>% group_by(class) %>% filter(n() > min_num_per_category)
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
  fit <- cv.glmnet(X, y, family="multinomial", type.measure="class")
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
  confusion_matrix <- gather(confusion_matrix, "predicted_class", "count", -true_class)
  confusion_matrix <- confusion_matrix %>% mutate(true_class = factor(true_class, levels=unique(true_class)),
                                                  predicted_class = factor(predicted_class, levels=unique(predicted_class)))
  
  p <- ggplot(confusion_matrix, aes(x=true_class, y=predicted_class, fill=count)) +
    geom_tile(color="white") +
    scale_fill_gradient(low="white", high="steelblue") +
    theme_light() +
    theme(axis.text.x = element_text(angle=45, hjust=1)) +
    labs(x="True class", y="Predicted class", fill="Count") 
  
  # add accuracy and all counts to plot
  p <- p + geom_text(aes(label=count), vjust=1, size=3) +
    ggtitle(paste0(metadata_label, " | Accuracy: ", round(acc, 2),
                   "\nSensitivity: ", round(sens, 2),
                   " | Specificity: ", round(spec, 2)))
  print(p)
  all_coef = rbind(all_coef, coef_df)
  cat("\n\n")
}

dev.off()

saveRDS(all_coef, file.path(temp_dir, "all_glmnet_coef.RDS"))
# write out the coefficients
# first filter for non-zero coefficients
# then join on the headers for the significant kmers
nonzero_coef <- all_coef %>% filter(!grepl("0,", coefficients)) %>% 
  mutate(cluster = str_extract(feature, "cluster_(\\d+)_", group=1)) %>% 
  mutate(cluster=as.double(cluster)) %>% 
  left_join(cluster_to_kmer_mapping, by="cluster") %>% 
  select(-cluster)

nonzero_coef %>% write_tsv(file=coefficients_out, col_names = T, quote = "needed", na = "")
