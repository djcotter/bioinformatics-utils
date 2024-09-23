# embeddings_glmnet.R
# Daniel Cotter
# 2024-09-15

# Take in a fully prepared embeddings file and metadata file and run glmnet
# output a pdf of confusion matrices and a tsv of the non-zero coefficients


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
  make_option(c("-p", "--output_prefix"), help="Output prefix.", type="character"),
  make_option(c("-s", "--min_samples_per_category"), help="Minimum number of samples per metadata category", 
              type="integer", default = 30),
  make_option(c("--even_classes"), help="Sample equal numbers of each class", action="store_true", default=FALSE),
  make_option(c("--temp_dir"), help="Temporary directory to store intermediate files", 
              type="character")
)

setDTthreads(threads=0L)


# parse command line arguments
opt <- parse_args(OptionParser(option_list = option_list))

# check that user specified all files
if (!file.exists(opt$metadata) | !file.exists(opt$embeddings) | is.null(opt$output_prefix)) {
  stop("Must provide metadata, embeddings, and output prefix")
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
if (grepl(".feather", opt$embeddings)) {
  cat("Using feather to read data frame...\n")
  main_dt <- feather::read_feather(opt$embeddings)
} else {
  cat("Copying data to scratch and using fread on temp file...\n")
  embeddings_temp <- file.path(temp_dir, "glmnet_embeddings_temp.tsv")
  system(paste("cp", opt$embeddings, embeddings_temp))
  main_dt <- fread(embeddings_temp, header = T)
}



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

dev.off()

saveRDS(all_coef, file.path(temp_dir, "all_glmnet_coef.RDS"))
# write out the coefficients
# first filter for non-zero coefficients
# then join on the headers for the significant kmers
nonzero_coef <- all_coef %>% filter(!grepl("0,", coefficients))

nonzero_coef %>% write_tsv(file=coefficients_out, col_names = T, quote = "needed", na = "")
