library(data.table)
library(glmnet)
library(tidyverse)
library(optparse)

## parse command line arguments
option_list <- list(
  make_option(c("-i", "--input"), type="character", default=NULL, help="Input file"),
  make_option(c("-m", "--metadata"), type="character", default=NULL, help="Metadata file"),
  make_option(c("--metadata_column"), type="character", default=NULL, help="Metadata column to use for classification"),
  make_option(c("--min_num_per_category"), type="integer", default=100, help="Minimum number of samples per category"),
  make_option(c("-o", "--output_prefix"), type="character", default=NULL, help="Output file")
)
opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input) || is.null(opt$metadata) || is.null(opt$output_prefix)) {
  stop("Missing required arguments")
} else {
  input_file <- opt$input
  metadata_file <- opt$metadata
  metadata_column <- opt$metadata_column
  output_prefix <- opt$output_prefix
  min_num_per_category <- 100
}


# Load in strep A embeddings
main_dt <- fread(input_file, header=T)

all_metadata <- fread(metadata_file, header=T)

# if sample_name is not present in metadata, use the first column as sample_name
if (!"sample_name" %in% colnames(all_metadata)) {
  all_metadata <- all_metadata %>% rename(sample_name=colnames(all_metadata)[1])
}

# check if metadata column is present in metadata file
if (!is.null(metadata_column)) {
  if (!metadata_column %in% colnames(all_metadata)) {
    stop("Metadata column not found in metadata file")
  }
}

# get all metadata labels that are categorical and have less than 10 unique values
if (!is.null(metadata_column)) {
  metadata_labels <- metadata_column
} else {
  # for each possible metadata column, check if it is categorical and has more than 1 unique value
  # with greater than 20 observations each to avoid overfitting
  metadata_labels <- NULL
  for (i in 1:ncol(all_metadata)) {
    if (colnames(all_metadata)[[i]] == "sample_name") {
      next
    }
    if (is.character(all_metadata[[i]]) && n_distinct(all_metadata[[i]]) > 1) {
      if (all_metadata[[i]] %>% table() %>% as_tibble() %>% filter(n >= min_num_per_category) %>% nrow() > 1) {
        metadata_labels <- c(metadata_labels, colnames(all_metadata)[i])
      }
    }
  }
}

all_coef <- NULL

pdf(paste0(output_prefix, "_confusion_matrices.pdf"), width=7, height=7)
for (i in 1:length(metadata_labels)) {
  print(metadata_labels[i])
  metadata <- all_metadata %>%
    select(sample_name, metadata_labels[i]) %>%
    rename(class=metadata_labels[i])
  
  # merge metadata with embeddings
  dt <- merge(main_dt, metadata, by.x="sample_name", by.y="sample_name")
  
  # filter out samples with missing serotype
  dt <- dt[!is.na(class),]
  
  # filter out serotypes with less than 100 samples
  dt <- dt[, .SD[.N >= min_num_per_category], by=class]
  
  # break if there are less than 2 classes remaining
  print("Data Breakdown...")
  print("")
  print(dt[, .N, by=class])
  if (nrow(dt[, .N, by=class]) < 2) {
    print("Less than 2 classes remaining. Skipping...")
    next
  }
  
  # decide on the number of samples to keep in the training set
  # such that the number of samples is approximately equal for each class
  n_samples <- min(dt[, .N, by=class]$N)
  n_train <- floor(n_samples * 0.8)
  
  # separate into training and testing
  # keep equal numbers of each class in training and testing
  set.seed(123)
  # sample n_train from each class
  train_idx <- dt[, sample(.I, n_train), by=class]$V1
  train <- dt[train_idx]
  test <- dt[-train_idx]
  
  # prepare data for glmnet
  X <- as.matrix(train[, -c("sample_name", "class")])
  y <- train$class
  
  # fit glmnet
  fit <- cv.glmnet(X, y, family="multinomial", type.measure="class")
  
  # predict on test set
  X_test <- as.matrix(test[, -c("sample_name", "class")])
  y_test <- test$class
  y_pred <- predict(fit, X_test, s="lambda.min", type="class")
  
  # get the accuracy
  acc <- sum(y_test == y_pred) / length(y_test)
  
  # extract non-zero coefficients
  coef <- coef(fit, s="lambda.min")
  
  # transform multinomial coefficients into data frame
  coef_df <- map2_dfc(coef,names(coef), \(x,y) as.matrix(x) %>% as.data.frame() %>% setnames(y)) %>% 
    rownames_to_column("feature") %>% 
    unite("coefficients", -feature, sep=",") %>%
    mutate(coefficients=paste0("[", coefficients, "]")) %>%
    mutate(classes=paste0("[", paste(names(coef), collapse=","), "]")) %>%
    mutate(metadata_category=metadata_labels[i], accuracy = acc) %>%
    select(metadata_category, feature, accuracy, classes, coefficients)
  
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
    ggtitle(paste(metadata_labels[i], "| Accuracy:", round(acc, 2)))
  print(p)
  
  all_coef <- rbind(all_coef, coef_df)
}
dev.off()
write_tsv(all_coef, paste0(output_prefix, "_data_summary.tsv"))