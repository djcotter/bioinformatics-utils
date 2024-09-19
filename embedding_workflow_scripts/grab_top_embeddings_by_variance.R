# grab_top_embeddings_by_variance.R
# Daniel Cotter
# 2024-09-18

# This script takes in one embeddings file and a paramater with the number of embedding dimensions, and
# outputs a trimmed embeddings file with only the top embeddings per anchor by variance.


## import packages --------
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(resample))


## parse arguments --------
# define command line arguments
option_list <- list(
  make_option(c("-e", "--embeddings"), help="Embeddings file", type="character"),
  make_option(c("-n", "--embeddings_dim"), help="Number of embeddings dimensions", type="integer"),
  make_option(c("--num_to_keep"), help="Number of embeddings to keep per anchor", type="integer", default =10),
  make_option(c("-p", "--output"), help="Output file", type="character"),
  make_option(c("--temp_dir"), help="Temporary directory to store intermediate files", 
              type="character")
)


# parse command line arguments
opt <- parse_args(OptionParser(option_list = option_list))

# check that user specified all files and parameters
if (!file.exists(opt$embeddings) | is.null(opt$output_file) | is.null(opt$embeddings_dim)) {
  stop("Must provide embeddings file, output file, and number of embeddings")
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

## print a summary of the arguments
cat("\n####################\n")
cat("Running grab_top_embeddings_by_variance.R with the following arguments:\n")
cat("Embeddings file: ", opt$embeddings, "\n")
cat("Output file: ", opt$output, "\n")
cat("Num embeddings dimensions: ", opt$embeddings_dim, "\n")
cat("Num embeddings to keep per anchor: ", opt$num_to_keep, "\n")
cat("Temporary directory: ", temp_dir, "\n")
cat("####################\n\n")

## load data --------
# load embeddings
cat("\nLoading embeddings...\n")
# copy the embeddings file to the temp directory to speed up I/O
embeddings_temp <- file.path(temp_dir, "embeddings_temp.tsv")
if (!file.exists(embeddings_temp)) {
  system(paste("cp", opt$embeddings, embeddings_temp))
}
embeddings <- fread(embeddings_temp, header = F)
# name the files as cluster_\d_embedding_[1-opt$embeddings_dim]
# determine the number of anchors by dividing the number of embeddings
# columns by the number of embeddings dimensions
num_anchors <- (ncol(embeddings)-1) / opt$embeddings_dim
if (num_anchors %% 1 != 0) {
  stop("Number of embeddings columns is not divisible by the number of embeddings dimensions")
}
new_col_names <- paste0("cluster_", rep(1:num_anchors, each=opt$embeddings_dim), "-embedding_", rep(1:opt$embeddings_dim, times=num_anchors))
new_col_names <- c("sample_name", new_col_names)

main_dt <- embeddings
colnames(main_dt) <- new_col_names

# calculate the variance of each embedding column
cat("Grabbing the top", opt$num_to_keep, "embeddings by variance per cluster...\n")
col_variances <- resample::colVars(main_dt %>% select(starts_with("cluster"))) %>% 
  enframe() %>% separate(name, into=c("cluster", "embedding"), sep="-")

# grab the top embeddings by variance
top_var_cols <- col_variances %>% 
  group_by(cluster) %>% 
  slice_max(value, n=opt$num_to_keep, with_ties = FALSE) %>%
  unite(name, cluster, embedding, sep="-") %>%
  pull(name)

# select the top embeddings by variance
main_dt <- main_dt %>% select("sample_name", all_of(top_var_cols))

# write the output file
cat("Writing the output file...\n")
write_tsv(main_dt, opt$output, col_names = T, quote = F)

cat("Done!\n")