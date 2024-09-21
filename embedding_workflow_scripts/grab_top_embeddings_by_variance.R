# grab_top_embeddings_by_variance.R
# Daniel Cotter
# 2024-09-18

# This script takes in one embeddings file and the ordering file and outputs a new
# embeddings tsv with samples as rows and the top 10 embeddings by variance per cluster as columns.


## import packages --------
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(resample))


## parse arguments --------
# define command line arguments
# define command line arguments
option_list <- list(
  make_option(c("-e", "--embeddings"), help="Embeddings file", type="character"),
  make_option(c("-o", "--ordering"), help="Ordering file", type="character"),
  make_option(c("-n", "--num_to_keep"), help="Number of components to keep per cluster", 
              type="integer", default = 10),
  make_option(c("-p", "--output_prefix"), help="Output prefix.", type="character"),
  make_option(c("--temp_dir"), help="Temporary directory to store intermediate files", 
              type="character")
)



# parse command line arguments
opt <- parse_args(OptionParser(option_list = option_list))

# check that user specified all files
if (!file.exists(opt$embeddings) | !file.exists(opt$ordering) | is.null(opt$output_prefix)) {
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

## print a summary of the arguments
cat("\n####################\n")
cat("Running grab_top_embeddings_by_variance.R with the following arguments:\n")
cat("Ordering file: ", opt$ordering, "\n")
cat("Embeddings file: ", opt$embeddings, "\n")
cat("Output file: ", opt$output, "\n")
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

# write out the cluster to kmer mapping
cluster_to_kmer_mapping_file <- paste0(opt$output_prefix, "_cluster_to_kmer_mapping.tsv")
cat("Writing cluster to kmer mapping to ", cluster_to_kmer_mapping_file, "\n")
write_tsv(cluster_to_kmer_mapping, cluster_to_kmer_mapping_file)

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
embeddings_topVar_file <- paste0(opt$output_prefix, "_top_variance_embeddings.tsv")
cat("Writing top variance embeddings to ", embeddings_topVar_file, "\n")
write_tsv(main_dt, embeddings_topVar_file)