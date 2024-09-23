# grab_all_embeddings_from_clusters.R
# Daniel Cotter
# 2024-09-18

# This script takes in one embeddings file and the ordering file and outputs a new
# embeddings tsv with samples as rows and all embeddings for the first N clusters defined.


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
  make_option(c("-e", "--embeddings"), help="Embeddings file", type="character"),
  make_option(c("-o", "--ordering"), help="Ordering file", type="character"),
  make_option(c("-n", "--num_clusters_to_keep"), help="Number of components to keep per cluster", 
              type="integer", default = 200),
  make_option(c("-p", "--output_prefix"), help="Output prefix.", type="character"),
  make_option(c("--temp_dir"), help="Temporary directory to store intermediate files", 
              type="character"),
  make_option(c("--num_threads"), help="Number of threads to use for parallel operations",
              type="integer", default=1)
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

# define future plans
setDTthreads(opt$num_threads)
plan(multicore, workers=opt$num_threads)
options(future.globals.maxSize=1000*1024^2)

## print a summary of the arguments
cat("\n####################\n")
cat("Running grab_top_embeddings_by_variance.R with the following arguments:\n")
cat("Ordering file: ", opt$ordering, "\n")
cat("Embeddings file: ", opt$embeddings, "\n")
cat("Output prefix: ", opt$output_prefix, "\n")
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


cluster_to_kmer_mapping <- ordering %>% select(cluster, kmer) %>% distinct()
cluster_to_kmer_mapping <- cluster_to_kmer_mapping %>% 
  arrange(cluster) %>% group_by(cluster) %>% 
  summarise(kmers = str_c(kmer, collapse=",")) %>% ungroup()

# write out the cluster to kmer mapping
cluster_to_kmer_mapping_file <- paste0(opt$output_prefix, "_cluster_to_kmer_mapping.tsv")
cat("Writing cluster to kmer mapping to ", cluster_to_kmer_mapping_file, "\n")
write_tsv(cluster_to_kmer_mapping, cluster_to_kmer_mapping_file)

cat("Formatting the embeddings for downstream use...\n")

# define a function to grab
grab_all_embedding_columns <- function(in_file) {
  cluster_num = str_extract(in_file, "cluster_(\\d+).csv", group=1) %>% as.integer()
  temp_dt <- fread(in_file, header=T, nThread = 1) %>% 
    select(sample_name, starts_with("embedding")) # first filter for only one cluster
  colnames(temp_dt) <- ifelse(grepl("embedding", colnames(temp_dt)),
                              yes=paste0("cluster_", cluster_num, "_", colnames(temp_dt)),
                              no = colnames(temp_dt))
  temp_dt <- temp_dt %>% arrange(sample_name)
  if (cluster_num!=0) {
    temp_dt <- temp_dt %>% select(-sample_name)
  }
  return(temp_dt)
}


cluster_dt <- future_map_dfc(cluster_files[1:opt$num_clusters_to_keep],
                             \(x) grab_all_embedding_columns(x),
                             .progress = T)
cluster_dt <- cluster_dt %>% relocate(sample_name)

# write out the embeddings matrix to a temp file
embeddings_cluster_file <- paste0(opt$output_prefix, "_embeddings_for_", opt$num_clusters_to_keep, "_clusters.tsv")
cat("Writing embeddings for first ", opt$num_clusters_to_keep, " clusters to ", embeddings_cluster_file, "\n")
write_tsv(cluster_dt, embeddings_cluster_file, col_names = T, quote="needed")
embeddings_cluster_feather <- file.path(temp_dir, 
                                    basename(paste0(opt$output_prefix, 
                                                    "_embeddings_for_", 
                                                    opt$num_clusters_to_keep,
                                                    "_clusters.feather")))
feather::write_feather(cluster_dt, embeddings_cluster_feather)
