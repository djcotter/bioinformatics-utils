# select_important_anchors.R
# Daniel Cotter
# 2024-09-12

# This script selects the most important anchors from the output of SPLASH
# it prioritizes anchors that have a high effect size and also have a large 
# number of nonzero samples. It also filters out anchors that have lookup
# table hits to artifcats. The output is a list of anchor sequences that can
# be used in downstream analyses.

## import packages --------
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(Biostrings))

## parse arguments --------
# define command line arguments
option_list <- list(
  make_option(c("-i", "--input"), "Input file", type="character"),
  make_option(c("-o", "--output"), "Output file", type="character"),
  make_option(c("-n", "--num_anchors"), "Number of anchors to select",
              type="integer", default = 500000),
  make_option(c("--num_clusters"), "Number of clusters to create",
              type="integer", default = 1000),
  make_option(c("-e", "--effect_size"), "Effect size threshold",
              type="numeric", default=0.7),
  make_option(c("-l", "--lookup_table"), "Lookup table file", type="character"),
  make_option(c("--splash_bin"), "Path to SPLASH binary folder",
              type="character", default="/oak/stanford/groups/horence/dcotter1/splash-2.6.1/")
)

# parse command line arguments
opt <- parse_args(OptionParser(option_list = option_list))

# check that input file exists
if (!file.exists(opt$input) | is.null(opt$output)) {
  stop("Must provide input and output files")
}
if (!file.exists(opt$lookup_table)) {
  stop("Must provide path to lookup table")
}

# create a temporary directory to store intermediate files
temp_dir <- file.path(dirname(opt$output), "tmp/")
cat(paste("Creating temp directory:", temp_dir))
cat("\n")
system(paste("mkdir -p", temp_dir))

## load the data
# read in the headers of the input file to identify the effect_size_bin column
headers <- fread(opt$input, nrows = 1, header=T)
effect_size_bin_col <- grep("effect_size_bin", names(headers))

# load the input file using awk to filter out rows with effect size < 0.7
load_cmd <- paste0("cat ", opt$input, " | awk '{OFS=\"\t\"}{if ($", effect_size_bin_col, " >= ", opt$effect_size, ") print $0}'")
if (grepl(".gz$", opt$input)) {
  load_cmd = paste0("z", load_cmd)
}

# only select the first 18 columns
cat(paste("Reading in anchors from", opt$input, "with effect size >=", opt$effect_size))
cat("\n")
dt <- fread(cmd=load_cmd, header=TRUE, select = 1:18)

## run lookup table to filter out artifacts -----------
# write all anchors to a file
anchors_to_keep <- dt %>% select(anchor) %>% distinct() 
anchor_file <- file.path(temp_dir, "all_anchors.txt") %>% gsub("//", "/", .)
anchors_to_keep %>% write.table(anchor_file, row.names=FALSE, col.names=FALSE, quote=FALSE)

# write anchor to a fasta file as well (for lookuptable)
anchors_fasta_file <- file.path(temp_dir, "all_anchors.fasta") %>% gsub("//", "/", .)
anchors_dna <- DNAStringSet(anchors_to_keep$anchor)
names(anchors_dna) <- anchors_to_keep$anchor
writeXStringSet(anchors_dna, anchors_fasta_file)

# run lookup table
out_lookup_stats <- file.path(temp_dir, "lookup_stats.txt") %>% gsub("//", "/", .)
lookup_cmd <- paste0(file.path(opt$splash_bin, "lookup_table"), 
                     " query --kmer_skip 1 --truncate_paths --stats_fmt with_stats ", 
                     anchors_fasta_file, " ", opt$lookup_table, " ", out_lookup_stats)
if (file.exists(out_lookup_stats)) {
  cat("Lookup stats already in tmp directory. Delete tmp directory to force lookup table to run again...\n")
  cat("\nFiltering anchors for artifacts...")
} else {
  cat("Running lookup table...\n")
  system(lookup_cmd)
  cat("\nFinished Querying lookup table. Filtering anchors for artifacts...")
}

# read in the lookup table stats
lookup_stats <- fread(out_lookup_stats, header=F, col.names=c("query", "stats"))
lookup_stats <- lookup_stats %>% mutate(anchor = anchors_to_keep$anchor)

# filter out anchors that have lookup table hits to artifacts
artifact_pattern <- "plas|illum|syn|arp|RF|JUNK|Ral|purge|P,|Univec"
anchors_to_keep <- lookup_stats %>% filter(!grepl(artifact_pattern, query, ignore.case=T)) %>% select(anchor)

## select the most important anchors -----------
anchors_to_keep <- anchors_to_keep %>% 
  left_join(dt %>% select(anchor, number_nonzero_samples), by="anchor") %>%
  arrange(desc(number_nonzero_samples)) %>%
  head(opt$num_anchors)

## cluster the anchors by sequence similarity -----------
cat("Clustering anchors by sequence similarity...\n")
dist_matrix <- anchors_to_keep$anchor %>% DNAStringSet() %>% stringDist(method="levenshtein") %>% as.matrix()

# use spectral clustering to cluster the anchors
cat(paste("Creating", opt$num_clusters, "clusters of anchors.\n"))
suppressPackageStartupMessages(library(spectralClustering))
clust <- spectralClustering(dist_matrix, k=opt$num_clusters, method="discretize")

# add the cluster assignments to the anchors
anchors_to_keep <- anchors_to_keep %>% mutate(cluster=clust)
cat("Finished clustering anchors.\n")


## write the output -----------
cat(paste("Writing", nrow(anchors_to_keep), "anchors to", opt$output))
cat("\n")
anchors_to_keep %>% write.table(opt$output, row.names=FALSE, col.names=FALSE, quote=FALSE)