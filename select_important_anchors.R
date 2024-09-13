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
  make_option(c("-o", "--output_prefix"), "Output prefix.", type="character"),
  make_option(c("-n", "--num_anchors"), "Number of anchors to select",
              type="integer", default = 150000),
  make_option(c("--num_clusters"), "Number of clusters to create",
              type="integer", default = 3000),
  make_option(c("-e", "--effect_size"), "Effect size threshold",
              type="numeric", default=0.7),
  make_option(c("-l", "--lookup_table"), "Lookup table file", type="character"),
  make_option(c("--splash_bin"), "Path to SPLASH binary folder",
              type="character", default="/oak/stanford/groups/horence/dcotter1/splash-2.6.1/"),
  make_option(c("--temp_dir"), "Temporary directory to store intermediate files", 
              type="character")
)

# parse command line arguments
opt <- parse_args(OptionParser(option_list = option_list))

# check that input file exists
if (!file.exists(opt$input) | is.null(opt$output_prefix)) {
  stop("Must provide input and output files")
}
if (!file.exists(opt$lookup_table)) {
  stop("Must provide path to lookup table")
}

# define output files 
anchors_only_out = paste0(opt$output_prefix, "_anchor_list.txt")
anchor_clusters_out = paste0(opt$output_prefix, "_anchor_clusters.tsv")

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

# cat to screen the input files, output files, temp dir and paramaters
cat("\n###################################################################\n")
cat("Running select_important_anchors.R with the following parameters:\n")
cat(paste("Input file:", opt$input))
cat("\n")
cat(paste("Output anchors:", anchors_only_out))
cat("\n")
cat(paste("Output anchor clusters:", anchor_clusters_out))
cat("\n")
cat(paste("Number of anchors to select:", opt$num_anchors))
cat("\n")
cat(paste("Number of clusters to create:", opt$num_clusters))
cat("\n")
cat(paste("Effect size threshold:", opt$effect_size))
cat("\n")
cat(paste("Lookup table file:", opt$lookup_table))
cat("\n")
cat(paste("SPLASH binary folder:", opt$splash_bin))
cat("\n")
cat(paste("Temporary directory:", temp_dir))
cat("\n")
cat("###################################################################\n\n")


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
cat("\n\n")
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
  cat("Filtering anchors for artifacts...\n")
} else {
  cat("Running lookup table...\n")
  system(lookup_cmd)
  cat("Finished Querying lookup table. Filtering anchors for artifacts...\n\n")
}

# read in the lookup table stats
lookup_stats <- fread(out_lookup_stats, header=F, col.names=c("query", "stats"))
lookup_stats <- lookup_stats %>% mutate(anchor = anchors_to_keep$anchor)

# filter out anchors that have lookup table hits to artifacts
artifact_pattern <- "plas|illum|syn|arp|RF|JUNK|Ral|purge|P,|Univec"
anchors_to_keep <- lookup_stats %>% filter(!grepl(artifact_pattern, query, ignore.case=T)) %>% select(anchor)

cat(paste0("Finished filtering. Kept ", nrow(anchors_to_keep), " anchors out of ", nrow(lookup_stats), " total anchors.\n"))
cat(paste0("Keeping the top ", opt$num_anchors, " by number_nonzero_samples for further analysis...\n\n"))

## select the most important anchors -----------
anchors_to_keep <- anchors_to_keep %>% 
  left_join(dt %>% select(anchor, number_nonzero_samples), by="anchor") %>%
  arrange(desc(number_nonzero_samples)) %>%
  head(opt$num_anchors)

# cluster the anchors using kmer distances
cat(paste("Clustering anchors into", opt$num_clusters, "clusters by kmer similarity...\n"))
# time
start_time <- Sys.time()
suppressPackageStartupMessages(library(ape))
suppressPackageStartupMessages(library(kmer))
suppressPackageStartupMessages(library(dendextend))

dna <- anchors_to_keep$anchor %>% DNAStringSet(.)
names(dna) <- anchors_to_keep$anchor
dna <- as.DNAbin(dna)
dna.tree <- kmer::cluster(dna)
anchor_clusters <- dendextend::cutree(dna.tree, 
                                      k=opt$num_clusters) %>% 
  enframe() %>% 
  dplyr::rename(anchor=name, cluster_id=value)
cat(paste("Finished clustering anchors in", round(Sys.time() - start_time, 2), "seconds.\n\n"))

## format anchor clusters ----------
cat("Arranging clusters by effect size...\n")
anchor_clusters <- anchor_clusters %>% left_join(dt, by="anchor") %>% 
  arrange(cluster_id, desc(effect_size_bin)) %>% select(cluster_id, anchor)
anchors_only <- anchor_clusters %>% select(anchor)

## write the output -----------
cat(paste("Writing", nrow(anchors_to_keep), "anchors to", anchors_only_out))
cat("\n")
anchors_only %>% write.table(anchors_only_out, row.names=FALSE, col.names=FALSE, quote=FALSE)

cat(paste("Writing", opt$num_clusters, "clusters and their anchors to", anchor_clusters_out))
cat("\n")
anchor_clusters %>% write.table(anchor_clusters_out, row.names=FALSE, col.names=FALSE, quote=FALSE)

