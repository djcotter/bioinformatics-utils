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

## parse arguments --------
# define command line arguments
option_list <- list(
  make_option(c("-i", "--input"), "Input file", type="character"),
  make_option(c("-o", "--output"), "Output file", type="character"),
  make_option(c("-n", "--num_anchors"), "Number of anchors to select",
              type="integer", default = 100000),
  make_option(c("-e", "--effect_size"), "Effect size threshold",
              type="numeric", default=0.7),
  make_option(c("-l", "--lookup_table"), "Lookup table file", type="character"),
  make_option(c("--splash_bin"), "Path to SPLASH binary folder",
              type="character", default="/oak/stanford/groups/horence/dcotter1/splash-2.6.1/")
)

# parse command line arguments
opt <- parse_args(OptionParser(option_list = option_list))

# check that input file exists
if (!file.exists(opt$input) | !file.exists(opt$output)) {
  stop("Must provide input and output files")
}
if (!file.exists(opt$lookup_table)) {
  stop("Must provide path to lookup table")
}

# create a temporary directory to store intermediate files
temp_dir <- file.path(dirname(opt$output), "tmp/")
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
print(paste("Reading in data from", opt$input, "with effect size >=", opt$effect_size))
dt <- fread(cmd=load_cmd, header=TRUE, select = 1:18)

## run lookup table to filter out artifacts -----------
# write all anchors to a file
anchors_to_keep <- dt %>% select(anchor) %>% distinct() 
anchor_file <- file.path(temp_dir, "all_anchors.txt")
anchors_to_keep %>% write.table(anchor_file, row.names=FALSE, col.names=FALSE, quote=FALSE)

# run lookup table
out_lookup_stats <- file.path(temp_dir, "lookup_stats.txt")
lookup_cmd <- paste0(file.path(opt$splash_bin, "lookup_table"), 
                     " query --kmer_skip 1 --truncate_paths --stats_fmt with_stats ", 
                     anchor_file, " ", opt$lookup_table, " ", out_lookup_stats)
print("Running lookup table")
system(lookup_cmd)

# read in the lookup table stats
lookup_stats <- fread(out_lookup_stats, header=F, col.names=c("query", "stats"))
lookup_stats <- lookup_stats %>% mutate(anchor = anchors_to_keep$anchor)

# filter out anchors that have lookup table hits to artifacts
artifact_pattern <- ""
anchors_to_keep <- lookup_stats %>% filter(!grepl(artifact_pattern, stats)) %>% select(anchor)
print(anchors_to_keep)
## select the most important anchors -----------

