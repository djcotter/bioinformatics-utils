# reorder_anchor_clusters.R
# Daniel Cotter
# 2024-09-12

# This script takes in an anchor cluster file and a splash stats file and smartly 
# reorders the anchor clusters (both the clusters and the anchors within the clusters)
# based on the effect size and the number of nonzero samples. The output is a new 
# anchor cluster file that can be used in downstream analyses.

## import packages --------
suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(optparse))

## parse arguments --------
# define command line arguments
option_list <- list(
  make_option(c("-i", "--input_anchor_clusters"), "Input anchor clusters", type="character"),
  make_option(c("-s", "--splash_stats"), "SPLASH stats file", type="character"),
  make_option(c("-o", "--output"), "Output file", type="character"),
  make_option(c("--temp_dir"), "Temporary directory to store intermediate files", 
              type="character")
)

# parse command line arguments
opt <- parse_args(OptionParser(option_list = option_list))

# check that input file exists
if (!file.exists(opt$input_anchor_clusters) | !file.exists(opt$splash_stats) | is.null(opt$output)) {
  stop("Must provide input and output files")
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

# cat to screen the input files, output files, temp dir and paramaters
cat("\n###################################################################\n")
cat("Running select_important_anchors.R with the following parameters:\n")
cat(paste("Input clusters:", opt$input_anchor_clusters))
cat("\n")
cat(paste("SPLASH stats file:", opt$splash_stats))
cat("\n")
cat(paste("Output file:", opt$output))
cat("\n")
cat(paste("Number of anchors to select:", opt$num_anchors))
cat("\n")
cat(paste("Temporary directory:", temp_dir))
cat("\n")
cat("###################################################################\n\n")

# read in the anchor cluster file
cat("Reading in the anchor clusters file\n")
anchor_clusters <- fread(opt$input_anchor_clusters, 
                    header = F, col.names = c("cluster_id", "anchor"))

# write out the anchors column to a temporary file
anchor_file <- file.path(temp_dir, "anchors_to_select.txt")
write_tsv(anchor_clusters %>% select(anchor), anchor_file, col_names = F, quote="none")

# read in the splash stats file (grepping for the anchors in the anchor file)
cat("Reading in the splash stats file\n")
read_cmd = paste0("grep -f ", anchor_file, " ", opt$splash_stats)
splash_stats <- fread(cmd=read_cmd, header = T, select = 1:18)

# join the splash stats file with the anchor clusters file
cat("Joining the splash stats file with the anchor clusters file\n")
anchor_clusters_with_stats <- anchor_clusters %>%
  left_join(splash_stats, by = "anchor")

# reorder the cluster ids based on the mean effect size and number of nonzero samples
cat("Reordering the cluster ids based on the mean effect size and number of nonzero samples\n")
anchor_clusters_with_stats <- anchor_clusters_with_stats %>%
  group_by(cluster_id) %>%
  mutate(mean_effect_size = mean(effect_size),
         num_nonzero_samples = sum(num_nonzero_samples)) %>%
  arrange(desc(mean_effect_size), desc(num_nonzero_samples)) %>%
  mutate(new_cluster_id = cur_group_id()) %>% 
  ungroup() %>% group_by(new_cluster_id) %>%
  arrange(new_cluster_id, desc(effect_size_bin * num_nonzero_samples)) %>%
  ungroup() %>% select(new_cluster_id, anchor)

# write out the new anchor clusters file
cat("Writing out the new anchor clusters file\n")
write_tsv(anchor_clusters_with_stats, opt$output, col_names = F, quote="none")

cat("Done!\n")


