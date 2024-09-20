#!/bin/bash
#
#SBATCH --job-name=wrap_embeddings
#SBATCH -p horence,normal
#SBATCH -c 16
#SBATCH --mem=128G
#SBATCH --time=16:00:00

## load modules
ml python/3.9.0
ml R/4.3.2
source /home/users/dcotter1/my_env/bin/activate

## Define scripts
ANCHOR_SELECT_SCRIPT=/oak/stanford/groups/horence/dcotter1/bioinformatics-utils/embedding_workflow_scripts/select_anchors.R
CLUSTER_SCRIPT=/oak/stanford/groups/horence/dcotter1/bioinformatics-utils/embedding_workflow_scripts/deduplication_clustering.py
REORDER_CLUSTER_SCRIPT=/oak/stanford/groups/horence/dcotter1/bioinformatics-utils/embedding_workflow_scripts/reorder_anchor_clusters.R
PREPARE_SEQUENCES=/oak/stanford/groups/horence/dcotter1/bioinformatics-utils/embedding_workflow_scripts/process_splash_data_for_embedding.R
DECOMPOSE_KMERS=/oak/stanford/groups/horence/dcotter1/bioinformatics-utils/kmer_decomposition.py
TRANSLATE_SCRIPT=/oak/stanford/groups/horence/dcotter1/bioinformatics-utils/translate_fasta.py
EMBED_SBATCH=/oak/stanford/groups/horence/dcotter1/bioinformatics-utils/embedding_workflow_scripts/embed_kmers.sbatch
FORMAT_EMBEDDINGS_VARIANCE=/oak/stanford/groups/horence/dcotter1/bioinformatics-utils/embedding_workflow_scripts/grab_top_embeddings_by_variance.R
GLMNET_SCRIPT=/oak/stanford/groups/horence/dcotter1/bioinformatics-utils/embedding_workflow_scripts/embeddings_glmnet.R

## Define data paths
SPLASH_DIR=
METADATA=
OUTPUT_PREFIX=
LOOKUP_TABLE=
TEMP_DIR=

## Define paramaters
NUM_CPUS=16
TRANSLATION_TABLE=11

## select anchors
printf "Selecting anchors\n"
Rscript Rscript --vanilla ${SCRIPT} --input ${SPLASH_RESULTS}/result.after_correction.scores.tsv --output ${OUTPUT_PREFIX}_selected_anchors.txt --lookup_table ${LOOKUP_TABLE} --temp_dir ${TEMP_DIR}

## cluster_anchors
printf "Clustering anchors\n"
python ${CLUSTER_SCRIPT} --input ${OUTPUT_PREFIX}_selected_anchors.txt --output ${OUTPUT_PREFIX}_anchor_clusters.txt

## reorder_clusters
printf "Reordering clusters\n"
Rscript --vanilla ${REORDER_CLUSTER_SCRIPT} --input_anchor_clusters ${OUTPUT_PREFIX}_anchor_clusters.txt \
  --splash_stats ${SPLASH_RESULTS}/result.after_correction.scores.tsv \
  --output ${OUTPUT_PREFIX}_REORDERED_anchor_clusters.txt \
  --temp_dir ${TEMP_DIR}
  
## prepare_sequences
printf "Preparing sequences\n"
Rscript --vanilla ${PREPARE_SEQUENCES} --anchor_file ${OUTPUT_PREFIX}_selected_anchors.txt \
  --cluster_file ${OUTPUT_PREFIX}_REORDERED_anchor_clusters.txt \
  --id_mapping ${SPLASH_RESULTS}/sample_name_to_id.mapping.txt \
  --satc_files ${SPLASH_RESULTS}/result_satc \
  --output_prefix ${OUTPUT_PREFIX} --temp_dir ${TEMP_DIR} \
  --num_cores ${NUM_CPUS}
  
## decompose_kmers
printf "Decomposing kmers\n"
python ${DECOMPOSE_KMERS} -k 54 ${OUTPUT_PREFIX}_sample_sequences.fasta ${OUTPUT_PREFIX}_kmer_decomposition

## translate_kmers
printf "Translating kmers\n"
python ${TRANSLATE_SCRIPT} -t ${TRANSLATION_TABLE} \
  ${OUTPUT_PREFIX}_kmer_decomposition_k54_s54_unique_kmers.fasta \
  ${OUTPUT_PREFIX}_kmer_decomposition_k54_s54_TRANSLATED_AA_unique_kmers.fasta
  
## embed_kmers
printf "Embedding kmers\n"
sbatch --wait ${EMBED_SBATCH} ${OUTPUT_PREFIX}_kmer_decomposition_k54_s54_TRANSLATED_AA_unique_kmers.fasta ${TEMP_DIR} ${OUTPUT_PREFIX}

## format embeddings for glmnet -- variance
Rscript --vanilla ${FORMAT_EMBEDDINGS_VARIANCE} --embeddings ${OUTPUT_PREFIX}_unique_AA_esm_embeddings.tsv \
  --ordering ${OUTPUT_PREFIX}_kmer_decomposition_k54_s54_kmer_ordering.tsv \
  --output_prefix ${OUTPUT_PREFIX} \
  --temp_dir ${TEMP_DIR}
  
## format embeddings for glmnet -- PCs
# Rscript --vanilla ${FORMAT_EMBEDDINGS_PCS} --embeddings ${OUTPUT_PREFIX}_unique_AA_esm_embeddings.tsv \
#   --ordering ${OUTPUT_PREFIX}_kmer_decomposition_k54_s54_kmer_ordering.tsv \
#   --output_prefix ${OUTPUT_PREFIX} \
#   --temp_dir ${TEMP_DIR}
  
## submit glmnet job -- Variance
sbatch -p horence,normal -c 16 --mem=256G --time=4:00:00 Rscript --vanilla ${GLMNET_SCRIPT} \
  --embeddings ${OUTPUT_PREFIX}_top_variance_embeddings.tsv \
  --metadata ${METADATA} \
  --output_prefix ${OUTPUT_PREFIX}_topVar_glmnet \
  --even_classes \
  --temp_dir ${TEMP_DIR}
  
## submit glmnet job -- PCs
sbatch -p horence,normal -c 16 --mem=256G --time=4:00:00 Rscript --vanilla ${GLMNET_SCRIPT} \
  --embeddings ${OUTPUT_PREFIX}_top_pcs_embeddings.tsv \
  --metadata ${METADATA} \
  --output_prefix ${OUTPUT_PREFIX}_topPCs_glmnet \
  --even_classes \
  --temp_dir ${TEMP_DIR}