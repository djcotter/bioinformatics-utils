#!/bin/bash
#
#SBATCH --job-name=wrap_embeddings
#SBATCH -p horence,normal
#SBATCH -c 24
#SBATCH --mem=256G
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
SPLASH_RESULTS=
METADATA=
OUTPUT_PREFIX=
LOOKUP_TABLE=
TEMP_DIR=

## Define paramaters
NUM_CPUS=24
TRANSLATION_TABLE=11

## select anchors
if ! test -f ${OUTPUT_PREFIX}_selected_anchors.txt; then
  printf "Selecting anchors\n"
  Rscript --vanilla ${ANCHOR_SELECT_SCRIPT} --input ${SPLASH_RESULTS}/result.after_correction.scores.tsv --output ${OUTPUT_PREFIX}_selected_anchors.txt --lookup_table ${LOOKUP_TABLE} --temp_dir ${TEMP_DIR}
fi

## cluster_anchors
if ! test -f ${OUTPUT_PREFIX}_anchor_clusters.txt; then
  printf "Clustering anchors\n"
  python ${CLUSTER_SCRIPT} --input ${OUTPUT_PREFIX}_selected_anchors.txt --output ${OUTPUT_PREFIX}_anchor_clusters.txt
fi

## reorder_clusters
if ! test -f ${OUTPUT_PREFIX}_REORDERED_anchor_clusters.txt; then
  printf "Reordering clusters\n"
  Rscript --vanilla ${REORDER_CLUSTER_SCRIPT} --input_anchor_clusters ${OUTPUT_PREFIX}_anchor_clusters.txt \
    --splash_stats ${SPLASH_RESULTS}/result.after_correction.scores.tsv \
    --output ${OUTPUT_PREFIX}_REORDERED_anchor_clusters.txt \
    --temp_dir ${TEMP_DIR}
fi

if !${OUTPUT_PREFIX}_selected_anchors_top4000.txt; then
  # filter the clusters to only include the top 4000
  awk '$1 <= 4000' ${OUTPUT_PREFIX}_REORDERED_anchor_clusters.txt > ${OUTPUT_PREFIX}_REORDERED_anchor_clusters_top4000.txt
  cut -f2 ${OUTPUT_PREFIX}_REORDERED_anchor_clusters_top4000.txt > ${OUTPUT_PREFIX}_selected_anchors_top4000.txt
fi 

## prepare_sequences
if ! test -f ${OUTPUT_PREFIX}_sample_sequences.fasta; then
  printf "Preparing sequences\n"
  Rscript --vanilla ${PREPARE_SEQUENCES} --anchor_file ${OUTPUT_PREFIX}_selected_anchors_top4000.txt \
    --cluster_file ${OUTPUT_PREFIX}_REORDERED_anchor_clusters_top4000.txt \
    --id_mapping ${SPLASH_RESULTS}/sample_name_to_id.mapping.txt \
    --satc_files ${SPLASH_RESULTS}/result_satc \
    --output_prefix ${OUTPUT_PREFIX} --temp_dir ${TEMP_DIR} \
    --num_cores ${NUM_CPUS}
fi
  
## decompose_kmers
if ! test -f ${OUTPUT_PREFIX}_kmer_decomposition_k54_s54_unique_kmers.fasta; then
  printf "Decomposing kmers\n"
  python ${DECOMPOSE_KMERS} -k 54 ${OUTPUT_PREFIX}_sample_sequences.fasta ${OUTPUT_PREFIX}_kmer_decomposition
fi

## translate_kmers
if ! test -f ${OUTPUT_PREFIX}_kmer_decomposition_k54_s54_TRANSLATED_AA_unique_kmers.fasta; then
  printf "Translating kmers\n"
  python ${TRANSLATE_SCRIPT} -t ${TRANSLATION_TABLE} \
    ${OUTPUT_PREFIX}_kmer_decomposition_k54_s54_unique_kmers.fasta \
    ${OUTPUT_PREFIX}_kmer_decomposition_k54_s54_TRANSLATED_AA_unique_kmers.fasta
fi

## embed_kmers
if ! test -f ${OUTPUT_PREFIX}_unique_AA_esm_embeddings.tsv; then
  printf "Embedding kmers\n"
  sbatch --wait ${EMBED_SBATCH} ${OUTPUT_PREFIX}_kmer_decomposition_k54_s54_TRANSLATED_AA_unique_kmers.fasta ${TEMP_DIR} ${OUTPUT_PREFIX}
fi

## format embeddings for glmnet -- variance
if ! test -f ${OUTPUT_PREFIX}_top_variance_embeddings.tsv; then
  printf "Formatting embeddings for glmnet\n"
  Rscript --vanilla ${FORMAT_EMBEDDINGS_VARIANCE} --embeddings ${OUTPUT_PREFIX}_unique_AA_esm_embeddings.tsv \
    --ordering ${OUTPUT_PREFIX}_kmer_decomposition_k54_s54_kmer_ordering.tsv \
    --output_prefix ${OUTPUT_PREFIX} \
    --temp_dir ${TEMP_DIR}
fi

## format embeddings for glmnet -- PCs
# Rscript --vanilla ${FORMAT_EMBEDDINGS_PCS} --embeddings ${OUTPUT_PREFIX}_unique_AA_esm_embeddings.tsv \
#   --ordering ${OUTPUT_PREFIX}_kmer_decomposition_k54_s54_kmer_ordering.tsv \
#   --output_prefix ${OUTPUT_PREFIX} \
#   --temp_dir ${TEMP_DIR}
  
## submit glmnet job -- Variance
if ! test -f ${OUTPUT_PREFIX}_topVar_glmnet_nonzero_coefficients.tsv; then
  printf "Submitting batch job for top variance glmnet\n"
  sbatch -p horence,normal -c 16 --mem=256G --time=4:00:00 Rscript --vanilla ${GLMNET_SCRIPT} \
    --embeddings ${OUTPUT_PREFIX}_top_variance_embeddings.tsv \
    --metadata ${METADATA} \
    --output_prefix ${OUTPUT_PREFIX}_topVar_glmnet \
    --even_classes \
    --temp_dir ${TEMP_DIR}
fi
  
# ## submit glmnet job -- PCs
# sbatch -p horence,normal -c 16 --mem=256G --time=4:00:00 Rscript --vanilla ${GLMNET_SCRIPT} \
#   --embeddings ${OUTPUT_PREFIX}_top_pcs_embeddings.tsv \
#   --metadata ${METADATA} \
#   --output_prefix ${OUTPUT_PREFIX}_topPCs_glmnet \
#   --even_classes \
#   --temp_dir ${TEMP_DIR}