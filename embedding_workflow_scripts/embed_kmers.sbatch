#!/bin/bash
#
#SBATCH --job-name=embed_kmers
#SBATCH -p gpu,owners
#SBATCH -c 8
#SBATCH -G 1
#SBATCH -C GPU_GEN:AMP|GPU_GEN:VLT|GPU_GEN:TUR
#SBATCH --mem=32G
#SBATCH --time=2:00:00

ml python/3.9.0
source /oak/stanford/groups/horence/dcotter1/envs/esm_venv/bin/activate

## declare path to torch directory
TORCH_DIR=/scratch/users/dcotter1/torch_home

## get the file paths from the command line 
INPUT_FASTA=$1
TEMP_DIR=$2
OUTPUT_PREFIX=$3

## calculate the embeddings
export TORCH_HOME=${TORCH_DIR}
esm-extract esm2_t33_650M_UR50D ${INPUT_FASTA} ${TEMP_DIR}/esm_embeddings --include mean per_tok

EXTRACT_EMBEDDINGS_SCRIPT=/oak/stanford/groups/horence/dcotter1/bioinformatics-utils/extract_embeddings_from_torch_objects.py
OUTPUT_EMBEDDINGS=${OUTPUT_PREFIX}_unique_AA_esm_embeddings.tsv

## extract the embeddings
python ${EXTRACT_EMBEDDINGS_SCRIPT} ${TEMP_DIR}/esm_embeddings ${OUTPUT_EMBEDDINGS}




