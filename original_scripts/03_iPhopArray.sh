#!/bin/bash

#SBATCH --mem=50GB
#SBATCH --time=100:00:00
#SBATCH --cpus-per-task=12
#SBATCH --output=iphopArray1_%A_%a.out
#SBATCH --array=1-394%40

source /ref/sahlab/software/miniforge3/bin/activate
conda activate /ref/sahlab/software/miniforge3/envs/iphop_env

BASE_DIR="/scratch/sahlab/Luis/RC2_hecatomb/2023_04_07_RC2_Freeze4_HecatombOut_WholeRun_MMseqsFast/test_newPhageWorkflowCeliac"
CLUSTERING_DIR="${BASE_DIR}/06_phages_Cluster"
PHAGES="${CLUSTERING_DIR}/vOTU_repSeqs.fasta"
CONTIG_INFO_DIR="${BASE_DIR}/07_phageGenomicInfo"


ID=$(sed -n ${SLURM_ARRAY_TASK_ID}p ${CONTIG_INFO_DIR}/iphop/lookup_magsSplit.txt) 
NAME=$(basename ${ID} .fasta); 

# Run iphop predict
iphop predict --fa_file ${CONTIG_INFO_DIR}/iphop/MAGs_split/${ID} \
--db_dir /ref/sahlab/data/viral_analysis_DBs/iphop_DBs/Aug_2023_pub_rw \
--out_dir ${CONTIG_INFO_DIR}/iphop/outputPerSplit/${NAME} \
--num_threads ${SLURM_CPUS_PER_TASK}