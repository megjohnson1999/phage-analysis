#!/bin/bash

#SBATCH --mem=10GB
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=4
#SBATCH --output=phactsArray_%A_%a.out
#SBATCH --array=1-20000%40


source /ref/sahlab/software/miniforge3/bin/activate

BASE_DIR="/scratch/sahlab/Luis/RC2_hecatomb/2023_04_07_RC2_Freeze4_HecatombOut_WholeRun_MMseqsFast/test_newPhageWorkflowCeliac"
CLUSTERING_DIR="${BASE_DIR}/06_phages_Cluster"
PHAGES="${CLUSTERING_DIR}/vOTU_repSeqs.fasta"
CONTIG_INFO_DIR="${BASE_DIR}/07_phageGenomicInfo"

proteome=$(sed -n ${SLURM_ARRAY_TASK_ID}p ${CONTIG_INFO_DIR}/phacts/lookup_04_SplitedByMAG_PredictedORFsFromMags.txt)

NAME=$(basename ${proteome} .faa); 

python /home/luisalberto/Softwares/PHACTS/phacts.py ${CONTIG_INFO_DIR}/phacts/MAGs_split_predictedORFs/${proteome} -o ${CONTIG_INFO_DIR}/phacts/predictionPerProteome/${NAME}