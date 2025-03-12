#!/bin/bash

#SBATCH --mem=26GB
#SBATCH --cpus-per-task=4
#SBATCH --time=300:00:00
#SBATCH --output=RC2_vOTU.o%j 

BASE_DIR="/scratch/sahlab/Luis/RC2_hecatomb/2023_04_07_RC2_Freeze4_HecatombOut_WholeRun_MMseqsFast/test_newPhageWorkflowCeliac"
PHAGES="${BASE_DIR}/06_phages_cluster/phageContigs.fasta"
CLUSTERING_DIR="${BASE_DIR}/06_phages_cluster"

source /ref/sahlab/software/miniforge3/bin/activate

# Filter genome pairs with at least 20 common 25-mers and 95% identity
vclust prefilter -i ${PHAGES} \
-o ${CLUSTERING_DIR}/vclust_fltr.txt \
--min-kmers 20 \
--min-ident 0.95 

# Align genome pairs filtered by the prefiltering command.
vclust align -i ${PHAGES} \
-o ${CLUSTERING_DIR}/vclust_ani.tsv \
--filter ${CLUSTERING_DIR}/vclust_fltr.txt 

# Cluster genomes based on ANI similarity measure, but connect them only
# if ANI ≥ 95% and query coverage ≥ 85%.
vclust cluster -i ${CLUSTERING_DIR}/vclust_ani.tsv \
-o ${CLUSTERING_DIR}/vclust_clusters.tsv \
--algorithm leiden \
--metric ani \
--ids ${CLUSTERING_DIR}/vclust_ani.ids.tsv \
--ani 0.95 \
--qcov 0.85 \
--rcov 0.85

#Extract lengths of the sequences
seqkit fx2tab --length --name ${PHAGES} > ${CLUSTERING_DIR}/seq_lengths.tsv

#Sort the clusters and lengths data files before joining
sort -k1,1 ${CLUSTERING_DIR}/vclust_clusters.tsv > ${CLUSTERING_DIR}/vclust_sorted_clusters.tsv
sort -k1,1 ${CLUSTERING_DIR}/seq_lengths.tsv > ${CLUSTERING_DIR}/sorted_seq_lengths.tsv

#Join the sorted files
join -1 1 -2 1 -t $'\t' ${CLUSTERING_DIR}/vclust_sorted_clusters.tsv ${CLUSTERING_DIR}/sorted_seq_lengths.tsv > ${CLUSTERING_DIR}/joined_Clusters_length.tsv

#Use awk to find the longest sequence in each cluster
awk -F'\t' '{
    if (!($2 in max_len) || max_len[$2] < $4) {
        max_len[$2] = $4;
        max_seq[$2] = $1;
    }
} END {
    for (c in max_seq)
        print max_seq[c];
}' ${CLUSTERING_DIR}/joined_Clusters_length.tsv > ${CLUSTERING_DIR}/vOTU_repSeqs.tsv


seqkit grep -f ${CLUSTERING_DIR}/vOTU_repSeqs.tsv ${PHAGES} > ${CLUSTERING_DIR}/vOTU_repSeqs.fasta