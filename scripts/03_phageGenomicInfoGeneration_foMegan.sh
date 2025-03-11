#!/bin/bash

#SBATCH --mem=700GB
#SBATCH --cpus-per-task=24
#SBATCH --time=300:00:00
#SBATCH --output=abx_genomeInfo_gen.o%j 

source /ref/sahlab/software/miniforge3/bin/activate

BASE_DIR="/scratch/sahlab/Luis/RC2_hecatomb/2023_04_07_RC2_Freeze4_HecatombOut_WholeRun_MMseqsFast/test_newPhageWorkflowCeliac"
CLUSTERING_DIR="${BASE_DIR}/06_phages_cluster"
PHAGES="${CLUSTERING_DIR}/vOTU_repSeqs.fasta"
CONTIG_INFO_DIR="${BASE_DIR}/07_phageGenomicInfo"


mkdir -p ${CONTIG_INFO_DIR}

###General workflow for taxonomy and lifestyle preiction: PHABOX###
conda activate /ref/sahlab/software/miniforge3/envs/phabox2
phabox2 --task end_to_end --dbdir /ref/sahlab/data/viral_analysis_DBs/phabox_db_v2 \
        --outpth  ${CONTIG_INFO_DIR}/Phabox_OUT \
        --contigs ${PHAGES} \
        --len 1000 \
        --threads ${SLURM_CPUS_PER_TASK}

echo "phabox2 done"

###Taxonomy###

#VC3
conda activate /ref/sahlab/software/miniforge3/envs/vcontct3

#vcontact3 prepare_databases --get-version "223" --set-location /ref/sahlab/data/viral_analysis_DBs/vcontact3_DB223
vcontact3 run --nucleotide ${PHAGES} \
--output ${CONTIG_INFO_DIR}/vc3_OUT \
--db-domain "prokaryotes" \
--db-version 223 \
-t ${SLURM_CPUS_PER_TASK} \
--db-path /ref/sahlab/data/viral_analysis_DBs/vcontact3_DB223

echo "vc3 done"

#mmseqs NR
conda activate /ref/sahlab/software/miniforge3/envs/mmseqs2_v15-6f452
DB=/ref/sahlab/data/nr/nr_mmseqs15_DB
TMP=/scratch/sahlab/Luis/temporal/

mkdir -p ${TMP}
mkdir -p ${CONTIG_INFO_DIR}/mmseqs

mmseqs easy-taxonomy ${PHAGES} ${DB} ${CONTIG_INFO_DIR}/mmseqs/mmseqs_vOTU_repSeqs ${TMP} \
--min-length 30 \
-e 1e-15 \
--search-type 2 \
-s 4.0 \
--shuffle 0 \
--lca-mode 2  \
-a \
--tax-lineage 2 \
--format-output "query,target,evalue,pident,fident,nident,mismatch,qcov,tcov,qstart,qend,qlen,tstart,tend,tlen,alnlen,bits,qheader,theader,taxid,taxname,taxlineage" \
--threads ${SLURM_CPUS_PER_TASK} \
--split-mode 0 \
--orf-filter 1

# echo "mmseqs done"
#crassUS for crassphage taxa ran interactively and manually#

###lifestyle###

#PHACTS
source /ref/sahlab/software/miniforge3/bin/activate

mkdir -p ${CONTIG_INFO_DIR}/phacts/
mkdir -p ${CONTIG_INFO_DIR}/phacts/MAGs_split
mkdir -p ${CONTIG_INFO_DIR}/phacts/MAGs_split_predictedORFs
mkdir -p ${CONTIG_INFO_DIR}/phacts/predictionPerProteome
mkdir -p ${CONTIG_INFO_DIR}/phacts/logs

#Split multifasta file into individual fasta files
seqkit split -s 1 -O ${CONTIG_INFO_DIR}/phacts/MAGs_split ${PHAGES}
for fasta in "${CONTIG_INFO_DIR}/phacts/MAGs_split"/*.fasta; do 
	header=$(grep '^>' "$fasta" | head -1); new_name=$(echo "$header" | sed 's/>//;s/ .*//'); 
	mv "$fasta" "${CONTIG_INFO_DIR}/phacts/MAGs_split/${new_name}.fasta"; 
done

#Predict proteomes for each fasta file
for fasta in ${CONTIG_INFO_DIR}/phacts/MAGs_split/*.fasta; do
    base=$(basename $fasta .fasta)
    prodigal -i $fasta -a ${CONTIG_INFO_DIR}/phacts/MAGs_split_predictedORFs/${base}.faa -p meta
done

#Generate lookup file. If the number of proteomes is higher that 20000, we need to split the files because 20000 is the array limiti per job 
find ${CONTIG_INFO_DIR}/phacts/MAGs_split_predictedORFs/ -name "*.faa" -exec basename {} \; > ${CONTIG_INFO_DIR}/phacts/lookup_04_SplitedByMAG_PredictedORFsFromMags.txt

split -l 20000 -d ${CONTIG_INFO_DIR}/phacts/lookup_04_SplitedByMAG_PredictedORFsFromMags.txt "${CONTIG_INFO_DIR}/phacts/lookup_04_SplitedByMAG_PredictedORFsFromMags_part" --additional-suffix=.txt

echo "PHACTS formating done"

#Run the array# 

#After array:  Concatenate all the files.

# for file in *; do tail -n +2 "$file" | awk -v fname="$file" '{print fname, $0}' >> ../PHACTS_concatenated_output.txt; done

###Host prediction###

#Iphop
mkdir -p ${CONTIG_INFO_DIR}/iphop
mkdir -p ${CONTIG_INFO_DIR}/iphop/MAGs_split
mkdir -p ${CONTIG_INFO_DIR}/iphop/outputPerSplit

seqkit split2 --by-size 100 -O ${CONTIG_INFO_DIR}/iphop/MAGs_split ${PHAGES}

ls ${CONTIG_INFO_DIR}/iphop/MAGs_split > ${CONTIG_INFO_DIR}/iphop/lookup_magsSplit.txt

#Run the array# 

#After array: Find all CSV files and concatenate them

# find ${CONTIG_INFO_DIR}/iphop/outputPerSplit -name 'Host_prediction_to_genome_m90.csv' -exec cat {} + > ${CONTIG_INFO_DIR}/iphop/iphopPred_allCombined_toGenome.txt &
# find ${CONTIG_INFO_DIR}/iphop/outputPerSplit -name 'Host_prediction_to_genus_m90.csv' -exec cat {} + > ${CONTIG_INFO_DIR}/iphop/iphopPred_allCombined_toGenus.txt &

#Remove duplicated headers except the first one
# awk 'NR == 1 || !/Virus,Host genome,Host taxonomy,Main method,Confidence score,Additional methods/' ${CONTIG_INFO_DIR}/iphop/iphopPred_allCombined_toGenome.txt > ${CONTIG_INFO_DIR}/iphop/iphopPred_allCombined_toGenome_final.csv &
# awk 'NR == 1 || !/Virus,Host genome,Host taxonomy,Main method,Confidence score,Additional methods/' ${CONTIG_INFO_DIR}/iphop/iphopPred_allCombined_toGenus.txt > ${CONTIG_INFO_DIR}/iphop/iphopPred_allCombined_toGenus_final.csv &

