#!/bin/bash

#SBATCH --mem=600GB
#SBATCH --cpus-per-task=24
#SBATCH --time=300:00:00
#SBATCH --output=MTBharris_phageNewWF.o%j 

#Exit immediately if a command exits with a non-zero status
# set -e

#Activate environment
source /ref/sahlab/software/miniforge3/bin/activate

#Define main paths for easier maintenance
SCRIPT_DIR="/home/luisalberto/Scripts/phageAnalysisPipeline"
REF_DIR="/ref/sahlab/data"
DB="${REF_DIR}/nr/nr_mmseqs15_DB"
GENOMAD_DB="${REF_DIR}/genomad_db"
PHOLD_DB="${REF_DIR}/protein_Structure_DBs/phold"

###RC2###
BASE_DIR="/scratch/sahlab/Luis/RC2_hecatomb/2023_04_07_RC2_Freeze4_HecatombOut_WholeRun_MMseqsFast/test_newPhageWorkflowCeliac"

#Input files and directories
HECATOMBASSEMBLY="/lts/sahlab/data3/2023_04_07_RC2_Freeze4_HecatombOut_WholeRun_MMseqsFast/processing/assembly/CONTIG_DICTIONARY/FLYE/assembly_graph.gfa"
READS="/lts/sahlab/data4/luis/RC2_analysisProgress/01_cleanAllFastqForMappingCov"

#Output directories
RENEO_OUT="${BASE_DIR}/01_reneo"
MMSEQS_OUT="${BASE_DIR}/02_mmseqs"
TMP_DIR="/scratch/sahlab/Luis/temporal/"
FILTERING_WD="${BASE_DIR}/03_selectingViruses"
PHAGEPRED_WD="${BASE_DIR}/04_phagePred"
QUALPRED_WD="${BASE_DIR}/05_qualPred"
PHAGE_INFO_WD="${BASE_DIR}/06_phages_cluster"

#Step 1: Binning using reneo. 
#IMPORTANT: For large datasets, make sure to have the reneo and Koverage config files in the directory with proper resources. 
conda activate reneo
reneo run --input "${HECATOMBASSEMBLY}" \
          --reads "${READS}" \
          --minlength 1000 \
          --output "${RENEO_OUT}" \
          --threads ${SLURM_CPUS_PER_TASK}

seqkit seq --min-len 1000 -g \
"${RENEO_OUT}/genomes_and_unresolved_edges.fasta" > "${RENEO_OUT}/genomes_and_unresolved_edges_1KB.fasta"

echo "reneo done"

#Step 2: mmseqs LCA for selecting root, NAs, and viral contigs
mkdir -p "${MMSEQS_OUT}"
mkdir -p "${TMP_DIR}"

conda activate mmseqs2_v15-6f452
mmseqs easy-taxonomy "${RENEO_OUT}/genomes_and_unresolved_edges_1KB.fasta" \
       "${DB}" \
       "${MMSEQS_OUT}/genomes_and_unresolved_edges_mmseqs" \
       "${TMP_DIR}" \
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

echo "mmseqs2 done"

#Step 3: Selecting contigs and generating viral fasta file
mkdir -p "${FILTERING_WD}"

python "${SCRIPT_DIR}/02_filter_mmseqsLCA.py" --mmseqs_LCA_table "${MMSEQS_OUT}/genomes_and_unresolved_edges_mmseqs_lca.tsv" \
       --contigs "${RENEO_OUT}/genomes_and_unresolved_edges_1KB.fasta" \
       --o_filtered_LCA_table "${FILTERING_WD}/filtered_output.txt" \
       --o_passing_contig_ids "${FILTERING_WD}/passing_contig_ids.txt" \
       --o_missing_contig_ids "${FILTERING_WD}/contigsInFastaNotInLCA.txt"

Step 4: Get fasta for the contigs
seqkit grep -f "${FILTERING_WD}/passing_contig_ids.txt" "${RENEO_OUT}/genomes_and_unresolved_edges_1KB.fasta" > "${FILTERING_WD}/passing_Viralcontigs.fasta"

echo "viral fasta done" 

#Step 5: Phage prediction
mkdir -p "${PHAGEPRED_WD}"

conda activate jaeger
Jaeger -i "${FILTERING_WD}/passing_Viralcontigs.fasta"  \
       -o "${PHAGEPRED_WD}/jaeger" \
       -s 2.5 \
       --fsize 1000 \
       --stride 1000

echo "jaeger done" 

conda activate /ref/sahlab/software/anaconda3/envs/genomad
genomad end-to-end --min-score 0.6 \
       --cleanup \
       --threads ${SLURM_CPUS_PER_TASK} \
       "${FILTERING_WD}/passing_Viralcontigs.fasta" \
       "${PHAGEPRED_WD}/geNomad" \
       "${GENOMAD_DB}"

echo "genomad done" 

#Step 6: Phold from phage prediction
conda activate pholdENV
phold run -i "${FILTERING_WD}/passing_Viralcontigs.fasta" \
       -o "${PHAGEPRED_WD}/phold" \
       -d "${PHOLD_DB}" \
       -t ${SLURM_CPUS_PER_TASK} --cpu --force

echo "phold done" 

#If the last step fail beause of memory. Running juts phold compare 
# phold compare -i ${FILTERING_WD}/passing_Viralcontigs.fasta \
#      --predictions_dir ${PHAGEPRED_WD}/phold \
#      -o ${PHAGEPRED_WD}/phold/pholdCompare -t 24 \
#      -d "${PHOLD_DB}" --force

#Step 7: Run checkV on the passing_ViralContigs file
mkdir -p "${QUALPRED_WD}"

conda activate /ref/sahlab/software/anaconda3/envs/checkV 
checkv end_to_end "${FILTERING_WD}/passing_Viralcontigs.fasta" \
       "${QUALPRED_WD}/passing_Viralcontigs_checkV" \
       -d /ref/sahlab/data/viral_analysis_DBs/checkV_DB/checkv-db-v1.5 \
       -t ${SLURM_CPUS_PER_TASK}

echo "checkV passing_ViralContigs done" 

#step 8: do the parsing using all the prediction informations and functional annotation and create the final phage set
mkdir -p "${PHAGE_INFO_WD}"
conda activate R

Rscript "${SCRIPT_DIR}/02_phagePrediction.R" \
"${PHAGEPRED_WD}/phold/pholdCompare/phold_per_cds_predictions.tsv" \
"${PHAGEPRED_WD}/jaeger/passing_Viralcontigs_default_jaeger.tsv" \
"${PHAGEPRED_WD}/geNomad/passing_Viralcontigs_summary/passing_Viralcontigs_virus_summary.tsv" \
"${QUALPRED_WD}/passing_Viralcontigs_checkV/quality_summary.tsv" \
"${PHAGE_INFO_WD}"

seqkit grep -f "${PHAGE_INFO_WD}/contig_ids.txt" "${FILTERING_WD}/passing_Viralcontigs.fasta" > "${PHAGE_INFO_WD}/phageContigs.fasta"

echo "phageContigs identification done" 

#final checkV checking 
conda activate /ref/sahlab/software/anaconda3/envs/checkV 
checkv end_to_end "${PHAGE_INFO_WD}/phageContigs.fasta" \
       "${QUALPRED_WD}/phageContigs_checkV" \
       -d /ref/sahlab/data/viral_analysis_DBs/checkV_DB/checkv-db-v1.5 \
       -t ${SLURM_CPUS_PER_TASK}

echo "checkV phageContigs done" 