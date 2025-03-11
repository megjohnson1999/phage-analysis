import pandas as pd
import argparse
from Bio import SeqIO

def filter_mmseqs_output(input_file, output_file, ids_output_file, fasta_file, missing_contigs_output_file):
    #Load the table
    df = pd.read_csv(input_file, sep="\t", header=None, names=[
        "contig_id", "tax_id", "rank", "label", "fragments", "assigned_fragments",
        "label_fragments", "support", "lineage"
    ])

    #Condition 1: First taxID in lineage is 10239 (virus)
    df["first_taxid_in_lineage"] = df["lineage"].apply(lambda x: int(x.split(";")[0]) if pd.notnull(x) else None)
    condition_virus = df["first_taxid_in_lineage"] == 10239

    #Condition 2: taxID is 1
    condition_taxid_one = df["tax_id"] == 1

    #Condition 3: taxID is 0
    condition_taxid_zero = df["tax_id"] == 0

    #Condition 4: Entries not in the table (from FASTA file)
    contig_ids_in_table = set(df["contig_id"])
    contig_ids_in_fasta = {record.id for record in SeqIO.parse(fasta_file, "fasta")}
    missing_contigs = contig_ids_in_fasta - contig_ids_in_table

    #Save missing contigs to a file
    with open(missing_contigs_output_file, "w") as mc_file:
        for contig_id in sorted(missing_contigs):
            mc_file.write(f"{contig_id}\n")

    #Apply conditions and filter the dataframe
    filtered_df = df[condition_virus | condition_taxid_one | condition_taxid_zero]
    passing_contig_ids = set(filtered_df["contig_id"]).union(missing_contigs)

    #Save filtered table
    filtered_df.to_csv(output_file, sep="\t", index=False, header=False)
    
    #Save passing contig IDs
    with open(ids_output_file, "w") as ids_file:
        for contig_id in sorted(passing_contig_ids):
            ids_file.write(f"{contig_id}\n")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Filter mmseqs2 easy-taxonomy output based on specific criteria.")
    parser.add_argument("--mmseqs_LCA_table", required=True, help="Path to the input mmseqs2 output table.")
    parser.add_argument("--o_filtered_LCA_table", required=True, help="Path to the output filtered table.")
    parser.add_argument("--o_passing_contig_ids", required=True, help="Path to the output file for filtered contig IDs.")
    parser.add_argument("--contigs", required=True, help="Path to the FASTA file of contigs.")
    parser.add_argument("--o_missing_contig_ids", required=True, help="Path to the output file for missing contig IDs.")

    args = parser.parse_args()
    filter_mmseqs_output(args.mmseqs_LCA_table, args.o_filtered_LCA_table, args.o_passing_contig_ids, args.contigs, args.o_missing_contig_ids)