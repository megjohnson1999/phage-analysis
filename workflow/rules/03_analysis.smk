"""
Rules for analyzing phage genomes, including host prediction and genomic characterization.
These rules can take input from either clustering or prediction steps.
"""

# 1. Split phage sequences for array-based processing
rule split_phage_sequences:
    input:
        phage_seqs = get_phage_input  # Function defined in 02_clustering.smk
    output:
        split_dir = directory(f"{config['output_dir']}/03_split_seqs"),
        split_list = f"{config['output_dir']}/03_split_seqs/split_file_list.txt"
    log:
        f"{config['output_dir']}/logs/split_phage_sequences.log"
    conda:
        config["conda_envs"]["seqkit"]
    shell:
        """
        # Create output directory
        mkdir -p {output.split_dir}
        
        # Split FASTA file into individual files
        seqkit split {input.phage_seqs} -O {output.split_dir} -s 1 > {log} 2>&1
        
        # Create list of split files
        find {output.split_dir} -name "*.fasta" > {output.split_list}
        """

# List all samples for iphop from split files
def get_iphop_samples():
    # After split_phage_sequences is run, this reads the split file list
    split_list = f"{config['output_dir']}/03_split_seqs/split_file_list.txt"
    if os.path.exists(split_list):
        with open(split_list, "r") as f:
            files = [line.strip() for line in f]
        return [os.path.splitext(os.path.basename(file))[0] for file in files]
    # Fallback to glob when file doesn't exist yet (for dry runs)
    elif os.path.exists(f"{config['output_dir']}/03_split_seqs"):
        return [os.path.splitext(os.path.basename(f))[0] 
                for f in glob.glob(f"{config['output_dir']}/03_split_seqs/*.fasta")]
    else:
        return []  # Return empty list if no directory exists yet

# 2. iphop workflow - Use checkpoint to wait for split files to be created
checkpoint wait_for_iphop_splits:
    input:
        split_list = f"{config['output_dir']}/03_split_seqs/split_file_list.txt"
    output:
        touch(f"{config['output_dir']}/03_iphop_results/.splits_ready")
    shell:
        "mkdir -p $(dirname {output})"

# Checkpoint to handle dynamic input files in a way that works with dry-run
checkpoint get_iphop_input_files:
    input:
        split_list = f"{config['output_dir']}/03_split_seqs/split_file_list.txt"
    output:
        flag = f"{config['output_dir']}/03_iphop_results/.input_files_found"
    shell:
        """
        mkdir -p $(dirname {output.flag})
        touch {output.flag}
        """

# 2a. Run iPhop for host prediction on a single split file
rule iphop_single_prediction:
    input:
        checkpoint = f"{config['output_dir']}/03_iphop_results/.splits_ready",
        flag = f"{config['output_dir']}/03_iphop_results/.input_files_found",
        phage_file = f"{config['output_dir']}/03_split_seqs/{{sample}}.fasta"
    output:
        prediction = f"{config['output_dir']}/03_iphop_results/tmp/{{sample}}/host_prediction_to_genus.csv"
    log:
        f"{config['output_dir']}/logs/iphop_prediction/{{sample}}.log"
    conda:
        config["conda_envs"]["iphop"]
    threads: 12
    shell:
        """
        # Create output directory
        mkdir -p $(dirname {output.prediction})
        
        # Run iPhop
        iphop predict --fa_file {input.phage_file} \
            --db_dir {config[databases][iphop][db]} \
            --out_dir $(dirname {output.prediction}) \
            --num_threads {threads} > {log} 2>&1
        """

# Helper rule to force running single predictions during dry run
rule run_all_iphop_predictions:
    input:
        checkpoint = f"{config['output_dir']}/03_iphop_results/.splits_ready",
        # For actual runs, get samples from the split files
        # For dry runs, this will be an empty list, which is fine
        samples = lambda wildcards: expand(
            f"{config['output_dir']}/03_iphop_results/tmp/{{sample}}/host_prediction_to_genus.csv",
            sample=get_iphop_samples()
        )
    output:
        touch(f"{config['output_dir']}/03_iphop_results/.all_predictions_done")

# 2b. Aggregate iPhop results
rule iphop_aggregate_results:
    input:
        # This is the key part that makes the parallelization work
        # Aggregation only happens after all individual predictions are done
        all_done = f"{config['output_dir']}/03_iphop_results/.all_predictions_done",
        predictions = lambda wildcards: expand(
            f"{config['output_dir']}/03_iphop_results/tmp/{{sample}}/host_prediction_to_genus.csv",
            sample=get_iphop_samples()
        )
    output:
        predictions = f"{config['output_dir']}/03_iphop_results/iphop_predictions_compiled.tsv"
    log:
        f"{config['output_dir']}/logs/iphop_aggregate_results.log"
    shell:
        """
        # Ensure results directory exists
        RESULTS_DIR=$(dirname {output.predictions})
        mkdir -p $RESULTS_DIR
        
        # Compile results
        echo "Compiling results" > {log}
        
        if [ -n "$(ls -A $RESULTS_DIR/tmp 2>/dev/null)" ]; then
            # If there are prediction files
            FIRST_FILE=$(find $RESULTS_DIR/tmp -name "host_prediction_to_genus.csv" | head -n 1)
            
            if [ -n "$FIRST_FILE" ]; then
                head -n 1 "$FIRST_FILE" > {output.predictions}.tmp
                find $RESULTS_DIR/tmp -name "host_prediction_to_genus.csv" | xargs cat | grep -v "query" >> {output.predictions}.tmp
            else
                # Create empty file with header
                echo "query,host,score,identity,coverage,kingdom,phylum,class,order,family,genus" > {output.predictions}.tmp
            fi
            
            # Convert to TSV
            tr ',' '\t' < {output.predictions}.tmp > {output.predictions}
            rm {output.predictions}.tmp
        else
            # Create empty file with header
            echo -e "query\thost\tscore\tidentity\tcoverage\tkingdom\tphylum\tclass\torder\tfamily\tgenus" > {output.predictions}
            echo "No prediction files found, created empty result file" >> {log}
        fi
        """

# 3. Run Prodigal for ORF prediction on phage sequences
rule prodigal_orf_prediction:
    input:
        phage_seqs = get_phage_input
    output:
        proteins = f"{config['output_dir']}/03_orf_predictions/proteins.faa",
        genes = f"{config['output_dir']}/03_orf_predictions/genes.fna"
    log:
        f"{config['output_dir']}/logs/prodigal_orf_prediction.log"
    conda:
        config["conda_envs"]["phacts"]
    shell:
        """
        # Run Prodigal for ORF prediction
        prodigal -i {input.phage_seqs} \
            -a {output.proteins} \
            -d {output.genes} \
            -p meta > {log} 2>&1
        """

# 4. Split protein files for PHACTS array processing
rule split_protein_files:
    input:
        proteins = f"{config['output_dir']}/03_orf_predictions/proteins.faa"
    output:
        split_dir = directory(f"{config['output_dir']}/03_split_proteins"),
        split_list = f"{config['output_dir']}/03_split_proteins/split_protein_list.txt"
    log:
        f"{config['output_dir']}/logs/split_protein_files.log"
    conda:
        config["conda_envs"]["seqkit"]
    shell:
        """
        # Create output directory
        mkdir -p {output.split_dir}
        
        # Extract protein files by contig/phage
        awk '/^>/ {{if (seqlen) print seqlen; print; seqlen=0; next}} {{seqlen+=length}} END {{print seqlen}}' {input.proteins} | 
        awk 'BEGIN {{count=0}} 
             /^>/ {{
                if (count>0) close(file); 
                count++; 
                match($0, />(\\S+)/, arr); 
                file="{output.split_dir}/"arr[1]".faa"; 
                print $0 > file; 
                next
             }} 
             {{print > file}}' >> {log} 2>&1
        
        # Create list of split files
        find {output.split_dir} -name "*.faa" > {output.split_list}
        """

# List all samples for phacts from split protein files
def get_phacts_samples():
    # After split_protein_files is run, this reads the split file list
    split_list = f"{config['output_dir']}/03_split_proteins/split_protein_list.txt"
    if os.path.exists(split_list):
        with open(split_list, "r") as f:
            files = [line.strip() for line in f]
        return [os.path.splitext(os.path.basename(file))[0] for file in files]
    # Fallback to glob when file doesn't exist yet (for dry runs)
    elif os.path.exists(f"{config['output_dir']}/03_split_proteins"):
        return [os.path.splitext(os.path.basename(f))[0] 
                for f in glob.glob(f"{config['output_dir']}/03_split_proteins/*.faa")]
    else:
        return []  # Return empty list if no directory exists yet

# 5. PHACTS workflow - Use checkpoint to wait for split files to be created
checkpoint wait_for_phacts_splits:
    input:
        split_list = f"{config['output_dir']}/03_split_proteins/split_protein_list.txt"
    output:
        touch(f"{config['output_dir']}/03_phacts_results/.splits_ready")
    shell:
        "mkdir -p $(dirname {output})"

# Checkpoint to handle dynamic input files in a way that works with dry-run
checkpoint get_phacts_input_files:
    input:
        split_list = f"{config['output_dir']}/03_split_proteins/split_protein_list.txt"
    output:
        flag = f"{config['output_dir']}/03_phacts_results/.input_files_found"
    shell:
        """
        mkdir -p $(dirname {output.flag})
        touch {output.flag}
        """

# 5a. Run PHACTS for lifestyle prediction on a single protein file
rule phacts_single_prediction:
    input:
        checkpoint = f"{config['output_dir']}/03_phacts_results/.splits_ready",
        flag = f"{config['output_dir']}/03_phacts_results/.input_files_found",
        protein_file = f"{config['output_dir']}/03_split_proteins/{{sample}}.faa"
    output:
        result = f"{config['output_dir']}/03_phacts_results/tmp/{{sample}}.phacts.out"
    threads: 4
    log:
        f"{config['output_dir']}/logs/phacts_prediction/{{sample}}.log"
    conda:
        config["conda_envs"]["phacts"]
    shell:
        """
        # Create output directory
        mkdir -p $(dirname {output.result})
        
        # Run PHACTS
        phacts.py {input.protein_file} {output.result} > {log} 2>&1
        """

# Helper rule to force running single phacts predictions during dry run
rule run_all_phacts_predictions:
    input:
        checkpoint = f"{config['output_dir']}/03_phacts_results/.splits_ready",
        # For actual runs, get samples from the split files
        # For dry runs, this will be an empty list, which is fine
        samples = lambda wildcards: expand(
            f"{config['output_dir']}/03_phacts_results/tmp/{{sample}}.phacts.out",
            sample=get_phacts_samples()
        )
    output:
        touch(f"{config['output_dir']}/03_phacts_results/.all_predictions_done")

# 5b. Aggregate PHACTS results
rule phacts_aggregate_results:
    input:
        # This is the key part that makes the parallelization work
        # Aggregation only happens after all individual predictions are done
        all_done = f"{config['output_dir']}/03_phacts_results/.all_predictions_done",
        predictions = lambda wildcards: expand(
            f"{config['output_dir']}/03_phacts_results/tmp/{{sample}}.phacts.out",
            sample=get_phacts_samples()
        )
    output:
        predictions = f"{config['output_dir']}/03_phacts_results/phacts_predictions_compiled.tsv"
    log:
        f"{config['output_dir']}/logs/phacts_aggregate_results.log"
    shell:
        """
        # Ensure results directory exists
        RESULTS_DIR=$(dirname {output.predictions})
        mkdir -p $RESULTS_DIR
        
        # Compile results
        echo "Compiling PHACTS results" > {log}
        echo -e "phage_id\\tlifestyle\\tprobability" > {output.predictions}
        
        if [ -n "$(ls -A $RESULTS_DIR/tmp 2>/dev/null)" ]; then
            # Process each phacts output file if it exists
            for file in $RESULTS_DIR/tmp/*.phacts.out 2>/dev/null; do
                if [ -f "$file" ]; then
                    phage_id=$(basename "$file" .phacts.out)
                    lifestyle=$(grep "Lifestyle:" "$file" | awk '{{print $2}}')
                    probability=$(grep "Probability:" "$file" | awk '{{print $2}}')
                    echo -e "$phage_id\\t$lifestyle\\t$probability" >> {output.predictions}
                fi
            done
            
            if [ "$(wc -l < {output.predictions})" -eq 1 ]; then
                echo "No valid PHACTS results were found, only header in output file" >> {log}
            else
                echo "Successfully compiled $(( $(wc -l < {output.predictions}) - 1 )) PHACTS results" >> {log}
            fi
        else
            echo "No PHACTS result files found, created empty result file" >> {log}
        fi
        """

# 6. Run mmseqs2 taxonomy assignment on phage genomes
rule mmseqs_phage_taxonomy:
    input:
        phage_seqs = get_phage_input
    output:
        results_dir = directory(f"{config['output_dir']}/03_genomic_info/mmseqs_output"),
        taxonomy = f"{config['output_dir']}/03_genomic_info/mmseqs_taxonomy.tsv"
    log:
        f"{config['output_dir']}/logs/mmseqs_phage_taxonomy.log"
    conda:
        config["conda_envs"]["mmseqs2"]
    threads: 24
    shell:
        """
        # Create temporary directory
        TMP_DIR=$(mktemp -d)
        
        # Run mmseqs2 for taxonomy assignment
        mmseqs easy-taxonomy {input.phage_seqs} {config[databases][mmseqs2][db]} \
            {output.results_dir} $TMP_DIR \
            --min-length 30 \
            -e 1e-15 \
            --search-type 2 \
            -s 4.0 \
            --shuffle 0 \
            --lca-mode 2 \
            -a \
            --tax-lineage 2 \
            --format-output "query,target,evalue,pident,fident,nident,mismatch,qcov,tcov,qstart,qend,qlen,tstart,tend,tlen,alnlen,bits,qheader,theader,taxid,taxname,taxlineage" \
            --threads {threads} \
            --split-mode 0 \
            --orf-filter 1 \
            > {log} 2>&1
            
        # Copy and format the taxonomy results
        cp {output.results_dir}/lca.tsv {output.taxonomy}
            
        # Clean up
        rm -rf $TMP_DIR
        """

# 7. Run Phabox2 for phage taxonomy and lifestyle prediction
rule phabox_prediction:
    input:
        phage_seqs = get_phage_input
    output:
        results_dir = directory(f"{config['output_dir']}/03_genomic_info/phabox_output"),
        taxonomy = f"{config['output_dir']}/03_genomic_info/phabox_output/taxonomy.tsv",
        lifestyle = f"{config['output_dir']}/03_genomic_info/phabox_output/lifestyle.tsv"
    log:
        f"{config['output_dir']}/logs/phabox_prediction.log"
    conda:
        config["conda_envs"]["phabox2"]
    threads: 24
    shell:
        """
        # Run Phabox2
        phabox2 --task end_to_end --dbdir {config[databases][phabox][db]} \
            --outpth {output.results_dir} \
            --contigs {input.phage_seqs} \
            --len 1000 \
            --threads {threads} > {log} 2>&1
        """

# 8. Run vContact3 for phage taxonomy based on gene content
rule vcontact3_taxonomy:
    input:
        proteins = f"{config['output_dir']}/03_orf_predictions/proteins.faa"
    output:
        results_dir = directory(f"{config['output_dir']}/03_genomic_info/vc3_output"),
        gene2genome = f"{config['output_dir']}/03_genomic_info/vc3_output/gene2genome.csv",
        clusters = f"{config['output_dir']}/03_genomic_info/vc3_output/genome_by_genome_overview.csv"
    log:
        f"{config['output_dir']}/logs/vcontact3_taxonomy.log"
    conda:
        config["conda_envs"]["vcontact3"]
    threads: 24
    shell:
        """
        # Create gene2genome file
        echo -e "protein_id\\tcontig_id" > {output.gene2genome}
        grep ">" {input.proteins} | sed 's/>//g' | 
        awk -F " # " '{{split($1,a,"_"); print $1"\\t"a[1]}}' >> {output.gene2genome}
        
        # Run vContact3
        vcontact3 run --nucleotide $(dirname {input.proteins})/../{get_phage_input(None).split('/')[-1]} \
            --output {output.results_dir} \
            --db-domain "prokaryotes" \
            --db-version 223 \
            --db-path {config[databases][vcontact3][db]} \
            -t {threads} > {log} 2>&1
        """
