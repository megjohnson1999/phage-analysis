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
        "bioconda::seqkit"
    shell:
        """
        # Create output directory
        mkdir -p {output.split_dir}
        
        # Split FASTA file into individual files
        seqkit split {input.phage_seqs} -O {output.split_dir} -s 1 > {log} 2>&1
        
        # Create list of split files
        find {output.split_dir} -name "*.fasta" > {output.split_list}
        """

# 2. Run iPhop for host prediction on each split file
rule iphop_host_prediction:
    input:
        split_list = f"{config['output_dir']}/03_split_seqs/split_file_list.txt"
    output:
        results_dir = directory(f"{config['output_dir']}/03_iphop_results"),
        predictions = f"{config['output_dir']}/03_iphop_results/iphop_predictions_compiled.tsv"
    log:
        f"{config['output_dir']}/logs/iphop_host_prediction.log"
    conda:
        config["conda_envs"]["iphop"]
    resources:
        mem_mb = config["resources"]["iphop"]["mem_mb"],
        threads = config["resources"]["iphop"]["threads"],
        time = config["resources"]["iphop"]["time"]
    shell:
        """
        # Create output directory
        mkdir -p {output.results_dir}/tmp
        
        # Process each split file
        while read phage_file; do
            filename=$(basename "$phage_file")
            outdir="{output.results_dir}/tmp/$(basename "$filename" .fasta)"
            
            echo "Processing $filename" >> {log}
            
            # Run iPhop
            iphop predict --fa_file "$phage_file" \
                --out_dir "$outdir" \
                --num_threads {resources.threads} >> {log} 2>&1
        done < {input.split_list}
        
        # Compile results
        echo "Compiling results" >> {log}
        head -n 1 $(find {output.results_dir}/tmp -name "host_prediction_to_genus.csv" | head -n 1) > {output.predictions}.tmp
        find {output.results_dir}/tmp -name "host_prediction_to_genus.csv" | xargs cat | grep -v "query" >> {output.predictions}.tmp
        
        # Convert to TSV
        tr ',' '\\t' < {output.predictions}.tmp > {output.predictions}
        rm {output.predictions}.tmp
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
        "bioconda::prodigal"
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
        "bioconda::seqkit"
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

# 5. Run PHACTS for lifestyle prediction
rule phacts_lifestyle_prediction:
    input:
        split_list = f"{config['output_dir']}/03_split_proteins/split_protein_list.txt"
    output:
        results_dir = directory(f"{config['output_dir']}/03_phacts_results"),
        predictions = f"{config['output_dir']}/03_phacts_results/phacts_predictions_compiled.tsv"
    log:
        f"{config['output_dir']}/logs/phacts_lifestyle_prediction.log"
    conda:
        config["conda_envs"]["phacts"]
    resources:
        mem_mb = config["resources"]["phacts"]["mem_mb"],
        threads = config["resources"]["phacts"]["threads"],
        time = config["resources"]["phacts"]["time"]
    shell:
        """
        # Create output directory
        mkdir -p {output.results_dir}/tmp
        
        # Process each split file
        while read protein_file; do
            filename=$(basename "$protein_file")
            phage_id=$(basename "$filename" .faa)
            outfile="{output.results_dir}/tmp/$phage_id.phacts.out"
            
            echo "Processing $phage_id" >> {log}
            
            # Run PHACTS
            phacts.py "$protein_file" "$outfile" >> {log} 2>&1
        done < {input.split_list}
        
        # Compile results
        echo -e "phage_id\\tlifestyle\\tprobability" > {output.predictions}
        for file in {output.results_dir}/tmp/*.phacts.out; do
            phage_id=$(basename "$file" .phacts.out)
            lifestyle=$(grep "Lifestyle:" "$file" | awk '{{print $2}}')
            probability=$(grep "Probability:" "$file" | awk '{{print $2}}')
            echo -e "$phage_id\\t$lifestyle\\t$probability" >> {output.predictions}
        done
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
    resources:
        mem_mb = config["resources"]["genomic_info"]["mem_mb"],
        threads = config["resources"]["genomic_info"]["threads"],
        time = config["resources"]["genomic_info"]["time"]
    shell:
        """
        # Create temporary directory
        TMP_DIR=$(mktemp -d)
        
        # Run mmseqs2 for taxonomy assignment
        mmseqs easy-taxonomy {input.phage_seqs} {config[resources][mmseqs2][db]} \
            {output.results_dir} $TMP_DIR \
            --threads {resources.threads} \
            --lca-ranks species,genus,family,order,class,phylum,superkingdom \
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
        results_dir = directory(f"{config['output_dir']}/03_genomic_info/phabox_output")
    log:
        f"{config['output_dir']}/logs/phabox_prediction.log"
    conda:
        config["conda_envs"]["phabox2"]
    resources:
        mem_mb = config["resources"]["genomic_info"]["mem_mb"],
        threads = config["resources"]["genomic_info"]["threads"],
        time = config["resources"]["genomic_info"]["time"]
    shell:
        """
        # Run Phabox2
        phabox.py -i {input.phage_seqs} -o {output.results_dir} \
            -t {resources.threads} > {log} 2>&1
        """

# 8. Run vContact3 for phage taxonomy based on gene content
rule vcontact3_taxonomy:
    input:
        proteins = f"{config['output_dir']}/03_orf_predictions/proteins.faa"
    output:
        results_dir = directory(f"{config['output_dir']}/03_genomic_info/vc3_output")
    log:
        f"{config['output_dir']}/logs/vcontact3_taxonomy.log"
    conda:
        config["conda_envs"]["vcontact3"]
    resources:
        mem_mb = config["resources"]["genomic_info"]["mem_mb"],
        threads = config["resources"]["genomic_info"]["threads"],
        time = config["resources"]["genomic_info"]["time"]
    shell:
        """
        # Create gene2genome file
        echo -e "protein_id\\tcontig_id" > {output.results_dir}/gene2genome.csv
        grep ">" {input.proteins} | sed 's/>//g' | 
        awk -F " # " '{{split($1,a,"_"); print $1"\\t"a[1]}}' >> {output.results_dir}/gene2genome.csv
        
        # Run vContact3
        vContact3 --raw-proteins {input.proteins} \
            --rel-mode 'Diamond' \
            --proteins-fp {output.results_dir}/gene2genome.csv \
            --output-dir {output.results_dir} \
            --threads {resources.threads} > {log} 2>&1
        """