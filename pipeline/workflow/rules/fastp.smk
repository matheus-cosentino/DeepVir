###################################################################################
#                        workflow/rules/fastp.smk                                 #
#                         MSc. Matheus Cosentino                                  # 
###################################################################################
#                                                                                 #
#         ████  █████ █████ ████  █   █ ███ ████                                  #
#         █░░░█ █░░░░░█░░░░░█░░░█ █░  █░ █░░█░░░█                                 #
#         █░░░█░████░░████░░████░░█░░ █░░█░░████░░                                #
#         █░░ █░█░░░░ █░░░░ █░░░░ ░█░█ ░░█░░█░░█░ ░                               #
#         ████ ░█████░█████░█░░░░░  █ ░ ███░█░░░█░                                #
#          ░░░░ ░░░░░░ ░░░░░ ░░      ░ ░ ░░░ ░░  ░                                 #
#           ░░░░  ░░░░░ ░░░░░ ░       ░   ░░░ ░   ░                                #
#                                                                                 #
###################################################################################
#                         version: 1.0  |  09.2026                              #
###################################################################################

############################################
# --- 1. Filter fastq.gz Paired Files --- #
########################################### 

#rule to process paired data
rule fastp_paired:
    input:
        r1 = get_input_r1,
        r2 = get_input_r2
    output:
        # Define all possible output files explicitly and statically ({sample}_{rule} standard).
        r1 = os.path.join(OUT_DIR, "{sample}", "fastp", "{sample}_fastp_paired_1.fastq.gz"),
        r2 = os.path.join(OUT_DIR, "{sample}", "fastp", "{sample}_fastp_paired_2.fastq.gz"),
        orphans = os.path.join(OUT_DIR, "{sample}", "fastp", "{sample}_fastp_paired_orphans.fastq.gz"),
        html = os.path.join(OUT_DIR, "{sample}", "fastp", "{sample}_fastp_paired.html"),
        json = os.path.join(OUT_DIR, "{sample}", "fastp", "{sample}_fastp_paired.json")
    log:
        os.path.join(LOG_DIR, "{sample}", "{sample}_fastp_paired.log")
    benchmark:
        os.path.join(OUT_DIR, "benchmarks", "fastp_paired", "{sample}_fastp_paired.tsv")
    params:
        length_required=lambda wc: config['fastp']['length_required'][0] if isinstance(config['fastp']['length_required'], (list, tuple)) else config['fastp']['length_required'],
        quality=lambda wc: config['fastp']['qualified_quality_phred'][0] if isinstance(config['fastp']['qualified_quality_phred'], (list, tuple)) else config['fastp']['qualified_quality_phred']
    conda:
        FASTP
    shell:
        """
        fastp \
            --in1 {input.r1} --in2 {input.r2} --out1 {output.r1} --out2 {output.r2} --unpaired1 {output.orphans} --unpaired2 {output.orphans} \
            --html {output.html} \
            --json {output.json} \
            --thread {resources.threads} \
            --length_required {params.length_required} \
            --qualified_quality_phred {params.quality} --dedup \
            2> {log}
        """

#if fail to make the paired reads cleanup, follown immediatelly to the fastp_unpaired rule
#############################################
# --- 2. Filter fastq.gz unpaired Files --- #
#############################################

rule fastp_unpaired:
    input:
        reads = get_input_unp,
    output:
        r1 = os.path.join(OUT_DIR, "{sample}", "fastp", "{sample}_fastp_unpaired.fastq.gz"),
        html = os.path.join(OUT_DIR, "{sample}", "fastp", "{sample}_fastp_unpaired.html"),
        json = os.path.join(OUT_DIR, "{sample}", "fastp", "{sample}_fastp_unpaired.json")
    log:
        os.path.join(LOG_DIR, "{sample}", "{sample}_fastp_unpaired.log")
    benchmark:
        os.path.join(OUT_DIR, "benchmarks", "fastp_unpaired", "{sample}_fastp_unpaired.tsv")
    params:
        is_nanopore = lambda wc: (config.get("platform") == "nanopore") or (isinstance(config.get("nanopore"), dict) and config.get("nanopore", {}).get("qc", {}).get("tool", "") == "chopper"),
        length_required = lambda wc: (config.get("nanopore", {}).get("min_read_len", 500) if isinstance(config.get("nanopore"), dict) else 500) if config.get("platform") == "nanopore" else (config['fastp']['length_required'][0] if isinstance(config['fastp']['length_required'], (list, tuple)) else config['fastp']['length_required']),
        quality = lambda wc: (config.get("nanopore", {}).get("qc", {}).get("min_qual", 10) if isinstance(config.get("nanopore"), dict) else 10) if config.get("platform") == "nanopore" else (config['fastp']['qualified_quality_phred'][0] if isinstance(config['fastp']['qualified_quality_phred'], (list, tuple)) else config['fastp']['qualified_quality_phred'])
    conda:
        QC_UNPAIRED_ENV
    shell:
        """
        mkdir -p "$(dirname "{output.r1}")"
        mkdir -p "$(dirname "{log}")"
        exec > "{log}" 2>&1

        if [ "{params.is_nanopore}" = "True" ] || [ "{params.is_nanopore}" = "true" ]; then
            echo "[INFO] Nanopore platform detected: Running Chopper for Long-Read QC & Length/Quality Filtering..."
            echo "[INFO] Parameters: Min Length = {params.length_required} bp, Min Phred Quality = Q{params.quality}"
            
            if [[ "{input.reads}" == *.gz ]]; then
                gzip -cd "{input.reads}" | chopper -q {params.quality} -l {params.length_required} --threads {resources.threads} | gzip -c > "{output.r1}"
            else
                cat "{input.reads}" | chopper -q {params.quality} -l {params.length_required} --threads {resources.threads} | gzip -c > "{output.r1}"
            fi

            RAW_COUNT=$(gzip -cd "{input.reads}" 2>/dev/null | grep -c "^@" || wc -l < "{input.reads}" || echo "0")
            CLEAN_COUNT=$(gzip -cd "{output.r1}" 2>/dev/null | grep -c "^@" || echo "0")
            
            echo "{{\\"summary\\": {{\"tool\\": \\"chopper\\", \\"raw_reads\\": $RAW_COUNT, \\"clean_reads\\": $CLEAN_COUNT, \\"min_qual\\": {params.quality}, \\"min_len\\": {params.length_required}}}}}" > "{output.json}"
            echo "<html><body><h2>Chopper Nanopore QC Summary</h2><p>Sample: {wildcards.sample}</p><p>Tool: Chopper (Rust long-read filter)</p><p>Filtered reads saved to {output.r1}</p></body></html>" > "{output.html}"
            echo "[INFO] Chopper QC complete. Output: {output.r1}"
        else
            echo "[INFO] Running fastp for Illumina short-read quality trimming..."
            fastp \
                --in1 {input.reads} --out1 {output.r1} \
                --html {output.html} \
                --json {output.json} \
                --thread {resources.threads} \
                --length_required {params.length_required} \
                --qualified_quality_phred {params.quality} \
                --dedup --trim_poly_g --trim_poly_x
        fi
        """
