###################################################################################
#                       workflow/rules/get_data.smk                               #
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

###################
# --- HELPERS --- #
###################

PAIRED_SRA = [s for s, m in SAMPLE_META.items() if m['mode'] == 'SRA' and len(m['files']) == 2]
SINGLE_SRA = [s for s, m in SAMPLE_META.items() if m['mode'] == 'SRA' and len(m['files']) == 1]

##############################################
# --- 1. Download SRA Data (PAIRED) --- #
############################################# 

rule download_sra_data_paired:
    output:
        r1 = os.path.join(config["data_dir"], "{sample}_1.fastq.gz"),
        r2 = os.path.join(config["data_dir"], "{sample}_2.fastq.gz")
    log:
        os.path.join(LOG_DIR, "{sample}", "{sample}_download.log")
    benchmark:
        os.path.join(config["output_dir"], "benchmarks", "download_sra", "{sample}.tsv")
    conda:
        DOWNLOAD
    shadow: 
        "minimal" 
    params:
        # FIX 1: Corrected Syntax with f-string
        tmpdir = lambda wildcards: os.path.join(config["data_dir"], f"tmp_{wildcards.sample}")
    wildcard_constraints:
        sample = "|".join(PAIRED_SRA) if PAIRED_SRA else "NO_PAIRED_SAMPLES"
    shell:
        """
        exec &> {log}
        echo "Step 1: Prefetching SRA file for {wildcards.sample}..."
        rm -rf ~/ncbi/public/sra/{wildcards.sample}.sra*
        prefetch {wildcards.sample} --max-size 100G
        
        echo "Step 2: Validating SRA integrity..."
        vdb-validate {wildcards.sample}

        echo "Step 3: Running fasterq-dump from local SRA file..."
        rm -rf {params.tmpdir}
        mkdir -p {params.tmpdir}
        fasterq-dump {wildcards.sample} \
            --split-files \
            --threads {threads} \
            --temp {params.tmpdir} \
            --outdir {params.tmpdir} \
            --skip-technical

        echo "Step 4: Compressing with pigz (faster) or gzip..."
        if command -v pigz >/dev/null 2>&1; then
            pigz -f -p {threads} "{params.tmpdir}/{wildcards.sample}_1.fastq"
            pigz -f -p {threads} "{params.tmpdir}/{wildcards.sample}_2.fastq"
        else
            gzip -f "{params.tmpdir}/{wildcards.sample}_1.fastq"
            gzip -f "{params.tmpdir}/{wildcards.sample}_2.fastq"
        fi
        
        echo "Step 5: Moving files and cleaning up..."
        mv "{params.tmpdir}/{wildcards.sample}_1.fastq.gz" {output.r1}
        mv "{params.tmpdir}/{wildcards.sample}_2.fastq.gz" {output.r2}
        
        # Remove the downloaded .sra file to save disk space in NCBI directory
        rm -rf ~/ncbi/public/sra/{wildcards.sample}.sra
        rm -rf {params.tmpdir}
        """

##############################################
# --- 2. Download SRA Data (SINGLE) --- #
############################################# 

rule download_sra_single:
    output:
        r1 = os.path.join(config["data_dir"], "{sample}_1.fastq.gz")
    log:
        os.path.join(LOG_DIR, "{sample}", "{sample}_download_single.log")
    benchmark:
        os.path.join(config["output_dir"], "benchmarks", "download_sra_single", "{sample}.tsv")
    conda:
        DOWNLOAD
    wildcard_constraints:
        sample = "|".join(SINGLE_SRA) if SINGLE_SRA else "NO_SINGLE_SAMPLES"
    shadow: 
        "minimal"
    # FIX 4: Added missing params section
    params:
        tmpdir = lambda wildcards: os.path.join(config["data_dir"], f"tmp_{wildcards.sample}")
    shell:
        """
        echo "Creating temp dir for download {wildcards.sample}..." > {log}
        mkdir -p {params.tmpdir}

        echo "Starting SingleEnd download for {wildcards.sample}..." >> {log}
        fasterq-dump --split-files --threads {resources.threads} -O {params.tmpdir} {wildcards.sample} >> {log} 2>&1

        echo "Compressing and renaming..." >> {log}
        
        # Logic to handle naming variations (sample.fastq vs sample_1.fastq)
        # Note: We must compress BEFORE moving to .gz output
        
        if command -v pigz >/dev/null 2>&1; then
            COMPRESS="pigz -f -p {resources.threads}"
        else
            COMPRESS="gzip -f"
        fi

        if [ -f "{params.tmpdir}/{wildcards.sample}_1.fastq" ]; then
            $COMPRESS "{params.tmpdir}/{wildcards.sample}_1.fastq"
            mv "{params.tmpdir}/{wildcards.sample}_1.fastq.gz" {output.r1}
        elif [ -f "{params.tmpdir}/{wildcards.sample}.fastq" ]; then
            $COMPRESS "{params.tmpdir}/{wildcards.sample}.fastq"
            mv "{params.tmpdir}/{wildcards.sample}.fastq.gz" {output.r1}
        fi

        echo "Cleaning up temp dir..." >> {log}
        rm -rf {params.tmpdir}
        """