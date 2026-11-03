###################################################################################
#                       workflow/rules/kraken2.smk                                #
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

_KRAKEN2_DB_DIR = _os.path.join(_PROJECT_ROOT, config["resources"]["kraken2"][0])
_KRAKEN2_CORE_FILES = [
    os.path.join(_KRAKEN2_DB_DIR, "hash.k2d"),
    os.path.join(_KRAKEN2_DB_DIR, "opts.k2d"),
    os.path.join(_KRAKEN2_DB_DIR, "taxo.k2d"),
]


rule download_default_kraken2_db:
    """
    Automatically downloads and extracts the base Kraken2 RefSeq Viral database
    (k2_viral_20240112.tar.gz from the First Time Local Usage tutorial) on first
    run when no external or pre-existing Kraken2 database is present in resources/kraken2/.
    """
    output:
        hash_k2d = os.path.join(_KRAKEN2_DB_DIR, "hash.k2d"),
        opts_k2d = os.path.join(_KRAKEN2_DB_DIR, "opts.k2d"),
        taxo_k2d = os.path.join(_KRAKEN2_DB_DIR, "taxo.k2d")
    params:
        db_dir = _KRAKEN2_DB_DIR,
        tar_gz = os.path.join(_KRAKEN2_DB_DIR, "k2_viral_20240112.tar.gz"),
        url = "https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20240112.tar.gz"
    log:
        os.path.join(LOG_DIR, "workflow", "download_default_kraken2_db.log")
    benchmark:
        os.path.join(OUT_DIR, "benchmarks", "download_default_kraken2_db", "k2_viral.tsv")
    conda:
        DOWNLOAD
    shell:
        """
        mkdir -p "$(dirname "{log}")"
        exec > "{log}" 2>&1
        if [ -L "{params.db_dir}" ] || [ -f "{params.db_dir}" ]; then
            if [ ! -d "{params.db_dir}" ]; then
                echo "[INFO] Removing broken symlink or file at {params.db_dir}..."
                rm -f "{params.db_dir}"
            fi
        fi
        mkdir -p "{params.db_dir}"
        echo "[INFO] No external Kraken2 DB found. Downloading base RefSeq Viral database (k2_viral_20240112.tar.gz)..."
        rm -f "{params.tar_gz}"
        if command -v wget >/dev/null 2>&1; then
            wget "{params.url}" -O "{params.tar_gz}"
        else
            curl -sSL "{params.url}" -o "{params.tar_gz}"
        fi
        echo "[INFO] Extracting Kraken2 RefSeq Viral database into {params.db_dir}..."
        tar --no-same-owner --no-same-permissions -zxvf "{params.tar_gz}" -C "{params.db_dir}"
        rm -f "{params.tar_gz}"
        echo "[INFO] Base Kraken2 RefSeq Viral database ready in {params.db_dir}"
        """


rule kraken2_contigs:
    # Permite ferramentas normais E ferramentas virtuais de k-mer (ex: spades_k33)
    wildcard_constraints:
        tool = r"spades_k[\w]+|spades|megahit|flye|raven|medaka_flye|medaka_raven|pre_assembled"
    input:
        contigs = get_contigs_path,
        k2_db = _KRAKEN2_CORE_FILES
    output:
        report = os.path.join(OUT_DIR, "{sample}", "kraken2_{tool}", "{sample}_kraken2_contigs_{tool}_report.txt"),
        out    = os.path.join(OUT_DIR, "{sample}", "kraken2_{tool}", "{sample}_kraken2_contigs_{tool}_output.txt")
    params:
        db = _os.path.join(_PROJECT_ROOT, config["resources"]["kraken2"][0]),        
        confidence = config["kraken2"]["confidence"]
    log:
        os.path.join(LOG_DIR, "{sample}", "{sample}_kraken2_contigs_{tool}.log")
    benchmark:
        os.path.join(OUT_DIR, "benchmarks", "kraken2_contigs", "{sample}_kraken2_contigs_{tool}.tsv")
    conda:
        KRAKEN2 # Verifique se o nome do env está correto no seu config
    shell:
        """
        # O input.contigs agora é garantido ser um arquivo único por job
        kraken2 \
        --db {params.db} \
        --confidence {params.confidence} \
        --report {output.report} \
        --threads {resources.threads} \
        --memory-mapping \
        --output {output.out} \
        {input.contigs} \
        > {log} 2>&1
        """

rule kraken_biom_contig:
    # A mesma constraint deve ser aplicada aqui para o Snakemake "linkar" as regras
    wildcard_constraints:
        tool = r"spades_k[\w]+|spades|megahit|flye|raven|medaka_flye|medaka_raven|pre_assembled"
    input:
        report = os.path.join(OUT_DIR, "{sample}", "kraken2_{tool}", "{sample}_kraken2_contigs_{tool}_report.txt"),
    output:
        biom = os.path.join(OUT_DIR, "{sample}", "kraken2_{tool}", "{sample}_kraken_biom_contig_{tool}.txt"),
    log:
        os.path.join(LOG_DIR, "{sample}", "{sample}_kraken_biom_contig_{tool}.log")
    benchmark:
        os.path.join(OUT_DIR, "benchmarks", "kraken_biom_contig", "{sample}_kraken_biom_contig_{tool}.tsv")
    params:
        maximun = config["kraken_biom"]["max"],
        minimum = config["kraken_biom"]["min"],
        out_format = config["kraken_biom"]["format"]
    conda:
        KRAKEN2_BIOM 
    shell:
        """
        kraken-biom \
        {input.report} \
        --max {params.maximun} \
        --min {params.minimum} \
        --fmt {params.out_format} \
        -o {output.biom} \
        > {log} 2>&1
        """

rule kraken2_reads_paired:
    input:
     r1 = get_denovo_r1,
     r2 = get_denovo_r2,
     extra = get_denovo_unpaired,
     k2_db = _KRAKEN2_CORE_FILES
    output:
      report = os.path.join(OUT_DIR, "{sample}", "kraken2_reads", "{sample}_kraken2_reads_paired_report.txt"),
      out = os.path.join(OUT_DIR, "{sample}", "kraken2_reads", "{sample}_kraken2_reads_paired_output.txt")
    params:
      db = _os.path.join(_PROJECT_ROOT, config["resources"]["kraken2"][0]),        
      confidence = config["kraken2"]["confidence"],
      int_file = os.path.join(OUT_DIR, "{sample}", "kraken2_reads", "{sample}_kraken2_reads_paired_merged.fastq.gz")
    log:
      os.path.join(LOG_DIR, "{sample}", "{sample}_kraken2_reads_paired.log")
    benchmark:
      os.path.join(OUT_DIR, "benchmarks", "kraken2_reads_paired", "{sample}_kraken2_reads_paired.tsv")
    conda:
      KRAKEN2        
    shell:
      """
      mkdir -p $(dirname {params.int_file})
      cat {input.r1} {input.r2} {input.extra} > {params.int_file}
      kraken2 \
      --db {params.db} \
      --confidence {params.confidence} \
      --report {output.report} \
      --threads {resources.threads} \
      --memory-mapping \
      --output {output.out} \
      {params.int_file} \
      > {log} 2>&1
      rm -f {params.int_file}
      """

rule kraken2_reads_unpaired:
    input:   
      r1 = get_denovo_unpaired,
      k2_db = _KRAKEN2_CORE_FILES
    output:
      report = os.path.join(OUT_DIR, "{sample}", "kraken2_reads", "{sample}_kraken2_reads_unpaired_report.txt"),
      out = os.path.join(OUT_DIR, "{sample}", "kraken2_reads", "{sample}_kraken2_reads_unpaired_output.txt")
    params:
      db = _os.path.join(_PROJECT_ROOT, config["resources"]["kraken2"][0]),        
      confidence = config["kraken2"]["confidence"]
    log:
      os.path.join(LOG_DIR, "{sample}", "{sample}_kraken2_reads_unpaired.log")
    benchmark:
      os.path.join(OUT_DIR, "benchmarks", "kraken2_reads_unpaired", "{sample}_kraken2_reads_unpaired.tsv")
    conda:
      KRAKEN2        
    shell:
      """
      kraken2 \
      --db {params.db} \
      --confidence {params.confidence} \
      --report {output.report} \
      --threads {resources.threads} \
      --memory-mapping \
      --output {output.out} \
      {input.r1} \
      > {log} 2>&1
      """

rule kraken_biom_reads:
    input:
      report = os.path.join(OUT_DIR, "{sample}", "kraken2_reads", "{sample}_kraken2_reads_{paired}_report.txt")
    output:
      biom = os.path.join(OUT_DIR, "{sample}", "kraken2_reads", "{sample}_kraken_biom_reads_{paired}.txt"),
    params:
      maximun = config["kraken_biom"]["max"],
      minimum = config["kraken_biom"]["min"],
      out_format = config["kraken_biom"]["format"]
    log:
        os.path.join(LOG_DIR, "{sample}", "{sample}_kraken_biom_reads_{paired}.log")
    benchmark:
        os.path.join(OUT_DIR, "benchmarks", "kraken_biom_reads", "{sample}_kraken_biom_reads_{paired}.tsv")
    conda:
        KRAKEN2_BIOM
    shell:
        """
        kraken-biom \
        {input.report} \
        --max {params.maximun} \
        --min {params.minimum} \
        --fmt {params.out_format} \
        -o {output.biom} \
        > {log} 2>&1
        """
        
rule kraken_biom_merge_all:
    input:
        reports = get_all_kraken_reports
    output:
        biom = os.path.join(OUT_DIR, "kraken2_all", "all_samples_kraken_biom_merge_all.biom")
    params:
        fmt = "json"
    conda:
        R_RAREFACTION
    log:
        os.path.join(LOG_DIR, "workflow", "kraken_biom_merge_all.log")
    benchmark:
        os.path.join(OUT_DIR, "benchmarks", "kraken_biom_merge_all", "all_samples_kraken_biom_merge_all.tsv")
    shell:
        """
        kraken-biom {input.reports} --fmt {params.fmt} -o {output.biom} > {log} 2>&1
        """

# Rule 2: Run the R Script
rule kraken_rarefaction_plot:
    input:
        biom = os.path.join(OUT_DIR, "kraken2_all", "all_samples_kraken_biom_merge_all.biom")
    output:
        pdf = os.path.join(OUT_DIR, "kraken2_all", "all_samples_kraken_rarefaction_plot.pdf"),
        table = os.path.join(OUT_DIR, "kraken2_all", "all_samples_kraken_rarefaction_plot_otu_table.tab")
    conda:
        R_RAREFACTION
    log:
        os.path.join(LOG_DIR, "workflow", "kraken_rarefaction_plot.log")
    script:
        "../scripts/plot_rarefaction.R"