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

rule kraken2_contigs:
    # Permite ferramentas normais E ferramentas virtuais de k-mer (ex: spades_k33)
    wildcard_constraints:
        tool = r"spades_k[\w]+|spades|megahit|flye|raven|medaka_flye|medaka_raven|pre_assembled"
    input:
        contigs = get_contigs_path
    output:
        report = os.path.join(OUT_DIR, "{sample}", "kraken2_{tool}", "{sample}_{tool}_contig_report.txt"),
        out    = os.path.join(OUT_DIR, "{sample}", "kraken2_{tool}", "{sample}_{tool}_contig_output.txt")
    params:
        db = _os.path.join(_PROJECT_ROOT, config["resources"]["kraken2"][0]),        
        confidence = config["kraken2"]["confidence"]
    log:
        os.path.join(LOG_DIR, "{sample}", "kraken2_contigs_{tool}_{sample}.log")
    benchmark:
        os.path.join(OUT_DIR, "benchmarks", "kraken2_contigs", "{sample}_{tool}.tsv")
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
        report = os.path.join(OUT_DIR, "{sample}", "kraken2_{tool}", "{sample}_{tool}_contig_report.txt"),
    output:
        biom = os.path.join(OUT_DIR, "{sample}", "kraken2_{tool}", "{sample}_{tool}_contig_biom.txt"),
    log:
        os.path.join(LOG_DIR, "{sample}", "kraken2_contigs_{tool}_{sample}_biom.log")
    benchmark:
        os.path.join(OUT_DIR, "benchmarks", "kraken_biom_contig", "{sample}_{tool}.tsv")
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
     r1 = os.path.join(OUT_DIR, "{sample}", "fastp", "{sample}_1.fastq.gz"),
     r2 = os.path.join(OUT_DIR, "{sample}", "fastp" , "{sample}_2.fastq.gz"),
     extra = os.path.join(OUT_DIR, "{sample}", "fastp", "{sample}_orphans.fastq.gz")
    output:
      report= os.path.join(OUT_DIR, "{sample}", "kraken2_reads", "{sample}_paired_reads_report.txt"),
      out=os.path.join(OUT_DIR, "{sample}", "kraken2_reads", "{sample}_paired_reads_output.txt")
    params:
      #db=f"{workflow.basedir}/{config["resources"]["kraken2"]}",
      db=_os.path.join(_PROJECT_ROOT, config["resources"]["kraken2"][0]),        
      confidence=config["kraken2"]["confidence"],
      int_file=temp(os.path.join(OUT_DIR, "{sample}", "fastp", "{sample}.fastq.gz"))
    log:
      os.path.join(LOG_DIR, "{sample}", "kraken2_paired_reads_{sample}.log")
    benchmark:
      os.path.join(OUT_DIR, "benchmarks", "kraken2_reads_paired", "{sample}.tsv")
    conda:
      KRAKEN2        
    shell:
      """
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
      """

rule kraken2_reads_unpaired:
    input:   
      r1= os.path.join(OUT_DIR, "{sample}", "fastp", "{sample}_unp.fastq.gz")
    output:
      report= os.path.join(OUT_DIR, "{sample}", "kraken2_reads", "{sample}_unpaired_reads_report.txt"),
      out=os.path.join(OUT_DIR, "{sample}", "kraken2_reads", "{sample}_unpaired_reads_output.txt")
    params:
      #db=f"{workflow.basedir}/{config["resources"]["kraken2"]}",
      db=_os.path.join(_PROJECT_ROOT, config["resources"]["kraken2"][0]),        
      confidence=config["kraken2"]["confidence"]
    log:
      os.path.join(LOG_DIR, "{sample}", "kraken2_paired_reads_{sample}.log")
    benchmark:
      os.path.join(OUT_DIR, "benchmarks", "kraken2_reads_unpaired", "{sample}.tsv")
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
      report= os.path.join(OUT_DIR, "{sample}", "kraken2_reads", "{sample}_{paired}_reads_report.txt")
    output:
      biom= os.path.join(OUT_DIR, "{sample}", "kraken2_reads", "{sample}_{paired}_reads_biom.txt"),
    params:
      maximun=config["kraken_biom"]["max"],
      minimum=config["kraken_biom"]["min"],
      out_format=config["kraken_biom"]["format"]
    log:
        os.path.join(LOG_DIR, "{sample}", "kraken2_biom_{paired}_reads_{sample}.log")
    benchmark:
        os.path.join(OUT_DIR, "benchmarks", "kraken_biom_reads", "{sample}_{paired}.tsv")
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
        biom = os.path.join(OUT_DIR, "kraken2_all", "all_samples.biom")
    params:
        fmt = "json"
    conda:
        R_RAREFACTION
    log:
        os.path.join(LOG_DIR, "workflow", "kraken_biom_merge_all.log")
    benchmark:
        os.path.join(OUT_DIR, "benchmarks", "kraken_biom_merge_all", "all_samples.tsv")
    shell:
        """
        kraken-biom {input.reports} --fmt {params.fmt} -o {output.biom} > {log} 2>&1
        """

# Rule 2: Run the R Script
rule kraken_rarefaction_plot:
    input:
        biom = os.path.join(OUT_DIR, "kraken2_all", "all_samples.biom")
    output:
        pdf = os.path.join(OUT_DIR, "kraken2_all", "Rarefaction_Curve.pdf"),
        table = os.path.join(OUT_DIR, "kraken2_all", "OTU_table.tab")
    conda:
        R_RAREFACTION
    log:
        os.path.join(LOG_DIR, "workflow", "kraken_rarefaction_plot.log")
    script:
        "../scripts/plot_rarefaction.R"