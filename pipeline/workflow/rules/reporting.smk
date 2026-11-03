###################################################################################
#                      workflow/rules/reportin.smk                                #
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

#######################################################################
# --- Statistics Generation for MultiQC Integration ---              #
#######################################################################

# --- Generate Assembly Statistics ---
rule generate_assembly_stats:
    """Generate contig statistics from assembly output"""
    wildcard_constraints:
        tool = r"spades_k[\w]+|spades|megahit|flye|raven|medaka_flye|medaka_raven|pre_assembled"
    input:
        fasta = get_contigs_path
    output:
        json = os.path.join(OUT_DIR, "{sample}", "stats", "{sample}_generate_assembly_stats_{tool}.json"),
        tsv = os.path.join(OUT_DIR, "{sample}", "stats", "{sample}_generate_assembly_stats_{tool}.tsv")
    log:
        os.path.join(LOG_DIR, "{sample}", "{sample}_generate_assembly_stats_{tool}.log")
    conda:
        CORE
    script:
        "../scripts/generate_assembly_stats.py"


# --- Generate Diamond Statistics ---
rule generate_diamond_stats:
    """Generate Diamond hit statistics for contigs"""
    wildcard_constraints:
        tool = r"spades_k[\w]+|spades|megahit|flye|raven|medaka_flye|medaka_raven|pre_assembled"
    input:
        diamond = os.path.join(OUT_DIR, "{sample}", "diamond_{tool}", "{sample}_diamond_view_contigs_{tool}.txt")
    output:
        json = os.path.join(OUT_DIR, "{sample}", "stats", "{sample}_generate_diamond_stats_{tool}.json"),
        tsv = os.path.join(OUT_DIR, "{sample}", "stats", "{sample}_generate_diamond_stats_{tool}.tsv")
    log:
        os.path.join(LOG_DIR, "{sample}", "{sample}_generate_diamond_stats_{tool}.log")
    conda:
        CORE
    script:
        "../scripts/generate_diamond_stats.py"


# --- Generate Diamond Statistics (Reads) ---
rule generate_diamond_reads_stats:
    """Generate Diamond hit statistics for reads"""
    input:
        diamond = os.path.join(OUT_DIR, "{sample}", "diamond_reads", "{sample}_diamond_view_reads.txt")
    output:
        json = os.path.join(OUT_DIR, "{sample}", "stats", "{sample}_generate_diamond_reads_stats.json"),
        tsv = os.path.join(OUT_DIR, "{sample}", "stats", "{sample}_generate_diamond_reads_stats.tsv")
    log:
        os.path.join(LOG_DIR, "{sample}", "{sample}_generate_diamond_reads_stats.log")
    conda:
        CORE
    script:
        "../scripts/generate_diamond_stats.py"


# --- Generate Kraken2 Statistics (Contigs) ---
rule generate_kraken2_stats:
    """Generate Kraken2 classification statistics for contigs"""
    wildcard_constraints:
        tool = r"spades_k[\w]+|spades|megahit|flye|raven|medaka_flye|medaka_raven|pre_assembled"
    input:
        report = os.path.join(OUT_DIR, "{sample}", "kraken2_{tool}", "{sample}_kraken2_contigs_{tool}_report.txt")
    output:
        json = os.path.join(OUT_DIR, "{sample}", "stats", "{sample}_generate_kraken2_stats_{tool}.json"),
        tsv = os.path.join(OUT_DIR, "{sample}", "stats", "{sample}_generate_kraken2_stats_{tool}.tsv")
    log:
        os.path.join(LOG_DIR, "{sample}", "{sample}_generate_kraken2_stats_{tool}.log")
    conda:
        CORE
    script:
        "../scripts/generate_kraken2_stats.py"


# --- Generate Kraken2 Statistics (Reads) ---
rule generate_kraken2_reads_stats:
    """Generate Kraken2 classification statistics for reads"""
    input:
        report = lambda wc: os.path.join(OUT_DIR, wc.sample, "kraken2_reads", 
                            f"{wc.sample}_kraken2_reads_{'paired' if (SAMPLE_META.get(wc.sample, {}).get('mode') == 'PAIRED' or (SAMPLE_META.get(wc.sample, {}).get('mode') == 'SRA' and len(SAMPLE_META.get(wc.sample, {}).get('files', [])) == 2)) else 'unpaired'}_report.txt")
    output:
        json = os.path.join(OUT_DIR, "{sample}", "stats", "{sample}_generate_kraken2_reads_stats.json"),
        tsv = os.path.join(OUT_DIR, "{sample}", "stats", "{sample}_generate_kraken2_reads_stats.tsv")
    log:
        os.path.join(LOG_DIR, "{sample}", "{sample}_generate_kraken2_reads_stats.log")
    conda:
        CORE
    script:
        "../scripts/generate_kraken2_stats.py"


# --- Generate Dark Matter Statistics ---
rule generate_darkmatter_stats:
    """Generate dark matter (no-hit contigs) statistics"""
    wildcard_constraints:
        tool = r"spades_k[\w]+|spades|megahit|flye|raven|medaka_flye|medaka_raven|pre_assembled"
    input:
        nohits = os.path.join(OUT_DIR, "{sample}", "darkmatter_{tool}", "{sample}_get_nohit_fasta_{tool}.fasta")
    output:
        json = os.path.join(OUT_DIR, "{sample}", "stats", "{sample}_generate_darkmatter_stats_{tool}.json"),
        tsv = os.path.join(OUT_DIR, "{sample}", "stats", "{sample}_generate_darkmatter_stats_{tool}.tsv")
    log:
        os.path.join(LOG_DIR, "{sample}", "{sample}_generate_darkmatter_stats_{tool}.log")
    conda:
        CORE
    script:
        "../scripts/generate_darkmatter_stats.py"

rule multiqc_aggregate:
    message:
        """
        Generate a multiqc report in HTML format all samples
        """
    conda:
        MULTIQC
    input:
        files = [f for s in SAMPLE for f in get_multiqc_inputs(sample=s)],
        final_outputs = get_final_outputs()
    output:
        report = os.path.join(OUT_DIR, "all_samples_multiqc_aggregate.html"),
        data_dir = directory(os.path.join(OUT_DIR, "all_samples_multiqc_aggregate_data"))
    params:
        files = [f for s in SAMPLE for f in get_multiqc_inputs(sample=s)],
        config = os.path.join(workflow.basedir, "multiqc_config.yaml"),
        extra = "--title 'DeepVir Aggregate Report'",
        outdir = lambda w, output: os.path.dirname(output.report)
    conda:
        MULTIQC    
    log:
        os.path.join(LOG_DIR, "workflow", "multiqc_aggregate.log")
    shell:
        """
        multiqc \
        --quiet \
        --export \
        --force \
        --config {params.config} \
        --outdir {params.outdir} \
        --filename $(basename {output.report}) \
        {params.extra} \
        {input.files} > {log} 2>&1 
        """

# --- Rule 2: Per-Sample Report ---
rule multiqc_sample:
    message:
        "Generate MultiQC report for sample: {wildcards.sample}"
    conda:
        MULTIQC
    input:
        # We pass the wildcards object to get inputs ONLY for this sample
        files = lambda wc: get_multiqc_inputs(wildcards=wc)
    output:
        report = os.path.join(OUT_DIR, "{sample}", "multiqc", "{sample}_multiqc_sample.html"),
        data_dir = directory(os.path.join(OUT_DIR, "{sample}", "multiqc", "{sample}_multiqc_sample_data"))
    params:
        config = os.path.join(workflow.basedir, "multiqc_config.yaml"),
        extra = "--title 'Report for {sample}'",
        outdir = lambda w, output: os.path.dirname(output.report)
    log:
        os.path.join(LOG_DIR, "{sample}", "{sample}_multiqc_sample.log")
    shell:
        """
        multiqc \
        --quiet \
        --export \
        --force \
        --config {params.config} \
        --outdir {params.outdir} \
        --filename $(basename {output.report}) \
        {params.extra} \
        {input.files} > {log} 2>&1 
        """