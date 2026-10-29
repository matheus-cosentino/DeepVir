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
        # Define all possible output files explicitly and statically.
        #r1= temp(os.path.join(OUT_DIR, "{sample}", "fastp", "{sample}_1.fastq.gz")),
        #r2= temp(os.path.join(OUT_DIR, "{sample}", "fastp" , "{sample}_2.fastq.gz")),
        r1= os.path.join(OUT_DIR, "{sample}", "fastp", "{sample}_1.fastq.gz"),
        r2= os.path.join(OUT_DIR, "{sample}", "fastp" , "{sample}_2.fastq.gz"),
        orphans=os.path.join(OUT_DIR, "{sample}", "fastp", "{sample}_orphans.fastq.gz"),
       #orphans=temp(os.path.join(OUT_DIR, "{sample}", "fastp", "{sample}_orphans.fastq.gz")),
        html = os.path.join(OUT_DIR, "{sample}", "fastp", "{sample}_paired.html"),
        json= os.path.join(OUT_DIR, "{sample}", "fastp", "{sample}_paired.json")
    log:
        os.path.join(LOG_DIR, "{sample}", "{sample}_fastp_paired.log")
    benchmark:
        os.path.join(OUT_DIR, "benchmarks", "fastp_paired", "{sample}.tsv")
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
        #r1= temp(os.path.join(OUT_DIR, "{sample}", "fastp", "{sample}_unp.fastq.gz")),
        r1= os.path.join(OUT_DIR, "{sample}", "fastp", "{sample}_unp.fastq.gz"),
        html= os.path.join(OUT_DIR, "{sample}", "fastp", "{sample}_unp.html"),      
        json= os.path.join(OUT_DIR, "{sample}", "fastp", "{sample}_unp.json")
    log:
        os.path.join(LOG_DIR, "{sample}", "{sample}_fastp_unpaired.log")
    benchmark:
        os.path.join(OUT_DIR, "benchmarks", "fastp_unpaired", "{sample}.tsv")
    params:
        length_required=lambda wc: config['fastp']['length_required'][0] if isinstance(config['fastp']['length_required'], (list, tuple)) else config['fastp']['length_required'],
        quality=lambda wc: config['fastp']['qualified_quality_phred'][0] if isinstance(config['fastp']['qualified_quality_phred'], (list, tuple)) else config['fastp']['qualified_quality_phred']
    conda:
        FASTP
    shell:
        """
        fastp \
            --in1 {input.reads} --out1 {output.r1} \
            --html {output.html} \
            --json {output.json} \
            --thread {resources.threads} \
            --length_required {params.length_required} \
            --qualified_quality_phred {params.quality} \
            --dedup --trim_poly_g --trim_poly_x \
            2> {log}
        """
