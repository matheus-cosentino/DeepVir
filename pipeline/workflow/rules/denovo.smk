###################################################################################
#                       workflow/rules/denovo.smk                                 #
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

###############################
# --- 1. Spades Assembly --- #
##############################

rule spades:
    input:
        r1 = get_denovo_r1,
        r2 = get_denovo_r2,
        extra = get_denovo_unpaired
    output:
        contigs = os.path.join(OUT_DIR, "{sample}", "spades", "kmer_{kmer_val}", "{sample}_spades_k{kmer_val}.fasta")
    log:
        os.path.join(LOG_DIR, "{sample}", "{sample}_spades_k{kmer_val}.log")
    benchmark:
        os.path.join(OUT_DIR, "benchmarks", "spades", "{sample}_spades_k{kmer_val}.tsv")
    params:
        outdir = lambda w, output: os.path.dirname(output.contigs),
        input_args = get_spades_params,
        extra = config["spades"]["algorithm"],
        #kmer = lambda wildcards: wildcards.kmer_val
        kmer = lambda wildcards: wildcards.kmer_val.replace("_", ","),
        mem_gb = lambda wildcards, resources: max(1, int(resources.mem_mb // 1024))
    conda:
        SPADES_ENV
    shadow:
        "minimal"
    shell:
        """
        export OMP_NUM_THREADS={resources.threads}
        spades.py \
            {params.extra} \
            --threads {resources.threads} \
            --memory {params.mem_gb} \
            {params.input_args} \
            -o {params.outdir} \
            --tmp-dir {resources.tmpdir} \
            -k {params.kmer} \
            --only-assembler \
            > {log} 2>&1
        if [ -f "{params.outdir}/contigs.fasta" ]; then
            cp -f "{params.outdir}/contigs.fasta" "{output.contigs}"
        else
            touch "{output.contigs}"
        fi
        """

###############################
# --- 2. Megahit Assembly --- #
##############################

rule megahit:
    input:
        r1 = get_denovo_r1,
        r2 = get_denovo_r2,
        extra = get_denovo_unpaired
    output:
        contigs = os.path.join(OUT_DIR, "{sample}", "megahit", "{sample}_megahit.fasta")
    log:
        os.path.join(LOG_DIR, "{sample}", "{sample}_megahit.log")
    benchmark:
        os.path.join(OUT_DIR, "benchmarks", "megahit", "{sample}_megahit.tsv")
    params:
        outdir = lambda w, output: os.path.dirname(output.contigs),
        input_args = get_megahit_params,
    conda:
        MEGAHIT_ENV
    shell:
        """
        tmp_megahit_out=$(mktemp -d -p "{resources.tmpdir}" megahit_{wildcards.sample}_XXXXXX)
        mkdir -p "{params.outdir}"
        
        megahit \
          {params.input_args} \
          --tmp-dir "{resources.tmpdir}" \
          -o "$tmp_megahit_out" \
          -t {resources.threads} \
          -f \
          > {log} 2>&1
        
        cp -rf "$tmp_megahit_out"/* "{params.outdir}/"
        rm -rf "$tmp_megahit_out"
        if [ -f "{params.outdir}/final.contigs.fa" ]; then
            cp -f "{params.outdir}/final.contigs.fa" "{output.contigs}"
        else
            touch "{output.contigs}"
        fi
        """


#############################
# --- 3. Flye Assembly --- #
############################

rule flye:  
    """
    WARNING: Flye is designed for Long Reads (Nanopore/PacBio).
    Do not use Illumina reads here.
    """
    input:
        reads = get_ONP_input
    output:
        contigs = os.path.join(OUT_DIR, "{sample}", "flye", "{sample}_flye.fasta")
    log:
        os.path.join(LOG_DIR, "{sample}", "{sample}_flye.log")
    benchmark:
        os.path.join(OUT_DIR, "benchmarks", "flye", "{sample}_flye.tsv")
    params:
        outdir = lambda w, output: os.path.dirname(output.contigs),
        extra = config["flye"]["type"]
    conda:
        FLYE_ENV
    shell:
        """
        # Flye fails if output directory exists
        rm -rf {params.outdir}

        # Run commands 

        flye --{params.extra} {input.reads} \
         --out-dir {params.outdir} \
         --threads {resources.threads} \
         --meta \
         > {log} 2>&1
        if [ -f "{params.outdir}/assembly.fasta" ]; then
            cp -f "{params.outdir}/assembly.fasta" "{output.contigs}"
        else
            touch "{output.contigs}"
        fi
        """

#############################
# --- 4. Raven Assembly --- #
############################

rule raven:  
    """
    WARNING: Flye is designed for Long Reads (Nanopore/PacBio).
    Do not use Illumina reads here.
    """
    input:
        reads = get_ONP_input
    output:
        contigs = os.path.join(OUT_DIR, "{sample}", "raven", "{sample}_raven.fasta")
    log:
        os.path.join(LOG_DIR, "{sample}", "{sample}_raven.log")
    benchmark:
        os.path.join(OUT_DIR, "benchmarks", "raven", "{sample}_raven.tsv")
    conda:
        RAVEN_ENV
    shell:
        """
        # Run commands 
        raven --threads {resources.threads} {input.reads} > {output.contigs} 2> {log}
        """

###############################
# --- 5. Medaka Correcton --- #
###############################

rule medaka_polish:
    input:
        reads = get_ONP_input,
        draft = os.path.join(OUT_DIR, "{sample}", "{assembler}", "{sample}_{assembler}.fasta")
    output:
        consensus = os.path.join(OUT_DIR, "{sample}", "medaka_{assembler}", "{sample}_medaka_polish_{assembler}.fasta")
    params:
        # Defina o modelo aqui ou no config.yaml
        model = config["medaka_model"],
        outdir = lambda w, output: os.path.dirname(output.consensus)
    log:
        os.path.join(LOG_DIR, "{sample}", "{sample}_medaka_polish_{assembler}.log")
    benchmark:
        os.path.join(OUT_DIR, "benchmarks", "medaka_polish", "{sample}_medaka_polish_{assembler}.tsv")
    conda:
        MEDAKA_ENV
    shell:
        """
        rm -rf {params.outdir}
        medaka_consensus -i {input.reads} \
                         -d {input.draft} \
                         -o {params.outdir} \
                         -t {resources.threads} \
                         -m {params.model} \
        > {log} 2>&1
        if [ -f "{params.outdir}/consensus.fasta" ]; then
            cp -f "{params.outdir}/consensus.fasta" "{output.consensus}"
        else
            touch "{output.consensus}"
        fi
        """