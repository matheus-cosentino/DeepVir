###################################################################################
#                       workflow/rules/diamond.smk                                #
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

rule download_default_diamond_db:
  """
  Automatically downloads and builds the base E. coli K-12 DIAMOND database
  (from the First Time Local Usage tutorial) on first run when no external
  or pre-existing DIAMOND database is present in resources/diamond/.
  """
  output:
    dmnd = os.path.join(RESOURCES_DIR, "diamond", "ecoli_db.dmnd")
  params:
    db_dir = os.path.join(RESOURCES_DIR, "diamond"),
    db_prefix = os.path.join(RESOURCES_DIR, "diamond", "ecoli_db"),
    raw_fasta = os.path.join(RESOURCES_DIR, "diamond", "ecoli_proteome.fasta"),
    fmt_fasta = os.path.join(RESOURCES_DIR, "diamond", "ecoli_formatted.fasta"),
    url = "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=(organism_id:83333)"
  log:
    os.path.join(LOG_DIR, "workflow", "download_default_diamond_db.log")
  benchmark:
    os.path.join(OUT_DIR, "benchmarks", "download_default_diamond_db", "ecoli_db.tsv")
  conda:
    DIAMOND
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
    echo "[INFO] No external DIAMOND DB found. Downloading base E. coli K-12 proteome from UniProt..."
    if command -v wget >/dev/null 2>&1; then
        wget "{params.url}" -O "{params.raw_fasta}"
    else
        curl -sSL "{params.url}" -o "{params.raw_fasta}"
    fi
    echo "[INFO] Formatting FASTA headers to keep protein accession IDs..."
    awk -F'|' '/^>/ {{print ">"$2}} !/^>/ {{print $0}}' "{params.raw_fasta}" > "{params.fmt_fasta}"
    echo "[INFO] Building base DIAMOND database: {output.dmnd}"
    diamond makedb --in "{params.fmt_fasta}" -d "{params.db_prefix}"
    rm -f "{params.raw_fasta}" "{params.fmt_fasta}"
    echo "[INFO] Base E. coli K-12 DIAMOND database ready at {output.dmnd}"
    """

rule filter_contigs_diamond:
  wildcard_constraints:
    tool = r"spades_k[\w]+|spades|megahit|flye|raven|medaka_flye|medaka_raven|pre_assembled"
  input:
    contigs = get_contigs_path
  output:
    filtered = temp(os.path.join(OUT_DIR, "{sample}", "diamond_{tool}", "{sample}_filter_contigs_diamond_{tool}.fasta"))
  params:
    min_len = config["diamond"]["min_contig_len"]
  log:
    os.path.join(LOG_DIR, "{sample}", "{sample}_filter_contigs_diamond_{tool}.log")
  benchmark:
    os.path.join(OUT_DIR, "benchmarks", "filter_contigs_diamond", "{sample}_filter_contigs_diamond_{tool}.tsv")
  conda:
    DIAMOND

  shell:
    """
    seqkit seq -m {params.min_len} {input.contigs} > {output.filtered} 2> {log}
    """

rule diamond_blastx_contigs:
  wildcard_constraints:
    tool = r"spades_k[\w]+|spades|megahit|flye|raven|medaka_flye|medaka_raven|pre_assembled"
  input:
    contigs = rules.filter_contigs_diamond.output.filtered,
    db_files = get_diamond_db_input
  output:
    daa = os.path.join(OUT_DIR, "{sample}", "diamond_{tool}", "{sample}_diamond_blastx_contigs_{tool}.daa"),
  params:
    db = get_diamond_db_name,
    db_dir = lambda wildcards: os.path.join(RESOURCES_DIR, "diamond"),
    max_target_seqs = config["diamond"]["max_target_seqs"],
    sensitivity = config["diamond"]["sensitivity"],
    block_size   = config["diamond"]["block_size"],
    index_chunks = config["diamond"]["index_chunks"]
  log:
    os.path.join(LOG_DIR, "{sample}", "{sample}_diamond_blastx_contigs_{tool}.log")
  benchmark:
    os.path.join(OUT_DIR, "benchmarks", "diamond_blastx_contigs", "{sample}_diamond_blastx_contigs_{tool}.tsv")
  conda:
    DIAMOND
  shadow:
    "minimal"
  shell:
    """
    exec > {log} 2>&1    
    echo "[INFO] Starting Diamond BlastX for Contigs (outfmt 100 / DAA)..."
    echo "[INFO] Input: {input.contigs}"
    
    DB_PATH="{params.db_dir}/{params.db}"
    
    echo "[INFO] Diamond Version:"
    diamond --version
    
    if [ ! -f "$DB_PATH" ] && [ ! -f "$DB_PATH.dmnd" ]; then
        echo "[ERROR] Diamond database not found at: $DB_PATH"
        exit 1
    fi

    # Garante que o TMPDIR exportado pelo DeepVir.sh ou o do sistema será usado
    # Fallback para o shadow dir do Snakemake, ou /tmp local do nó caso a variável falhe
    DIAMOND_TMP="${{TMPDIR:-${{SNAKEMAKE_SHADOW_DIR:-/tmp}}}}"
    mkdir -p "$DIAMOND_TMP"
    
    echo "[INFO] DB: $DB_PATH"
    echo "[INFO] Using Temporary Directory: $DIAMOND_TMP"
    
    diamond blastx \
      --query {input.contigs} \
      --db "$DB_PATH" \
      --out {output.daa} \
      --threads {resources.threads} \
      --outfmt 100 \
      --max-target-seqs {params.max_target_seqs} \
      {params.sensitivity} \
      --block-size {params.block_size} \
      --index-chunks {params.index_chunks} \
      --tmpdir "$DIAMOND_TMP"
    """

rule diamond_view_contigs:
  wildcard_constraints:
    tool = r"spades_k[\w]+|spades|megahit|flye|raven|medaka_flye|medaka_raven|pre_assembled"
  input:
    daa = rules.diamond_blastx_contigs.output.daa
  output:
    hits = os.path.join(OUT_DIR, "{sample}", "diamond_{tool}", "{sample}_diamond_view_contigs_{tool}.txt"),
  params:
    outfmt = config["diamond"]["outfmt"],
  log:
    os.path.join(LOG_DIR, "{sample}", "{sample}_diamond_view_contigs_{tool}.log")
  benchmark:
    os.path.join(OUT_DIR, "benchmarks", "diamond_view_contigs", "{sample}_diamond_view_contigs_{tool}.tsv")
  conda:
    DIAMOND
  shell:
    """
    exec > {log} 2>&1
    echo "[INFO] Converting DAA to tabular format (outfmt {params.outfmt})..."
    echo "[INFO] Input DAA: {input.daa}"

    diamond view \
      --daa {input.daa} \
      --outfmt {params.outfmt} \
      --out {output.hits} \
      --threads {resources.threads}

    echo "[INFO] Done. Output: {output.hits}"
    """


rule diamond_blastx_reads:
  input:
    r1 = get_denovo_r1,
    r2 = get_denovo_r2,
    extra = get_denovo_unpaired,
    db_files = get_diamond_db_input
  output:
    daa = os.path.join(OUT_DIR, "{sample}", "diamond_reads", "{sample}_diamond_blastx_reads.daa"),
  params:
    db = get_diamond_db_name,
    db_dir = lambda wildcards: os.path.join(RESOURCES_DIR, "diamond"),
    max_target_seqs = config["diamond"]["max_target_seqs"],
    sensitivity   = config["diamond"]["sensitivity"],
    block_size    = config["diamond"]["block_size"],
    index_chunks  = config["diamond"]["index_chunks"]
  log:
    os.path.join(LOG_DIR, "{sample}", "{sample}_diamond_blastx_reads.log")
  benchmark:
    os.path.join(OUT_DIR, "benchmarks", "diamond_blastx_reads", "{sample}_diamond_blastx_reads.tsv")
  conda:
    DIAMOND
  shadow:
    "shallow"
  shell:
    """
    exec > {log} 2>&1    
    echo "[INFO] Starting Diamond BlastX for Reads (outfmt 100 / DAA)..."
    echo "[INFO] Input R1: {input.r1}"
    echo "[INFO] Input R2: {input.r2}"
    echo "[INFO] Input Extra: {input.extra}"
    
    DB_PATH="{params.db_dir}/{params.db}"
    
    echo "[INFO] Diamond Version:"
    diamond --version
    
    if [ ! -f "$DB_PATH" ] && [ ! -f "$DB_PATH.dmnd" ]; then
        echo "[ERROR] Diamond database not found at: $DB_PATH"
        exit 1
    fi

    # Garante que o TMPDIR exportado pelo DeepVir.sh ou o do sistema será usado
    # Fallback para o shadow dir do Snakemake, ou /tmp local do nó caso a variável falhe
    DIAMOND_TMP="${{TMPDIR:-${{SNAKEMAKE_SHADOW_DIR:-/tmp}}}}"
    mkdir -p "$DIAMOND_TMP"

    echo "[INFO] DB: $DB_PATH"
    echo "[INFO] Using Temporary Directory: $DIAMOND_TMP"
    
    # Diamond supports a maximum of 2 query files in blastx mode, so we concatenate the inputs.
    # To avoid issues with concatenated gzip streams (multiple members), we decompress them on the fly.
    QUERY_FILE="$DIAMOND_TMP/merged_queries.fastq"
    if [[ "{input.r1}" == *.gz ]]; then
        gzip -cd {input.r1} {input.r2} {input.extra} > "$QUERY_FILE"
    else
        cat {input.r1} {input.r2} {input.extra} > "$QUERY_FILE"
    fi
    
    diamond blastx \
      --query "$QUERY_FILE" \
      --db "$DB_PATH" \
      --out {output.daa} \
      --threads {resources.threads} \
      --outfmt 100 \
      --max-target-seqs {params.max_target_seqs} \
      {params.sensitivity} \
      --block-size {params.block_size} \
      --index-chunks {params.index_chunks} \
      --tmpdir "$DIAMOND_TMP"
    """

rule diamond_view_reads:
  input:
    daa = rules.diamond_blastx_reads.output.daa
  output:
    hits = os.path.join(OUT_DIR, "{sample}", "diamond_reads", "{sample}_diamond_view_reads.txt"),
  params:
    outfmt = config["diamond"]["outfmt"],
  log:
    os.path.join(LOG_DIR, "{sample}", "{sample}_diamond_view_reads.log")
  benchmark:
    os.path.join(OUT_DIR, "benchmarks", "diamond_view_reads", "{sample}_diamond_view_reads.tsv")
  conda:
    DIAMOND
  shell:
    """
    exec > {log} 2>&1
    echo "[INFO] Converting DAA to tabular format (outfmt {params.outfmt})..."
    echo "[INFO] Input DAA: {input.daa}"

    diamond view \
      --daa {input.daa} \
      --outfmt {params.outfmt} \
      --out {output.hits} \
      --threads {resources.threads}

    echo "[INFO] Done. Output: {output.hits}"
    """