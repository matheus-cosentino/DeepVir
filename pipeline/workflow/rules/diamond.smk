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

#########################################################################
# --- Default DIAMOND DB (RVDB protein FASTA): download > unpack > cd-hit > makedb --- #
#########################################################################
# Runs on first use only, when no external/pre-existing DIAMOND DB is present
# in resources/diamond/. RVDB headers: >acc|p_bank|p_acc|n_bank|n_acc|descr[sp]

rule download_rvdb_prot_fasta:
  """Step 1/4: download the RVDB protein FASTA (xz-compressed)."""
  output:
    xz = temp(os.path.join(RESOURCES_DIR, "diamond", "U-RVDBv32.0-prot.fasta.xz"))
  params:
    db_dir = os.path.join(RESOURCES_DIR, "diamond"),
    url = "https://rvdb-prot.pasteur.fr/files/U-RVDBv32.0-prot.fasta.xz"
  log:
    os.path.join(LOG_DIR, "workflow", "download_rvdb_prot_fasta.log")
  benchmark:
    os.path.join(OUT_DIR, "benchmarks", "download_rvdb_prot_fasta", "download_rvdb_prot_fasta.tsv")
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
    echo "[INFO] No external DIAMOND DB found. Downloading default RVDB protein FASTA..."
    rm -f "{output.xz}"
    if command -v wget >/dev/null 2>&1; then
        wget "{params.url}" -O "{output.xz}"
    else
        curl -sSL "{params.url}" -o "{output.xz}"
    fi
    echo "[INFO] Download complete: {output.xz}"
    """

rule unpack_rvdb_prot_fasta:
  """Step 2/4: decompress the RVDB protein FASTA (.xz -> .fasta)."""
  input:
    xz = rules.download_rvdb_prot_fasta.output.xz
  output:
    fasta = temp(os.path.join(RESOURCES_DIR, "diamond", "U-RVDBv32.0-prot.fasta"))
  log:
    os.path.join(LOG_DIR, "workflow", "unpack_rvdb_prot_fasta.log")
  benchmark:
    os.path.join(OUT_DIR, "benchmarks", "unpack_rvdb_prot_fasta", "unpack_rvdb_prot_fasta.tsv")
  conda:
    DOWNLOAD
  shell:
    """
    mkdir -p "$(dirname "{log}")"
    exec > "{log}" 2>&1
    echo "[INFO] Decompressing {input.xz}..."
    xz -dc "{input.xz}" > "{output.fasta}"
    echo "[INFO] Unpacked FASTA: {output.fasta}"
    """

rule cd_hit_rvdb_prot_fasta:
  """Step 3/4: cluster RVDB protein FASTA at 90% identity using cd-hit to reduce redundancy."""
  input:
    fasta = rules.unpack_rvdb_prot_fasta.output.fasta
  output:
    fasta = temp(os.path.join(RESOURCES_DIR, "diamond", "U-RVDBv32.0-prot_nr90.fasta")),
    clstr = temp(os.path.join(RESOURCES_DIR, "diamond", "U-RVDBv32.0-prot_nr90.fasta.clstr"))
  params:
    c = 0.90,      # 90% sequence identity threshold
#    aS = 0.90,     # alignment coverage for the shorter sequence
    d = 0          # keep full fasta header description (>acc|p_bank|p_acc|n_bank|n_acc|descr[sp])
  log:
    os.path.join(LOG_DIR, "workflow", "cd_hit_rvdb_prot_fasta.log")
  benchmark:
    os.path.join(OUT_DIR, "benchmarks", "cd_hit_rvdb_prot_fasta", "cd_hit_rvdb_prot_fasta.tsv")
  conda:
    PALM
  shell:
    """
    mkdir -p "$(dirname "{log}")"
    exec > "{log}" 2>&1
    echo "[INFO] Running cd-hit on {input.fasta} (identity: {params.c})..."
    cd-hit -i "{input.fasta}" -o "{output.fasta}" -c {params.c} -M {resources.mem_mb} -d {params.d} -T {resources.threads}
    echo "[INFO] Redundancy reduction complete: {output.fasta}"
    """

rule makedb_rvdb_prot_diamond:
  """Step 4/4: build the default DIAMOND database (rvdb_prot.dmnd) from non-redundant FASTA."""
  input:
    fasta = rules.cd_hit_rvdb_prot_fasta.output.fasta
  output:
    dmnd = os.path.join(RESOURCES_DIR, "diamond", "rvdb_prot.dmnd")
  params:
    db_prefix = os.path.join(RESOURCES_DIR, "diamond", "rvdb_prot")
  log:
    os.path.join(LOG_DIR, "workflow", "makedb_rvdb_prot_diamond.log")
  benchmark:
    os.path.join(OUT_DIR, "benchmarks", "makedb_rvdb_prot_diamond", "makedb_rvdb_prot_diamond.tsv")
  conda:
    DIAMOND
  shell:
    """
    mkdir -p "$(dirname "{log}")"
    exec > "{log}" 2>&1
    echo "[INFO] Building default DIAMOND database: {output.dmnd}"
    diamond makedb --in "{input.fasta}" -d "{params.db_prefix}" --threads {threads}
    echo "[INFO] Default RVDB protein DIAMOND database ready at {output.dmnd}"
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