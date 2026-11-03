###################################################################################
#                        workflow/rules/RVDB.smk                                  #
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

rule get_rvdb:
  output:
    hmm = os.path.join(RVDB_DIR[0], "U-RVDBv31.0-prot.hmm.xz"),
    sql =  os.path.join(RVDB_DIR[0], "U-RVDBv31.0-prot-hmm.sqlite.xz")
  log:
    os.path.join(LOG_DIR, "workflow", "get_rvdb.log")
  conda:
    DOWNLOAD
  shell:
    """
    curl -L https://rvdb-prot.pasteur.fr/files/U-RVDBv31.0-prot.hmm.xz -o  {output.hmm} 2>> {log}
    curl -L https://rvdb-prot.pasteur.fr/files/U-RVDBv31.0-prot-hmm.sqlite.xz -o  {output.sql} 2>> {log}
    """

rule prepare_rvdb:
  input:
    hmm = os.path.join(RVDB_DIR[0], "U-RVDBv31.0-prot.hmm.xz"),
    sql = os.path.join(RVDB_DIR[0], "U-RVDBv31.0-prot-hmm.sqlite.xz")
  output:
    hmm = os.path.join(RVDB_DIR[0], "U-RVDBv31.0-prot.hmm"),
    sql = os.path.join(RVDB_DIR[0], "U-RVDBv31.0-prot-hmm.sqlite"),
    i1 = os.path.join(RVDB_DIR[0], "U-RVDBv31.0-prot.hmm.h3m"),
    i2 = os.path.join(RVDB_DIR[0], "U-RVDBv31.0-prot.hmm.h3i"),
    i3 = os.path.join(RVDB_DIR[0], "U-RVDBv31.0-prot.hmm.h3f"),
    i4 = os.path.join(RVDB_DIR[0], "U-RVDBv31.0-prot.hmm.h3p")
  log:
    os.path.join(LOG_DIR, "workflow", "prepare_rvdb.log")
  conda:
    HMMER
  shell:
    """
    xz -dc {input.hmm} > {output.hmm} 2>> {log}
    xz -dc {input.sql} > {output.sql} 2>> {log}
    hmmpress -f {output.hmm} >> {log} 2>&1
    """

rule rvdb_search:
  input:
    orfs = os.path.join(OUT_DIR, "{sample}", "darkmatter_{tool}", "{sample}_cd_hit_{tool}.fasta"),
    hmm = os.path.join(RVDB_DIR[0], "U-RVDBv31.0-prot.hmm"),
    indices = multiext(os.path.join(RVDB_DIR[0], "U-RVDBv31.0-prot.hmm"), ".h3m", ".h3i", ".h3f", ".h3p")
  output:
    result = os.path.join(OUT_DIR, "{sample}", "rvdb_{tool}", "{sample}_rvdb_search_{tool}.tsv")
  params:
    evalue=config['hmmscan']['evalue']
  log:
    os.path.join(LOG_DIR, "{sample}", "{sample}_rvdb_search_{tool}.log")
  conda:
    HMMER
  shell:
    """
    hmmscan --cpu {resources.threads} -E {params.evalue} --tblout {output.result} {input.hmm} {input.orfs} 2>> {log}
    """   

rule rvdb_summarize:
  input:
    tbl = os.path.join(OUT_DIR, "{sample}", "rvdb_{tool}", "{sample}_rvdb_search_{tool}.tsv"),
    sqlite = os.path.join(RVDB_DIR[0], "U-RVDBv31.0-prot-hmm.sqlite")
  output:
    csv = os.path.join(OUT_DIR, "{sample}", "rvdb_{tool}", "{sample}_rvdb_summarize_{tool}.tsv")
  log:
    os.path.join(LOG_DIR, "{sample}", "{sample}_rvdb_summarize_{tool}.log")
  conda:
    CORE
  script:
    "../scripts/rvdb_summarize.py"

rule rvdb_validate_structural:
  input:
    fasta = get_contigs_path,
    orfs  = os.path.join(OUT_DIR, "{sample}", "darkmatter_{tool}", "{sample}_cd_hit_{tool}.fasta"),
    tsv   = os.path.join(OUT_DIR, "{sample}", "rvdb_{tool}", "{sample}_rvdb_summarize_{tool}.tsv")
  output:
    raw_rvdb_orfs    = os.path.join(OUT_DIR, "{sample}", "rvdb_{tool}", "{sample}_rvdb_validate_structural_{tool}_orfs.fasta"),
    raw_rvdb_contigs = os.path.join(OUT_DIR, "{sample}", "rvdb_{tool}", "{sample}_rvdb_validate_structural_{tool}_contigs.fasta")
  log:
    os.path.join(LOG_DIR, "{sample}", "{sample}_rvdb_validate_structural_{tool}.log")
  conda:
    REPORT
  script:
    "../scripts/rvdb_validate_structural.py"

rule rvdb_validate_pol:
  input:
    fasta = get_contigs_path,
    orfs  = os.path.join(OUT_DIR, "{sample}", "darkmatter_{tool}", "{sample}_cd_hit_{tool}.fasta"),
    tsv   = os.path.join(OUT_DIR, "{sample}", "rvdb_{tool}", "{sample}_rvdb_summarize_{tool}.tsv")
  output:
    raw_rvdb_orfs    = os.path.join(OUT_DIR, "{sample}", "rvdb_{tool}", "{sample}_rvdb_validate_pol_{tool}_orfs.fasta"),
    raw_rvdb_contigs = os.path.join(OUT_DIR, "{sample}", "rvdb_{tool}", "{sample}_rvdb_validate_pol_{tool}_contigs.fasta")
  log:
    os.path.join(LOG_DIR, "{sample}", "{sample}_rvdb_validate_pol_{tool}.log")
  conda:
    REPORT
  script:
    "../scripts/rvdb_validate_pol.py"