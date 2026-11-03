###################################################################################
#                      workflow/rules/duskmatter.smk                              #
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


######################################################
# --- 1. Filter Contig Fasta With no Diamond Hit --- #
###################################################### 
rule get_nohit_fasta:
  """
    Get contigs within no hits against the nr Diamond database
  """
  # ADICIONAR ESTA CONSTRAINT
  wildcard_constraints:
    tool = r"spades_k[\w]+|spades|megahit|flye|raven|medaka_flye|medaka_raven|pre_assembled"
  input:
    fasta = get_contigs_path,
    diamond = os.path.join(OUT_DIR, "{sample}", "diamond_{tool}", "{sample}_append_lineage_{tool}.tsv")  
  output:
    nohits = os.path.join(OUT_DIR, "{sample}", "darkmatter_{tool}", "{sample}_get_nohit_fasta_{tool}.fasta")
  log:
    os.path.join(LOG_DIR, "{sample}", "{sample}_get_nohit_fasta_{tool}.log")
  benchmark:
    os.path.join(OUT_DIR, "benchmarks", "get_nohit_fasta", "{sample}_get_nohit_fasta_{tool}.tsv")
  conda:
    CORE
  script:
    "../scripts/filter_nohits.py"


#################################################
# --- 2. Find Orfs for the no Hits Contigs --- #
################################################ 
rule find_orfs:
  input: 
    fasta = os.path.join(OUT_DIR, "{sample}", "darkmatter_{tool}", "{sample}_get_nohit_fasta_{tool}.fasta")
  output:
    orfs = temp(os.path.join(OUT_DIR, "{sample}", "darkmatter_{tool}", "{sample}_find_orfs_{tool}.fasta.temp"))
  log:
    os.path.join(LOG_DIR, "{sample}", "{sample}_find_orfs_{tool}.log")
  benchmark:
    os.path.join(OUT_DIR, "benchmarks", "find_orfs", "{sample}_find_orfs_{tool}.tsv")
  conda:
    PALM
  shell:
    """
    # Ativa modo de debug para ver o erro no log
    set -x 
    
    # Limpa arquivo de saida
    > {output.orfs}

    # Evita falha do getorf caso o arquivo de entrada nao possua sequencias
    if [ ! -s {input.fasta} ]; then
        echo "[INFO] Nenhum nohit identificado. Gerando arquivo vazio." >> {log}
        exit 0
    fi

    # Lista explicita para evitar erro de expansao do shell
    for table in 1 3 4 5 6 11 16; do
        tmp_file="$(dirname {output.orfs})/tmp_Orfs_{wildcards.sample}_{wildcards.tool}_Table${{table}}.fasta"     
        getorf -sequence {input.fasta} -minsize 600 -table "$table" -find 0 -outseq "$tmp_file" 2>> {log} || true
        if [ -s "$tmp_file" ]; then
            sed -i.bak "s/^>/>gc_${{table}}_/g" "$tmp_file" 2>> {log}
            cat "$tmp_file" >> {output.orfs}
        fi
        
        rm -f "$tmp_file" "$tmp_file.bak"
    done
    """

#######################################################
# --- 3. Remove Redundancy from final Orfs Fasta --- #
###################################################### 
rule cd_hit:
  """
  Remove duplicated ORFs using CD-HIT.
  """
  input:
    fasta = os.path.join(OUT_DIR, "{sample}", "darkmatter_{tool}", "{sample}_find_orfs_{tool}.fasta.temp")
  output:  
    orfs = os.path.join(OUT_DIR, "{sample}", "darkmatter_{tool}", "{sample}_cd_hit_{tool}.fasta")
  log:
    os.path.join(LOG_DIR, "{sample}", "{sample}_cd_hit_{tool}.log")
  benchmark:
    os.path.join(OUT_DIR, "benchmarks", "cd_hit", "{sample}_cd_hit_{tool}.tsv")
  conda:
    PALM
  shell:
    """
    if [ -s {input.fasta} ]; then
        #remove duplicated orfs
        cd-hit -i {input.fasta} -o {output.orfs} -c 0.9 -d 1 -T {resources.threads} -M {resources.mem_mb} 2>> {log} 2>&1

        #substitute spaces per _
        awk 'substr($0, 1, 1) == ">" {{gsub(" ", "_")}} 1' {output.orfs} > {output.orfs}.tmp && mv {output.orfs}.tmp {output.orfs}
    else
        echo "[INFO] Nenhuma sequência para agrupar. Criando arquivo vazio." > {log}
        touch {output.orfs}
    fi
    """

###################################################
# --- 4. Iddentify putative RdRp within Data --- #
################################################## 
rule palm_annot:
  input:
    fasta = os.path.join(OUT_DIR, "{sample}", "darkmatter_{tool}", "{sample}_cd_hit_{tool}.fasta")
  output:
    fev  = os.path.join(OUT_DIR, "{sample}", "darkmatter_{tool}", "{sample}_palm_annot_{tool}.fev"),
    rdrp = os.path.join(OUT_DIR, "{sample}", "darkmatter_{tool}", "{sample}_palm_annot_{tool}.fasta")
  params:
    seqtype = config["palm_annot"]["seqtype"],
    minscore = config["palm_annot"]["minscore"],
    minpssmscore = config["palm_annot"]["minpssmscore"]
  conda:
    PALM
  log:
    os.path.join(LOG_DIR, "{sample}", "{sample}_palm_annot_{tool}.log")  
  benchmark:
    os.path.join(OUT_DIR, "benchmarks", "palm_annot", "{sample}_palm_annot_{tool}.tsv")
  script:
    "../scripts/palm_annot_run.py"
  
##################################
# --- 5. Convert FEV to TSV --- #
################################# 

rule fev2tsv_single:
  input:
    fev = os.path.join(OUT_DIR, "{sample}", "darkmatter_{tool}", "{sample}_palm_annot_{tool}.fev")
  output:
    tsv = os.path.join(OUT_DIR, "{sample}", "darkmatter_{tool}", "{sample}_fev2tsv_single_{tool}.tsv")
  conda:
    PALM
  log:
    os.path.join(LOG_DIR, "{sample}", "{sample}_fev2tsv_single_{tool}.log")
  benchmark:
    os.path.join(OUT_DIR, "benchmarks", "fev2tsv_single", "{sample}_fev2tsv_single_{tool}.tsv")
  script:
    "../scripts/fev2tsv_run.py"

##############################################
# --- 6. Summarize Dusk Matter Findings --- #
############################################## 
rule report_summarize:
    input:
        fasta = get_contigs_path,
        lineages = os.path.join(OUT_DIR, "{sample}", "diamond_{tool}", "{sample}_append_lineage_{tool}.tsv"),
        #basta_lca = os.path.join(OUT_DIR, "{sample}", "basta_{tool}", "{sample}_{tool}_lca.tsv"),
        dusk = os.path.join(OUT_DIR, "{sample}", "darkmatter_{tool}", "{sample}_fev2tsv_single_{tool}.tsv")
    output:
        html = os.path.join(OUT_DIR, "{sample}", "darkmatter_{tool}", "{sample}_report_summarize_{tool}.html")
    params:
        out_dir = lambda w, output: os.path.dirname(output.html),
        logos = os.path.join(RESOURCES_DIR, "logo")
    log:
        os.path.join(LOG_DIR, "{sample}", "{sample}_report_summarize_{tool}.log")
    benchmark:
        os.path.join(OUT_DIR, "benchmarks", "report_summarize", "{sample}_report_summarize_{tool}.tsv")
    conda:
        REPORT
    shell:
        """
        # Passamos o arquivo do BASTA no argumento --diamond_path.
        # O script R foi atualizado para detectar automaticamente que é um arquivo BASTA.
        Rscript {workflow.basedir}/scripts/generate_report.R \
            --sample_name {wildcards.sample} \
            --fasta_path {input.fasta} \
            --diamond_path {input.lineages} \
            --darkmatter_path {input.dusk} \
            --output_dir {params.out_dir} \
            --report_name $(basename {output.html}) \
            --logos {params.logos} \
            --input {workflow.basedir}/scripts/Report_Model.Rmd \
            > {log} 2>&1
        """

#########################################################
# --- 7. Extract ORFs Putative Raw RdRp and Contigs --- #
#########################################################
rule dm_validate:
  input:
    fasta = get_contigs_path,
    orfs  = os.path.join(OUT_DIR, "{sample}", "darkmatter_{tool}", "{sample}_cd_hit_{tool}.fasta"),
    dusk  = os.path.join(OUT_DIR, "{sample}", "darkmatter_{tool}", "{sample}_fev2tsv_single_{tool}.tsv")
  output:
    raw_rdrp    = os.path.join(OUT_DIR, "{sample}", "darkmatter_{tool}", "{sample}_dm_validate_{tool}_orfs.fasta"),
    raw_contigs = os.path.join(OUT_DIR, "{sample}", "darkmatter_{tool}", "{sample}_dm_validate_{tool}_contigs.fasta")
  params:
    minscore = config["palm_annot"]["minscore"]
  log:
    os.path.join(LOG_DIR, "{sample}", "{sample}_dm_validate_{tool}.log")
  benchmark:
    os.path.join(OUT_DIR, "benchmarks", "dm_validate", "{sample}_dm_validate_{tool}.tsv")
  conda:
    REPORT
  script:
    "../scripts/dm_validate.py"