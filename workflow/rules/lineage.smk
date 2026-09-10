###################################################################################
#                       workflow/rules/lineage.smk                                #
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

#################################################
# --- 1. Map Proteins Id Diamond for Taxid --- #
################################################ 
rule download_prot:
  """
  Download Protein acc to taxid data
  """
  output:
    gz = os.path.join(TAXONOMY_DIR[0], "prot.accession2taxid.gz"),
    nodes = os.path.join(TAXONOMY_DIR[0], "nodes.dmp"),
    names = os.path.join(TAXONOMY_DIR[0], "names.dmp")
  log:
    os.path.join(OUT_DIR, "log", "tax_download_mapping.log")
  benchmark:
    os.path.join(OUT_DIR, "benchmarks", "download_prot", "taxonomy.tsv")
  params:
    taxon_gz = os.path.join(TAXONOMY_DIR[0], "taxdump.tar.gz"),
    tax_dir = lambda w, output: os.path.dirname(output.gz)
  conda:
    DOWNLOAD
  shell:
    """
    wget -c https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/prot.accession2taxid.gz -O {output.gz} 2> {log}
    # Remove existing taxdump to avoid 'trailing garbage' errors with wget -c if the remote file changed size
    rm -f {params.taxon_gz}
    wget https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz -O {params.taxon_gz} 2>> {log}
    
    tar -zxvf {params.taxon_gz} -C {params.tax_dir} names.dmp nodes.dmp 2>> {log}
    """


#################################################
# --- 2. Map Proteins Id Diamond for Taxid --- #
################################################ 

rule map_accession_to_taxid:
    """
    Maps protein IDs (Subject ID, column 2 of DIAMOND) to TaxIDs.
    """
    wildcard_constraints:
      tool = r"spades_k[\w]+|spades|megahit|flye|raven|medaka_flye|medaka_raven|pre_assembled|reads"
    input:
      hit_file = os.path.join(OUT_DIR, "{sample}", "diamond_{tool}", "{sample}_{tool}_report.txt"),
      taxid_map = rules.download_prot.output.gz
    output:
      ids = os.path.join(OUT_DIR, "{sample}", "diamond_{tool}", "{sample}_{tool}_hits_with_taxid.tmp")
    shadow: 
      "minimal"
    log:
      os.path.join(OUT_DIR, "{sample}", "log", "{sample}_{tool}_map_acc_prot.log")
    benchmark:
      os.path.join(OUT_DIR, "benchmarks", "map_accession_to_taxid", "{sample}_{tool}.tsv")
    conda:
      CORE
    script:
      "../scripts/map_acc_to_taxid.py"


##################################################
# --- 3. Split Hits for Taxid and Not Found --- #
################################################# 
rule split_hits_by_taxid:
    """
    Splits the input file into two: one with valid taxids and one with "NOT_FOUND".
    """
    wildcard_constraints:
      tool = r"spades_k[\w]+|spades|megahit|flye|raven|medaka_flye|medaka_raven|pre_assembled|reads"
    input:
        os.path.join(OUT_DIR, "{sample}", "diamond_{tool}", "{sample}_{tool}_hits_with_taxid.tmp")
    output:
        valid_hits = os.path.join(OUT_DIR, "{sample}", "diamond_{tool}", "{sample}_{tool}_hits_with_header.tsv"),
        no_lineage = temp(os.path.join(OUT_DIR, "{sample}", "diamond_{tool}", "{sample}_{tool}_hits_no_lineage.temp"))
    params:
        header="qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\ttaxid"
    log:
        os.path.join(OUT_DIR, "{sample}", "log", "{sample}_{tool}_split_hits_by_taxid.log")
    benchmark:
        os.path.join(OUT_DIR, "benchmarks", "split_hits_by_taxid", "{sample}_{tool}.tsv")
    conda:
        CORE
    shell:
        """
        # Create output files with headers
        echo -e "{params.header}" > {output.valid_hits}
        echo -e "{params.header}" > {output.no_lineage}
        
        # CORREÇÃO AQUI: Adicionado BEGIN {{OFS="\\t"}}
        # Isso garante que ao modificar a linha, o awk mantenha os TABs e não troque por espaços.
        
        awk -F'\\t' 'BEGIN {{OFS="\\t"}} 
        {{
            # Limpa carriage return (\r) do Windows se existir
            gsub(/\\r/, "", $NF)
            
            if ($NF == "NOT_FOUND") {{
                print $0 >> "{output.no_lineage}"
            }} else {{
                print $0 >> "{output.valid_hits}"
            }}
        }}' {input}
        """


##################################################
# --- 3. Split Hits for Taxid and Not Found --- #
################################################# 
rule append_lineage:
    wildcard_constraints:
      tool = r"spades_k[\w]+|spades|megahit|flye|raven|medaka_flye|medaka_raven|pre_assembled|reads"
    input:
        valid_hits = rules.split_hits_by_taxid.output.valid_hits
    output:
        lineages = os.path.join(OUT_DIR, "{sample}", "diamond_{tool}", "{sample}_{tool}_hits_with_lineage.tsv")
    params:
        base_header="qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\tTaxid",
        lineage_header="Lineage\tCelular\tAcelular\tRealm\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies",
        nodes=os.path.join(TAXONOMY_DIR[0], "nodes.dmp")
    log:
        os.path.join(OUT_DIR, "{sample}", "log", "{sample}_{tool}_diamond_lineages.log")
    benchmark:
        os.path.join(OUT_DIR, "benchmarks", "append_lineage", "{sample}_{tool}.tsv")
    conda:
        TAXONKIT
    shell:
        """
        DB_DIR=$(dirname {params.nodes})
        
        # Skip header (criado na regra anterior) e extrai apenas a coluna TaxID (13)
        tail -n +2 {input.valid_hits} | cut -f 13 > {output}.taxids.tmp
        
        # 1. Get lineage information
        # TaxonKit retorna: TaxID <TAB> Lineage <TAB> Realm <TAB> Kingdom ...
        taxonkit lineage --threads {resources.threads} --data-dir "${{DB_DIR}}" {output}.taxids.tmp 2>> {log} | \
        taxonkit reformat --threads {resources.threads} --data-dir "${{DB_DIR}}" \
            -f "{{C}}\\t{{a}}\\t{{d}}\\t{{k}}\\t{{p}}\\t{{c}}\\t{{o}}\\t{{f}}\\t{{g}}\\t{{s}}" \\
            -F \
             2>> {log} > {output}.lineage.tmp
        
        # 2. Extract reformatted columns only (drop TaxID na col 1 do output do taxonkit)
        cut -f 2- {output}.lineage.tmp > {output}.reformat_only.tmp
        
        # 3. Create final output header
        echo -e "{params.base_header}\\t{params.lineage_header}" > {output.lineages}
        
        # 4. Paste original columns + reformatted lineage
        # Pegamos colunas 1-12 originais + coluna 13 (taxid) original + linhagem nova
        tail -n +2 {input.valid_hits} | cut -f 1-12 | \\
        paste - <(tail -n +2 {input.valid_hits} | cut -f 13) {output}.reformat_only.tmp >> {output.lineages}
        
        # Cleanup
        rm {output}.taxids.tmp {output}.lineage.tmp {output}.reformat_only.tmp
        """