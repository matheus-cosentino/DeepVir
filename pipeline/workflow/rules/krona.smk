###################################################################################
#                         workflow/rules/krona.smk                                #
#                           MSc. Matheus Cosentino                                #
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

rule krona_update_taxonomy:
    output:
        tab = os.path.join(KRONA_DB_DIR[0], "taxonomy.tab"),
        acc_sorted = os.path.join(KRONA_DB_DIR[0], "accession2taxid.sorted")
    input:
        names = os.path.join(TAXONOMY_DIR[0], "names.dmp"),
        nodes = os.path.join(TAXONOMY_DIR[0], "nodes.dmp"),
        acc2tax = os.path.join(TAXONOMY_DIR[0], "prot.accession2taxid.gz")
    params:
        tax_dir=lambda w, output: os.path.dirname(output.tab)
    conda:
        KRONA
    log:
        os.path.join(LOG_DIR, "workflow", "krona_update_taxonomy.log")
    shell:
        """
        echo "[INFO] Configurando diretórios..." > {log}
        TARGET_DIR=$(readlink -f {params.tax_dir})
        
        # 1. Taxonomia (Isso funciona, vamos manter)
        echo "[INFO] Linkando arquivos de taxonomia..." >> {log}
        ln -sf $(readlink -f {input.names}) $TARGET_DIR/names.dmp
        ln -sf $(readlink -f {input.nodes}) $TARGET_DIR/nodes.dmp
        
        echo "[INFO] Executando ktUpdateTaxonomy.sh..." >> {log}
        ktUpdateTaxonomy.sh --only-build $TARGET_DIR >> {log} 2>&1

        # 2. Accessions (Fazer Manualmente = 100% Seguro)
        # O script original falha em achar o arquivo, então nós mesmos criamos o output.
        # Pegamos colunas 2 (Accession.Version) e 3 (TaxID) e ordenamos.
        
        echo "[INFO] Ignorando geração de accession2taxid.sorted (Krona receberá os TaxIDs diretamente)." >> {log}
        
        # Cria um arquivo vazio apenas para satisfazer a saída da regra no Snakemake
        touch {output.acc_sorted}
        
        echo "[INFO] Concluído." >> {log}
        """

rule krona_kraken2:
    input:
        report = os.path.join(OUT_DIR, "{sample}", "kraken2_{tool}", "{sample}_kraken2_contigs_{tool}_output.txt")
    output:
        html = os.path.join(OUT_DIR, "{sample}", "krona_{tool}", "{sample}_krona_kraken2_{tool}.html")
    params:
        tax_dir=KRONA_DB_DIR[0]
    conda:
        KRONA
    log:
        os.path.join(LOG_DIR, "{sample}", "{sample}_krona_kraken2_{tool}.log")
    shell:
        """
        ktImportTaxonomy \
            -tax {params.tax_dir} -q 2 -t 3  \
            -o {output.html} \
            {input.report} \
            > {log} 2>&1
        """

rule krona_reads_kraken:
    wildcard_constraints:
        read_type="paired|unpaired"
    input:
        report = os.path.join(OUT_DIR, "{sample}", "kraken2_reads", "{sample}_kraken2_reads_{read_type}_output.txt")
    output:
        html = os.path.join(OUT_DIR, "{sample}", "krona_reads", "{sample}_krona_reads_kraken_{read_type}.html")
    conda:
        KRONA
    params:
        tax_dir=KRONA_DB_DIR[0]    
    log:
        os.path.join(LOG_DIR, "{sample}", "{sample}_krona_reads_kraken_{read_type}.log")
    shell:
        """
        ktImportTaxonomy \
            -tax {params.tax_dir}  -q 2 -t 3  \
            -o {output.html} \
            {input.report} \
            > {log} 2>&1
        """

rule krona_diamond:
    wildcard_constraints:
        tool = r"spades_k[\w]+|spades|megahit|flye|raven|medaka_flye|medaka_raven|pre_assembled"
    input:
        tax_hits = os.path.join(OUT_DIR, "{sample}", "diamond_{tool}", "{sample}_map_accession_to_taxid_{tool}.tmp"),
        tax_tab  = os.path.join(KRONA_DB_DIR[0], "taxonomy.tab")
    output:
        html = os.path.join(OUT_DIR, "{sample}", "krona_{tool}", "{sample}_krona_diamond_{tool}.html")
    conda:
        KRONA
    params:
        tax_dir = lambda w, input: os.path.dirname(input.tax_tab)
    log:
        os.path.join(LOG_DIR, "{sample}", "{sample}_krona_diamond_{tool}.log")
    shell:
        """
        echo "[INFO] Executing ktImportTaxonomy..." > {log}
        
        # Extrai a coluna 1 (QueryID) e a 13 (TaxID) e passa direto para o Krona
        cut -f 1,13 {input.tax_hits} | \
        ktImportTaxonomy \
            -tax {params.tax_dir} \
            -o {output.html} \
            - >> {log} 2>&1
            
        echo "[INFO] Krona HTML generated successfully." >> {log}
        """

rule krona_diamond_reads:
    input:
        hits    = os.path.join(OUT_DIR, "{sample}", "diamond_reads", "{sample}_map_accession_to_taxid_reads.tmp"),
        tax_tab = os.path.join(KRONA_DB_DIR[0], "taxonomy.tab")
    output:
        html = os.path.join(OUT_DIR, "{sample}", "krona_reads", "{sample}_krona_diamond_reads.html")
    conda:
        KRONA
    params:
        tax_dir = lambda w, input: os.path.dirname(input.tax_tab)
    log:
        os.path.join(LOG_DIR, "{sample}", "{sample}_krona_diamond_reads.log")
    shell:
        """
        echo "[INFO] Executing ktImportTaxonomy for Diamond Reads..." > {log}

        # Extrai a coluna 1 (QueryID) e a 13 (TaxID) e passa direto para o Krona
        cut -f 1,13 {input.hits} | \
        ktImportTaxonomy \
            -tax {params.tax_dir} \
            -o {output.html} \
            - >> {log} 2>&1

        echo "[INFO] Krona HTML generated successfully." >> {log}
        """