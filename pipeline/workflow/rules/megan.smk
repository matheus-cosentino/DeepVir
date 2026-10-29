###################################################################################
#                         workflow/rules/megan.smk                                #
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

#1. GET MAPPER
if config["modules"].get("reads_diamond", False):

  rule get_megan6_mapper:
    output:
      db = os.path.join(MEGAN_DIR[0], "megan-map-Feb2022-ue.db"),
    params:
      dir = lambda w, output: os.path.dirname(output.db)
    log:
      os.path.join(LOG_DIR, "workflow", "megan6_download_mapper.log")
    conda:
      DOWNLOAD
    shell:
      """
      wget https://software-ab.cs.uni-tuebingen.de/download/megan6/megan-map-Feb2022-ue.db.zip -O {output.db}.zip 2> {log}
      unzip -o {output.db}.zip -d {params.dir} >> {log} 2>&1
      rm {output.db}.zip
      """


  rule reads_meganizer:
    input:
      daa = rules.diamond_blastx_reads.output.daa,
      db  = rules.get_megan6_mapper.output.db
    output:
      daa = os.path.join(OUT_DIR, "{sample}", "megan_reads", "{sample}_reads_report_meganized.daa")
    log:
      os.path.join(LOG_DIR, "{sample}", "meganizer_{sample}.log")
    conda:
      MEGAN
    shell:
      """
      cp {input.daa} {output.daa}

      # Adicionando parâmetros explícitos do LCA para evitar perda de hits
      # -ms 40: Diminui o min-score aceitável (útil para reads curtos/divergentes)
      # -me 0.01: Max e-value aceitável
      # -top 10: Limiar de 10% do melhor score para o LCA
      # -sup 1: Suporte mínimo de 1 read para confirmar um táxon (CRÍTICO)

      daa-meganizer \
          -i {output.daa} \
          -mdb {input.db} \
          -ms 40 \
          -me 0.01 \
          -top 10 \
          -sup 1 \
          -t {threads} > {log} 2>&1
      """


  rule reads_daa2info:
    input:
      daa = rules.reads_meganizer.output.daa
    output:
      summary = os.path.join(OUT_DIR, "{sample}", "megan_reads", "{sample}_reads_summary.megan")
    log:
      os.path.join(LOG_DIR, "{sample}", "daa2info_{sample}.log")
    conda:
      MEGAN
    shell:
      """
      daa2info -i {input.daa} --extractSummaryFile {output.summary} > {log} 2>&1
      """