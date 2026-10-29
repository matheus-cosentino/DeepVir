import pandas as pd
import re
from Bio import SeqIO

# Snakemake makes snakemake object available magically in script execution
log_file = snakemake.log[0]
dusk_file = snakemake.input.dusk
orfs_file = snakemake.input.orfs
fasta_file = snakemake.input.fasta
raw_rdrp = snakemake.output.raw_rdrp
raw_contigs = snakemake.output.raw_contigs
sample = snakemake.wildcards.sample
tool = snakemake.wildcards.tool

with open(log_file, "w") as logf:
    logf.write(f"[INFO] dm_validate started for {sample} / {tool}\n")

    # Lê o TSV — pode estar vazio se nenhum RdRp foi encontrado
    try:
        df = pd.read_csv(dusk_file, sep='\t')
    except Exception as e:
        logf.write(f"[ERROR] Could not read RdRp TSV: {e}\n")
        df = pd.DataFrame()

    if df.empty or 'Label' not in df.columns:
        logf.write("[WARNING] No RdRp hits found. Creating empty output files.\n")
        open(raw_rdrp, 'w').close()
        open(raw_contigs, 'w').close()
    else:
        target_orfs = set(df['Label'].astype(str).tolist())
        target_contigs = set()

        # Lógica para extrair o ID do Contig a partir do Label do ORF
        for label in target_orfs:
            match = re.search(r'^gc_\d+_(.+?)_\d+_\[', label)
            if match:
                contig_id = match.group(1)
                target_contigs.add(contig_id)
            else:
                logf.write(f"[WARNING] Não foi possível parsear o contig ID do label: {label}\n")

        # Salva os ORFs identificados
        with open(raw_rdrp, "w") as out_orf:
            for record in SeqIO.parse(orfs_file, "fasta"):
                if record.id in target_orfs:
                    SeqIO.write(record, out_orf, "fasta")

        # Salva os Contigs originais correspondentes
        with open(raw_contigs, "w") as out_contig:
            for record in SeqIO.parse(fasta_file, "fasta"):
                if record.id in target_contigs:
                    SeqIO.write(record, out_contig, "fasta")

        logf.write(f"[INFO] Done. {len(target_orfs)} RdRp ORFs, {len(target_contigs)} contigs extracted.\n")
