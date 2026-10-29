import pandas as pd
import re
from Bio import SeqIO

# Variáveis globais providenciadas pelo Snakemake
log_file = snakemake.log[0]
input_tsv = snakemake.input.tsv
input_orfs = snakemake.input.orfs
input_fasta = snakemake.input.fasta
output_raw_rvdb_orfs = snakemake.output.raw_rvdb_orfs
output_raw_rvdb_contigs = snakemake.output.raw_rvdb_contigs

# Carrega o CSV sumarizado do RVDB
try:
    df = pd.read_csv(input_tsv, sep='\t')
    
    # --- FILTRO DE CONFIANÇA ALTA APENAS ---
    if 'Confidence' in df.columns:
        df = df[df['Confidence'] == 'High']
        
    target_orfs = set(df['Sequence_ID'].astype(str).tolist())
except (pd.errors.EmptyDataError, KeyError):
    target_orfs = set()
    
target_contigs = set()

# Logica para extrair o ID do Contig a partir do Label do ORF
for label in target_orfs:
    match = re.search(r'^gc_\d+_(.+?)_\d+_\[', label)
    if match:
        contig_id = match.group(1)
        target_contigs.add(contig_id)
    else:
        with open(log_file, "a") as f:
            f.write(f"[WARNING] Nao foi possivel parsear o contig ID do label: {label}\n")

# Salva os ORFs identificados
with open(output_raw_rvdb_orfs, "w") as out_orf:
    if len(target_orfs) > 0:
        for record in SeqIO.parse(input_orfs, "fasta"):
            if record.id in target_orfs:
                SeqIO.write(record, out_orf, "fasta")

# Salva os Contigs originais correspondentes 
with open(output_raw_rvdb_contigs, "w") as out_contig:
    if len(target_contigs) > 0:
        for record in SeqIO.parse(input_fasta, "fasta"):
            if record.id in target_contigs:
                SeqIO.write(record, out_contig, "fasta")
