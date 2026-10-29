import pandas as pd
import re
from Bio import SeqIO

# Variáveis globais do Snakemake
log_file = snakemake.log[0]
input_tsv = snakemake.input.tsv
input_orfs = snakemake.input.orfs
input_fasta = snakemake.input.fasta
output_raw_rvdb_orfs = snakemake.output.raw_rvdb_orfs
output_raw_rvdb_contigs = snakemake.output.raw_rvdb_contigs

# 1. Carrega o CSV tratando as aspas que aparecem nos seus dados 
try:
    # O quotechar='"' garante que o pandas remova as aspas dos IDs e anotações 
    df = pd.read_csv(input_tsv, sep='\t', quotechar='"')
    
    # 2. Filtro atualizado para incluir "reverse" (comum no RVDB para pols) 
    pol_pattern = r'pol|polymerase|transcriptase|polyprotein|rdrp|integrase|protease|retrotransposon|reverse|\brt\b'
    
    mask = pd.Series(False, index=df.index)
    for col in ['Annotation', 'Description', 'str']:
        if col in df.columns:
            mask |= df[col].str.contains(pol_pattern, case=False, na=False)
    
    df = df[mask]
    target_orfs = set(df['Sequence_ID'].astype(str).tolist())
except (pd.errors.EmptyDataError, KeyError):
    target_orfs = set()
    
target_contigs = set()

# 3. Extração do ID do Contig via Regex 
for label in target_orfs:
    # Captura o que está entre 'gc_X_' e o penúltimo '_X_[' 
    match = re.search(r'gc_\d+_(.+?)_\d+_\[', label)
    if match:
        contig_id = match.group(1)
        target_contigs.add(contig_id)
    else:
        with open(log_file, "a") as f:
            f.write(f"[WARNING] Falha no parse do ID: {label}\n") 

# 4. Salva os ORFs 
with open(output_raw_rvdb_orfs, "w") as out_orf:
    if len(target_orfs) > 0:
        for record in SeqIO.parse(input_orfs, "fasta"):
            if record.id in target_orfs: 
                SeqIO.write(record, out_orf, "fasta") 

# 5. Salva os Contigs (usando record.id completo) 
with open(output_raw_rvdb_contigs, "w") as out_contig:
    if len(target_contigs) > 0:
        for record in SeqIO.parse(input_fasta, "fasta"):
            # Se o seu FASTA tem o nome completo (NODE_10426_length...), 
            # o record.id deve bater com o contig_id extraído.
            if record.id in target_contigs: 
                SeqIO.write(record, out_contig, "fasta") 
