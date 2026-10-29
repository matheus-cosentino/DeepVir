import pandas as pd
import sqlite3
import re
import urllib.request
import json
import time
from datetime import datetime
import sys

# --- CONFIGURAÇÕES (Snakemake) ---
INPUT_TBL = snakemake.input.tbl
DATABASE = snakemake.input.sqlite
OUTPUT_FILE = snakemake.output.csv
LOG_FILE = snakemake.log[0]

sys.stdout = open(LOG_FILE, "w")
sys.stderr = sys.stdout

# Palavras-chave que não agregam valor biológico
BLACKLIST = {
    'protein', 'hypothetical', 'unnamed', 'predicted', 'unknown', 
    'virus', 'sequence', 'structure', 'homolog', 'orf', 'rna', 'type', 'like'
}

# Cache para evitar consultas repetidas ao NCBI
tax_cache = {}

def get_scientific_name(taxid):
    """Consulta o nome científico no NCBI via API"""
    if not taxid or taxid == "N/A": return "Unknown"
    if taxid in tax_cache: return tax_cache[taxid]
    
    try:
        url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=taxonomy&id={taxid}&retmode=json"
        with urllib.request.urlopen(url) as response:
            data = json.loads(response.read().decode())
            name = data['result'][str(taxid)]['scientificname']
            tax_cache[taxid] = name
            time.sleep(0.1) # Respeitar limites do NCBI
            return name
    except:
        return f"TaxID:{taxid}"

def process_results():
    print("1. A ler ficheiro HMMER...")
    # O HMMER usa espaços variáveis, por isso usamos sep='\s+'
    cols = ["Target_ID", "query_name", "E_value", "Score"]
    try:
        df = pd.read_csv(INPUT_TBL, comment='#', header=None, sep=r'\s+', 
                         usecols=[0, 2, 4, 5], names=cols)
    except Exception as e:
        print(f"Erro ao ler .tbl: {e}")
        return

    print(f"2. Cruzando com SQLite ({len(df)} hits encontrados)...")
    conn = sqlite3.connect(DATABASE)
    final_data = []

    for _, row in df.iterrows():
        # Extrair número do ID (FAM012345 -> 12345)
        match = re.search(r'FAM(\d+)', str(row['Target_ID']))
        if not match: continue
        fam_id = int(match.group(1))

        # Query 1: Melhor palavra-chave (ignorando a blacklist)
        kw_query = f"""
            SELECT k.str FROM keyword k
            JOIN fam_kw_seqnames fks ON k.id = fks.kwId
            WHERE fks.famId = {fam_id}
            ORDER BY fks.freq DESC
        """
        kws = pd.read_sql_query(kw_query, conn)
        desc = "Unknown Protein"
        for word in kws['str']:
            if word.lower() not in BLACKLIST:
                desc = word
                break

        # Query 2: TaxID (LCA - Last Common Ancestor)
        tax_query = f"SELECT LCAtaxid FROM family WHERE id = {fam_id}"
        tax_res = pd.read_sql_query(tax_query, conn)
        taxid = tax_res['LCAtaxid'].iloc[0] if not tax_res.empty else "N/A"

        # Calcular nível de confiança
        e_val = float(row['E_value'])
        conf = "High" if e_val < 1e-10 else ("Medium" if e_val < 1e-5 else "Low")

        final_data.append({
            'Sequence_ID': row['query_name'],
            'HMM_Family': row['Target_ID'],
            'E_value': e_val,
            'Bit_Score': row['Score'],
            'Annotation': desc,
            'TaxID': taxid,
            'Confidence': conf
        })

    conn.close()
    
    # Criar DataFrame final
    df_final = pd.DataFrame(final_data)

    print("3. Organizando hits (Todos os hits mantidos)...")
    # Mantém todos os hits ordenados por qualidade (descrescentemente por E-value)
    df_best = df_final.sort_values(['Sequence_ID', 'E_value'])

    print("4. Traduzindo TaxIDs para nomes científicos (pode demorar)...")
    df_best['Scientific_Name'] = df_best['TaxID'].apply(get_scientific_name)

    # Reorganizar colunas para o relatório
    ordem = ['Sequence_ID', 'Annotation', 'Scientific_Name', 'Confidence', 'E_value', 'Bit_Score', 'HMM_Family', 'TaxID']
    df_best = df_best[ordem]

    # Salvar
    df_best.to_csv(OUTPUT_FILE, sep='\t', index=False)
    print(f"--- SUCESSO ---\nRelatório final guardado em: {OUTPUT_FILE}")

if __name__ == "__main__":
    process_results()
    sys.stdout.close()
