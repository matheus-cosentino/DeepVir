#!/bin/bash
# Encerra o script se qualquer comando falhar, variável não definida for usada ou se um comando em pipeline falhar
set -euo pipefail

# --- Configurações do job SLURM ---
#SBATCH --job-name=DeepVir                       # Nome do job na fila
#SBATCH --cpus-per-task=1                        # CPUs do wrapper (o pipeline gerencia seus próprios threads)
#SBATCH --mem=32G                                 # Memória do wrapper
#SBATCH --time=3-00:00:00                        # Tempo máximo: 3 dias
#SBATCH --output=DeepVir_%j.out                  # Log de saída (%j = ID do job)

echo "
#####M#A#C#C##2#0#2#6###############I#R#D##&##U#F#R#J#####
## Name ___________________ DeepVir.slurm              ##
## Version ________________ v.2026.07                   ##
## Author _________________ Matheus Cosentino           ##
## Affiliation ____________ IRD TransVIHMI & UFRJ LDDV  ##
#####M#A#C#C##2#0#2#6###############I#R#D##&##U#F#R#J#####
"

echo "##### DeepVir SLURM Wrapper #####"
echo "[INFO] Início: $(date)"

# --- Valida argumentos: exige ao menos o diretório de dados ---
if [[ $# -lt 1 ]]; then
    echo "[ERROR] Uso: sbatch $(basename "$0") <DATA_DIR> [OUT_DIR]"
    echo "[ERROR] Exemplo: sbatch slurm_scripts/crab_novavir.slurm /data/meu_projeto /home/user/resultados"
    exit 1
fi

# -------- KEEP ---------
# --- 1. DEFINIÇÃO DE CAMINHOS ---

# Diretório do projeto (de onde o job foi submetido); deve estar no sistema de arquivos compartilhado
PROJECT_HOME="${SLURM_SUBMIT_DIR}"

# Variáveis dinâmicas (mudam a cada execução)
DATA_DIR="$1"        # Dados de entrada (1º argumento)
FINAL_HOME_DIR="${2:-/data04/projects04/MarceloSoares/viromapn/results/DeepVir/results_${SLURM_JOB_ID}}"      # Destino final dos resultados (2º argumento, ou padrão)

# Bancos de dados fixos
DIAMOND_BANK=/data03/databases/sources/ncbi/blast/dbV5/nr    # Banco DIAMOND (proteínas NR do NCBI)
KRAKEN2_BANK=/data04/databases/sources/KRAKEN2_DB/PlusPF     # Banco Kraken2 (classificação taxonômica)

# Diretórios no scratch local do nó (I/O rápido durante o processamento)
SCRATCH_DIR="/scratch/USERS/$USER/DeepVir_$SLURM_JOB_ID"
TEMP_DIR="/scratch/USERS/$USER/temp_DeepVir_$SLURM_JOB_ID"

mkdir -p "$SCRATCH_DIR" "$TEMP_DIR" "$FINAL_HOME_DIR"

# --- Handler de limpeza (executado automaticamente ao sair) ---
# Sucesso: remove o scratch. Falha: preserva o scratch para inspeção.
cleanup() {
    local exit_code=$?
    if [[ $exit_code -ne 0 ]]; then
        echo "[WARN] Job encerrado com código $exit_code. Scratch preservado para inspeção: $SCRATCH_DIR"
    else
        echo "[INFO] Limpando scratch..."
        rm -rf "${SCRATCH_DIR}" "${TEMP_DIR}"
        echo "[INFO] Scratch limpo com sucesso."
    fi
    echo "[INFO] Fim: $(date)"
}
trap cleanup EXIT

# --- 2. ATIVAÇÃO DO AMBIENTE ---
module load miniconda3
eval "$(conda shell.bash hook)"   # Inicializa o shell para uso do conda

cd "$PROJECT_HOME" || { echo "Error: Could not find $PROJECT_HOME"; exit 1; }
echo "[INFO] Executando a partir de: $(pwd)"
echo "[INFO] Dados de entrada:     ${DATA_DIR}"
echo "[INFO] Processamento em:     ${SCRATCH_DIR}"
echo "[INFO] Resultados finais em: ${FINAL_HOME_DIR}"

# -------- CHANGE ---------
# --- 3. EXECUÇÃO DO PIPELINE ---
bash DeepVir.sh \
    --input "${DATA_DIR}" \
    --output "${SCRATCH_DIR}" \
    --assembly \
    --diamond \
    --reads-diamond \
    --reads-kraken \
    --kraken2 \
    --darkmatter \
    --diamond_db "${DIAMOND_BANK}" \
    --kraken2_db "${KRAKEN2_BANK}" \
    --temp-dir "${TEMP_DIR}"

# --- 4. RECUPERAÇÃO DOS RESULTADOS ---
# Copia do scratch para o destino final antes que o job encerre
echo "[INFO] Copiando resultados do scratch..."
rsync -av --progress --checksum "${SCRATCH_DIR}/" "${FINAL_HOME_DIR}/"

echo "[INFO] Concluído. Resultados disponíveis em: ${FINAL_HOME_DIR}"