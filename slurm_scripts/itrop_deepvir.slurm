#!/bin/bash
# Exit immediately if a command fails, an undefined variable is used, or a pipeline fails
set -euo pipefail

# --- SLURM Job Settings ---
#SBATCH --job-name=DeepVir                      # Job name in the queue
#SBATCH --cpus-per-task=1                        # CPUs for the wrapper (the pipeline manages its own threads)
#SBATCH --mem=64G                                # Wrapper memory (measured: ~54GB MaxRSS via sacct; includes mamba solver + Snakemake DAG + SLURM job monitoring)
#SBATCH --time=7-00:00:00                        # Maximum runtime: 7 days
#SBATCH --constraint=infiniband                  # Node constraint


echo "
#####M#A#C#C##2#0#2#6###############I#R#D##&##U#F#R#J#####
## Name ___________________ DeepVir.slurm              ##
## Version ________________ v.2026.07                   ##
## Author _________________ Matheus Cosentino           ##
## Affiliation ____________ IRD TransVIHMI & UFRJ LDDV  ##
#####M#A#C#C##2#0#2#6###############I#R#D##&##U#F#R#J#####
"

echo "##### DeepVir SLURM Wrapper #####"
echo "[INFO] Start: $(date)"

# --- Validate arguments: requires at least the data directory ---
if [[ $# -lt 1 ]]; then
    echo "[ERROR] Usage: sbatch $(basename "$0") <DATA_DIR> [OUT_DIR]"
    echo "[ERROR] Example: sbatch slurm_scripts/itrop_deepvir.slurm /data/my_project /home/user/results"
    exit 1
fi

# -------- KEEP ---------
# --- 1. PATH DEFINITIONS ---

# Project directory (where the job was submitted from); must be on a shared filesystem
PROJECT_HOME="${SLURM_SUBMIT_DIR}"

# Dynamic variables (change with each run)
DATA_DIR="$1"        # Input data (1st argument)
FINAL_HOME_DIR="${2:-/scratch-ib/$USER/DeepVir/results_${SLURM_JOB_ID}}"      # Final results destination (2nd argument, or default)

# Fixed databases
DIAMOND_BANK=/share/banks/diamond/nr/21-09/diamond_nr.dmnd     # DIAMOND database (NCBI NR proteins)
KRAKEN2_BANK=/share/banks/kraken2/nt/21-09/nt/                 # Kraken2 database (taxonomic classification)

# Scratch directories on the node (fast I/O during processing)
SCRATCH_DIR="/scratch-ib/$USER/DeepVir_$SLURM_JOB_ID"
TEMP_DIR="/scratch-ib/$USER/temp_DeepVir_$SLURM_JOB_ID"

mkdir -p "$SCRATCH_DIR" "$TEMP_DIR" "$FINAL_HOME_DIR"

# --- Robust directory removal (handles NFS stale locks and .condatmp residues) ---
# Retries rm -rf up to 5 times with increasing delays.
# On NFS, conda (.condatmp) and Snakemake (.snakemake/log) can leave open handles
# that cause "directory not empty" even when no process is actively writing.
robust_rm() {
    local target="$1"
    local max_attempts=5
    local attempt=1
    while [[ -e "$target" && $attempt -le $max_attempts ]]; do
        rm -rf "$target" 2>/dev/null && return 0
        echo "[WARN] rm attempt $attempt/$max_attempts failed for: $target — retrying in ${attempt}s..."
        sleep "$attempt"
        (( attempt++ ))
    done
    if [[ -e "$target" ]]; then
        echo "[WARN] Could not fully remove: $target"
        echo "[WARN] Remaining contents:"
        find "$target" -maxdepth 3 2>/dev/null | head -20
        return 1
    fi
    return 0
}

# --- Cleanup handler (runs automatically on exit) ---
# Success: removes scratch. Failure: preserves scratch for inspection.
cleanup() {
    local exit_code=$?
    if [[ $exit_code -ne 0 ]]; then
        echo "[WARN] Job exited with code $exit_code. Scratch preserved for inspection: $SCRATCH_DIR"
    else
        echo "[INFO] Cleaning up scratch..."
        robust_rm "${SCRATCH_DIR}"
        robust_rm "${TEMP_DIR}"
        echo "[INFO] Scratch cleaned successfully."
    fi
    echo "[INFO] End: $(date)"
}
trap cleanup EXIT

# --- 2. ENVIRONMENT ACTIVATION ---
module load bioinfo-itrop
module load miniconda3
eval "$(conda shell.bash hook)"   # Initialize shell for conda usage

cd "$PROJECT_HOME" || { echo "Error: Could not find $PROJECT_HOME"; exit 1; }
echo "[INFO] Running from: $(pwd)"
echo "[INFO] Input data:       ${DATA_DIR}"
echo "[INFO] Processing in:    ${SCRATCH_DIR}"
echo "[INFO] Final results at: ${FINAL_HOME_DIR}"

# -------- CHANGE ---------
# --- 3. PIPELINE EXECUTION ---
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
    --temp-dir "${TEMP_DIR}" \
    --profile profiles/profile_slurm/

# --- 4. RESULTS RETRIEVAL ---
# Copy from scratch to final destination before the job ends
echo "[INFO] Copying results from scratch..."
rsync -av --progress --checksum "${SCRATCH_DIR}/" "${FINAL_HOME_DIR}/"

module purge
echo "[INFO] Done. Results available at: ${FINAL_HOME_DIR}"