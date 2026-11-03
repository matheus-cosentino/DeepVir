#!/bin/bash
#########################################################################
# DeepVir_gui.sh - DeepVir 2.0 1-Click Launcher for Linux & WSL         #
# Auto-detects Conda/Miniforge or System Python 3 on any machine        #
#########################################################################

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# 1. Source Conda if available on Linux / WSL / macOS
for _conda_sh in \
    "$HOME/miniforge3/etc/profile.d/conda.sh" \
    "$HOME/mambaforge/etc/profile.d/conda.sh" \
    "$HOME/miniconda3/etc/profile.d/conda.sh" \
    "$HOME/anaconda3/etc/profile.d/conda.sh" \
    "/usr/local/bioinfo/Miniforge3-24.9.2-0/etc/profile.d/conda.sh"; do
    if [ -f "$_conda_sh" ]; then
        . "$_conda_sh"
        conda activate DeepVir 2>/dev/null || true
        break
    fi
done

# 2. If active Conda env has headless tk-noxft (which locks X11 fonts to 6x13 fixed), upgrade to tk-xft once
if [ -n "$CONDA_PREFIX" ] && ls "$CONDA_PREFIX"/conda-meta/tk-*noxft*.json >/dev/null 2>&1; then
    echo "[INFO] Updating Conda Tk from headless 'noxft' to 'xft' for scalable GUI fonts (one-time setup)..."
    conda install -y -c conda-forge "tk=*=xft_*" fonts-conda-ecosystem
fi

# 3. Launch Universal Python Executor in GUI mode
if command -v python3 &>/dev/null; then
    exec python3 "$SCRIPT_DIR/DeepVir.py" --gui "$@"
elif command -v python &>/dev/null; then
    exec python "$SCRIPT_DIR/DeepVir.py" --gui "$@"
else
    echo "[ERROR] Python 3 was not found in PATH!"
    exit 1
fi
