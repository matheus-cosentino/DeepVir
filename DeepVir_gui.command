#!/bin/bash
#########################################################################
# DeepVir_gui.command - 1-Click Launcher for macOS (Finder)             #
# Recognizes Miniforge, Mambaforge, Miniconda, Anaconda, Homebrew & OS  #
#########################################################################

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

# 1. Source Conda if available on macOS (Intel or Apple Silicon)
for _conda_sh in \
    "$HOME/miniforge3/etc/profile.d/conda.sh" \
    "$HOME/mambaforge/etc/profile.d/conda.sh" \
    "$HOME/miniconda3/etc/profile.d/conda.sh" \
    "$HOME/anaconda3/etc/profile.d/conda.sh" \
    "/opt/homebrew/Caskroom/miniforge/base/etc/profile.d/conda.sh" \
    "/opt/homebrew/anaconda3/etc/profile.d/conda.sh" \
    "/usr/local/Caskroom/miniforge/base/etc/profile.d/conda.sh"; do
    if [ -f "$_conda_sh" ]; then
        . "$_conda_sh"
        conda activate DeepVir 2>/dev/null || true
        break
    fi
done

# 2. Ensure Homebrew paths are in PATH on macOS
export PATH="/opt/homebrew/bin:/usr/local/bin:$PATH"

# 3. Locate recognized Python 3 with Tkinter and launch DeepVir GUI
if command -v python3 &>/dev/null; then
    exec python3 "$SCRIPT_DIR/DeepVir.py" --gui "$@"
elif command -v python &>/dev/null; then
    exec python "$SCRIPT_DIR/DeepVir.py" --gui "$@"
else
    echo "[ERROR] Python 3 was not found on this Mac. Please install Miniforge or Python 3."
    read -r -p "Press Enter to exit..."
    exit 1
fi

