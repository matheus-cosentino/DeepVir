@echo off
TITLE DeepVir - Desktop GUI (Windows + WSL)
cd /d "%~dp0"

REM 1. Try native Windows Python first (zero WSLg/X11 window overhead)
where python >nul 2>&1
if %ERRORLEVEL% EQU 0 (
    python -c "import tkinter" >nul 2>&1
    if %ERRORLEVEL% EQU 0 (
        python DeepVir.py --gui
        goto :EOF
    )
)

REM 2. Try Windows Python Launcher (py -3)
where py >nul 2>&1
if %ERRORLEVEL% EQU 0 (
    py -3 -c "import tkinter" >nul 2>&1
    if %ERRORLEVEL% EQU 0 (
        py -3 DeepVir.py --gui
        goto :EOF
    )
)

REM 3. Fallback: Run via WSL (activating Conda 'DeepVir' and ensuring tk-xft is installed for scalable fonts)
wsl bash -c "for c in \"$HOME/mambaforge/etc/profile.d/conda.sh\" \"$HOME/miniforge3/etc/profile.d/conda.sh\" \"$HOME/miniconda3/etc/profile.d/conda.sh\" \"$HOME/anaconda3/etc/profile.d/conda.sh\" \"/usr/local/bioinfo/Miniforge3-24.9.2-0/etc/profile.d/conda.sh\"; do [ -f \"$c\" ] && . \"$c\" && break; done; conda activate DeepVir 2>/dev/null || true; if ls \"$CONDA_PREFIX\"/conda-meta/tk-*noxft*.json >/dev/null 2>&1; then echo '[INFO] Updating Conda Tk from headless noxft to xft for scalable GUI fonts (one-time setup)...'; conda install -y -c conda-forge 'tk=*=xft_*' fonts-conda-ecosystem; fi; python3 DeepVir.py --gui"
pause
