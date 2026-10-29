# scripts/palm_annot_run.py

import subprocess
import os
import sys

s = snakemake
# Discover palm_annot directory dynamically relative to this script
palm_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "palm_annot")
log_file = s.log[0]

# --- 0. Check for Empty Input ---
# Se o arquivo FASTA estiver vazio, cria os arquivos de saída em branco e sai com sucesso
if os.path.getsize(s.input.fasta) == 0:
    with open(log_file, "w") as log_handle:
        log_handle.write("[INFO] Input FASTA vazio. Pulando palm_annot.\n")
    open(s.output.fev, "w").close()
    open(s.output.rdrp, "w").close()
    sys.exit(0)

# --- 1. Determine the Correct Python Executable ---
# Use the python that is currently executing this script (which should be the Conda one)
python_exe = sys.executable 

# The script to run is the external PALM script, which we run *using* the Conda Python
palm_script = os.path.join(palm_dir, "py", "palm_annot.py") 

if not os.path.isfile(palm_script):
    sys.exit(f"[ERROR] PALM script not found at: {palm_script}\n"
             f"Please ensure the 'workflow/scripts/palm_annot' directory is fully copied to your working repository.")


# --- 2. Construct the Command List ---
# The command starts with the explicit Python executable
command = [
    python_exe,          # Use the explicit Conda Python executable path
    palm_script,         # The external PALM script
    "--input", s.input.fasta,
    "--seqtype", s.params.seqtype,
    "--fev", s.output.fev,
    "--rdrp", s.output.rdrp,
    "--minscore", str(s.params.minscore),
    "--threads", str(s.resources.threads),
    "--minpssmscore", str(s.params.minpssmscore)
]

new_path = f"{palm_dir}/bin:{palm_dir}/py:{os.environ.get('PATH', '')}"
env_vars = os.environ.copy()
env_vars['PATH'] = new_path


# --- 4. Execute the Command ---
with open(log_file, "w") as log_handle:
    try:
        subprocess.run(
            command,
            check=True,
            stderr=log_handle
        )
    except subprocess.CalledProcessError as e:
        print(f"PALM annotation failed with exit code: {e.returncode}", file=sys.stderr)
        sys.exit(1)
    except FileNotFoundError:
        print(f"Error: The Python executable ({python_exe}) or PALM script was not found.", file=sys.stderr)
        sys.exit(1)
