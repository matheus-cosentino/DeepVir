# scripts/fev2tsv_run.py

import subprocess
import os
import sys

# The 'snakemake' object is available
s = snakemake 
# Discover palm_annot directory dynamically relative to this script
palm_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "palm_annot")
log_file = s.log[0]

# --- 0. Check for Empty Input ---
# Se o arquivo FEV estiver vazio, cria o TSV em branco e sai com sucesso
if os.path.getsize(s.input.fev) == 0:
    with open(log_file, "w") as log_handle:
        log_handle.write("[INFO] Input FEV vazio. Pulando conversão para TSV.\n")
    open(s.output.tsv, "w").close()
    sys.exit(0)

# --- 1. Determine Executables ---
# Use the python that is currently executing this script (the Conda python)
python_exe = sys.executable 

# The fev2tsv.py script is located within the palm_annot directory structure
fev2tsv_script = os.path.join(palm_dir, "py", "fev2tsv.py")

# --- 2. Construct the Command List ---
# Execute the external script using the explicit Python interpreter
command = [
    python_exe,
    fev2tsv_script, 
    "--input", s.input.fev,
    "--output", s.output.tsv
]

# --- 3. Set the Custom PATH Environment Variable ---
# The tool might call other scripts/executables, so we must set the PATH.
new_path = f"{palm_dir}/bin:{palm_dir}/py:{os.environ.get('PATH', '')}"
env_vars = os.environ.copy()
env_vars['PATH'] = new_path

# --- 4. Execute the Command and Redirect Errors ---
with open(log_file, "w") as log_handle:
    try:
        # subprocess.run executes the command
        subprocess.run(
            command,
            check=True,       # Raise CalledProcessError on non-zero exit status
            stderr=log_handle,  # Redirect stderr (2>) to the log file
            env=env_vars      # Pass the updated environment variables
        )
    except subprocess.CalledProcessError as e:
        print(f"fev2tsv failed with exit code: {e.returncode}", file=sys.stderr)
        sys.exit(1)
    except FileNotFoundError:
        # This catches errors if the Python executable or the script itself can't be found
        print(f"Error: Executable not found. Python path: {python_exe}", file=sys.stderr)
        sys.exit(1)