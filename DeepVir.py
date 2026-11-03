#!/usr/bin/env python3
"""
DeepVir.py - Universal Multi-Platform Python Executor (Windows+WSL, macOS, Linux)
=================================================================================
Author: MSc. Matheus Cosentino / DeepVir Team
Version: 2.0 | 10.2026

Project Directory Architecture:
  1. pipeline/ -> Core Snakemake workflow, rules, envs, configs, profiles, GUI & resources
  2. docs/     -> Licenses (GPLv3) and Markdown documentation guides
  3. logs/     -> Centralized execution logs, Snakemake rule logs & run_overrides history
=================================================================================
"""

import argparse
import getpass
import os
import platform
import re
import shutil
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

# =============================================================================
# CORE 3-DIRECTORY ARCHITECTURE PATHS
# =============================================================================
ROOT_DIR = Path(__file__).resolve().parent
PIPELINE_DIR = ROOT_DIR / "pipeline"
DOCS_DIR = ROOT_DIR / "docs"
LOGS_DIR = ROOT_DIR / "logs"

ENV_NAME = "DeepVir"
ENV_FILE = PIPELINE_DIR / "workflow" / "envs" / "DeepVir.yaml"

# Regex for stripping sequencing file extensions and read pair suffixes
SAMPLE_SUFFIX_RE = re.compile(
    r"(_1|_2|_R1|_R2|_unp|_orphans)?(_001)?\.(fastq\.gz|fastq|fasta|fa|fas)$",
    re.IGNORECASE,
)
VALID_SEQ_EXTS = (".fastq.gz", ".fastq", ".fasta", ".fa", ".fas")


# =============================================================================
# PLATFORM DETECTION & PATH TRANSLATION (WINDOWS+WSL, MACOS, LINUX)
# =============================================================================
def detect_runtime_platform() -> str:
    """
    Detects the current operating system environment:
    - 'windows_wsl': Running on native Windows Python, bridging Snakemake to WSL
    - 'macos': Running natively on macOS (Darwin - Intel x86_64 or Apple Silicon arm64)
    - 'linux': Running natively on Linux or inside a WSL Linux distribution
    """
    sys_name = platform.system()
    if sys_name == "Windows":
        return "windows_wsl"
    if sys_name == "Darwin":
        return "macos"
    return "linux"


def to_engine_path(path_str: str, runtime: str) -> str:
    """
    Normalizes a file path for the target Snakemake execution engine:
    - On 'windows_wsl': converts 'D:\\DeepVir\\data' -> '/mnt/d/DeepVir/data'
    - On 'macos' / 'linux': returns standard POSIX path
    """
    if not path_str:
        return path_str
    normalized = path_str.strip().replace("\\", "/")
    if runtime == "windows_wsl":
        match = re.match(r"^([a-zA-Z]):/(.*)$", normalized)
        if match:
            drive = match.group(1).lower()
            rest = match.group(2)
            return f"/mnt/{drive}/{rest}".rstrip("/")
    return normalized


def to_host_path(path_str: str) -> Path:
    """
    Converts a path string (including '/mnt/d/...' WSL paths) into a native
    host `pathlib.Path` so Python file operations work on Windows, macOS, and Linux.
    """
    if platform.system() != "Windows":
        return Path(path_str).expanduser()
    normalized = path_str.strip().replace("\\", "/")
    match = re.match(r"^/mnt/([a-zA-Z])/(.*)$", normalized)
    if match:
        drive = match.group(1).upper()
        rest = match.group(2)
        return Path(f"{drive}:/{rest}")
    return Path(path_str)


# =============================================================================
# 1. SAMPLE DISCOVERY & SRA LIST GENERATION
# =============================================================================
def generate_sample_list(
    input_dir_str: str,
    output_dir_str: str,
    sra_file_str: str = "",
    sra_list_str: str = "",
) -> Path:
    """
    Scans `input_dir` for sequencing files (.fastq.gz, .fastq, .fasta, .fa, .fas),
    appends valid SRA accessions from `sra_file` (file path or direct IDs) and
    `sra_list_str` (direct comma/space/newline-separated SRA IDs), and writes
    `<output_dir>/samples_detected.txt`.
    """
    print("[INFO] Generating sample list (Python Executor)...", flush=True)
    input_dir = to_host_path(input_dir_str)
    output_dir = to_host_path(output_dir_str)
    output_dir.mkdir(parents=True, exist_ok=True)

    sample_list_file = output_dir / "samples_detected.txt"
    detected_samples = set()

    if input_dir.is_dir():
        for entry in input_dir.iterdir():
            if entry.is_file() and any(entry.name.lower().endswith(ext) for ext in VALID_SEQ_EXTS):
                sample_name = SAMPLE_SUFFIX_RE.sub("", entry.name)
                if sample_name:
                    detected_samples.add(sample_name)

    def _add_sra_token(raw_token: str) -> None:
        acc = raw_token.strip()
        if not acc or acc.startswith("#"):
            return
        if (input_dir / acc).is_dir():
            print(
                f"[WARNING] '{acc}' is a directory inside '{input_dir_str}' — skipping (not treated as SRA).",
                flush=True,
            )
        else:
            detected_samples.add(acc)

    if sra_file_str:
        sra_path = to_host_path(sra_file_str)
        if sra_path.is_file():
            print(f"[INFO] Appending SRA Accessions from file ({sra_path.name})...", flush=True)
            with open(sra_path, "r", encoding="utf-8") as sf:
                for raw_line in sf:
                    for token in re.split(r"[,;\s]+", raw_line.strip()):
                        _add_sra_token(token)
        else:
            # Allow passing direct accession IDs to --sra as well
            print("[INFO] Appending direct SRA Accessions from --sra...", flush=True)
            for token in re.split(r"[,;\s]+", sra_file_str.strip()):
                _add_sra_token(token)

    if sra_list_str:
        print("[INFO] Appending direct SRA Accessions (--sra-list)...", flush=True)
        for token in re.split(r"[,;\s]+", sra_list_str.strip()):
            _add_sra_token(token)

    sorted_samples = sorted(detected_samples)
    if not sorted_samples:
        print("[ERROR] FAILED: No samples found locally and no valid SRA accessions provided.", flush=True)
        sys.exit(1)

    with open(sample_list_file, "w", encoding="utf-8", newline="\n") as out_f:
        for s in sorted_samples:
            out_f.write(f"{s}\n")

    print(f"[INFO] OK ({len(sorted_samples)} unique samples: {', '.join(sorted_samples)})", flush=True)
    return sample_list_file


# =============================================================================
# 2. EXTERNAL DATABASE LINKING (INTO pipeline/resources/)
# =============================================================================
def setup_resources(diamond_db: str = "", kraken2_db: str = "", runtime: str = "linux") -> None:
    """
    Links external DIAMOND (.dmnd) and Kraken2 database directories into `pipeline/resources/`.
    Also cleans up any broken symlinks in `pipeline/resources/diamond` and `pipeline/resources/kraken2`
    so default database auto-download rules can create directories cleanly.
    """
    res_dir = PIPELINE_DIR / "resources"
    target_dia = res_dir / "diamond"
    target_k2 = res_dir / "kraken2"

    # Clean up broken symlinks for resources/diamond and resources/kraken2
    for target in (target_dia, target_k2):
        if not target.is_dir():
            print(f"[INFO] Cleaning up non-directory or broken symlink at: {target}", flush=True)
            try:
                target.unlink(missing_ok=True)
            except OSError:
                pass
            if runtime == "windows_wsl":
                subprocess.run(["wsl", "rm", "-rf", to_engine_path(str(target), runtime)], check=False)
            target.mkdir(parents=True, exist_ok=True)

    if target_dia.is_dir():
        for entry in list(target_dia.iterdir()):
            try:
                if entry.is_symlink() and not entry.exists():
                    print(f"[INFO] Removing broken symlink inside diamond/: {entry.name}", flush=True)
                    entry.unlink()
            except OSError:
                pass

    if diamond_db:
        dia_path = to_host_path(diamond_db)
        parent = dia_path.parent
        prefix = dia_path.name
        matches = list(parent.glob(f"{prefix}*")) if parent.exists() else []
        if matches:
            print("[INFO] Linking external Diamond DB files to pipeline/resources/diamond/...", flush=True)
            target_dia.mkdir(parents=True, exist_ok=True)
            for existing in target_dia.iterdir():
                if existing.is_file() or existing.is_symlink():
                    existing.unlink()
            for m in matches:
                dest = target_dia / m.name
                try:
                    dest.symlink_to(m.resolve())
                except OSError:
                    shutil.copy2(m, dest)
        else:
            print(f"[ERROR] Diamond DB files not found for prefix: {diamond_db}", flush=True)
            sys.exit(1)
    else:
        target_dia.mkdir(parents=True, exist_ok=True)

    if kraken2_db:
        k2_path = to_host_path(kraken2_db.rstrip("/\\"))
        if k2_path.is_dir():
            print("[INFO] Linking external Kraken2 DB to pipeline/resources/kraken2/...", flush=True)
            if target_k2.is_symlink() or target_k2.is_file():
                target_k2.unlink()
            elif target_k2.is_dir():
                shutil.rmtree(target_k2)
            try:
                target_k2.symlink_to(k2_path.resolve(), target_is_directory=True)
            except OSError:
                shutil.copytree(k2_path, target_k2)
        else:
            print(f"[ERROR] Kraken2 directory not found: {kraken2_db}", flush=True)
            sys.exit(1)
    else:
        if not target_k2.is_dir():
            try:
                target_k2.unlink(missing_ok=True)
            except OSError:
                pass
            target_k2.mkdir(parents=True, exist_ok=True)


# =============================================================================
# 3. PROFILE & DYNAMIC OVERRIDES GENERATOR (SAVED IN logs/)
# =============================================================================
def discover_profile_tmpdir(profile: str, user: str) -> str:
    """Reads `default-resources.tmpdir` from `pipeline/profiles/<profile>/config.yaml`."""
    candidates = [
        PIPELINE_DIR / "profiles" / profile / "config.yaml",
        PIPELINE_DIR / profile / "config.yaml",
        ROOT_DIR / "profiles" / profile / "config.yaml",
    ]
    for cfg in candidates:
        if cfg.is_file():
            try:
                content = cfg.read_text(encoding="utf-8")
                match = re.search(r"tmpdir:\s*['\"]?([^'\"\r\n#]+)['\"]?", content)
                if match:
                    val = match.group(1).strip().replace("${USER}", user).replace("$USER", user)
                    if val:
                        return val
            except Exception:
                pass
    return f"/tmp/{user}"


def resolve_profile_dir(profile: str) -> Path:
    """Resolves the Snakemake profile directory inside `pipeline/profiles/`."""
    candidates = [
        PIPELINE_DIR / "profiles" / profile,
        PIPELINE_DIR / profile,
        ROOT_DIR / "profiles" / profile,
        Path(profile),
    ]
    for c in candidates:
        if c.is_dir():
            return c.resolve()
    print(f"[ERROR] Profile directory not found: {profile}", flush=True)
    sys.exit(1)


def generate_run_overrides(
    args: argparse.Namespace,
    snakemake_workdir: str,
    temp_dir: str,
    input_posix: str,
    output_posix: str,
    logs_posix: str,
    sample_list_posix: str,
    runtime: str,
) -> Path:
    """
    Generates `logs/run_overrides.yaml` and a timestamped copy
    `logs/run_overrides_<timestamp>.yaml` in the centralized `logs/` directory.
    """
    LOGS_DIR.mkdir(parents=True, exist_ok=True)
    run_overrides = LOGS_DIR / "run_overrides.yaml"
    timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")
    backup_overrides = LOGS_DIR / f"run_overrides_{timestamp}.yaml"

    keep_download = not args.remove_download
    if args.download_only:
        keep_download = True

    lines = [
        "# Dynamic Overrides generated by DeepVir.py (Python Executor)",
        f'workdir: "{snakemake_workdir}"',
        f'output_dir: "{output_posix}"',
        f'log_dir: "{logs_posix}"',
        f'data_dir: "{input_posix}"',
        f'sample_list: "{sample_list_posix}"',
        "",
        "default-resources:",
        f'  tmpdir: "{temp_dir}"',
        "",
        "modules:",
        f"  keep_download: {str(keep_download).lower()}",
        f"  download_only: {str(args.download_only).lower()}",
        f"  host_filter: {str(args.host_filter).lower()}",
        f"  assembly: {str(args.assembly).lower()}",
        f"  kraken2: {str(args.kraken2).lower()}",
        f"  diamond: {str(args.diamond).lower()}",
        f"  darkmatter: {str(args.darkmatter).lower()}",
        f"  rvdb: {str(args.rvdb).lower()}",
        f"  reads_kraken2: {str(args.reads_kraken).lower()}",
        f"  reads_diamond: {str(args.reads_diamond).lower()}",
        "",
        "host_filter:",
        f'  accession: "{args.host_accession.strip() or "GCA_022682495.2"}"',
    ]

    if args.assembler:
        tools = [t.strip() for t in args.assembler.split(",") if t.strip()]
        if tools:
            lines.append("")
            lines.append("tool:")
            lines.append("  denovo:")
            for t in tools:
                lines.append(f"    - '{t}'")

    if args.fastp_length or args.fastp_quality:
        lines.append("")
        lines.append("fastp:")
        if args.fastp_length:
            lines.append(f"  length_required: {args.fastp_length}")
        if args.fastp_quality:
            lines.append(f"  qualified_quality_phred: {args.fastp_quality}")

    if args.rdrp_minscore is not None:
        minscore_val = int(args.rdrp_minscore) if float(args.rdrp_minscore).is_integer() else float(args.rdrp_minscore)
        lines.append("")
        lines.append("palm_annot:")
        lines.append('  seqtype: "aa"')
        lines.append(f"  minscore: {minscore_val}")
        lines.append("  minpssmscore: 10")

    yaml_text = "\n".join(lines) + "\n"
    run_overrides.write_text(yaml_text, encoding="utf-8", newline="\n")
    shutil.copy2(run_overrides, backup_overrides)

    print(f"[INFO] Generating run overrides: {to_engine_path(str(run_overrides), runtime)}", flush=True)
    print(f"[INFO] A copy will be saved at: {logs_posix}/run_overrides_{timestamp}.yaml", flush=True)
    return run_overrides


# =============================================================================
# 4. CROSS-PLATFORM CONDA & SNAKEMAKE EXECUTION ENGINE
# =============================================================================
def get_conda_shell_init() -> str:
    """
    Returns a shell snippet that locates and initializes Conda/Mamba across
    macOS (Homebrew/Miniforge/Anaconda), Linux, and WSL.
    """
    return (
        "for c in "
        "\"$HOME/mambaforge/etc/profile.d/conda.sh\" "
        "\"$HOME/miniforge3/etc/profile.d/conda.sh\" "
        "\"$HOME/miniconda3/etc/profile.d/conda.sh\" "
        "\"$HOME/anaconda3/etc/profile.d/conda.sh\" "
        "\"/opt/homebrew/Caskroom/miniforge/base/etc/profile.d/conda.sh\" "
        "\"/opt/homebrew/anaconda3/etc/profile.d/conda.sh\" "
        "\"/usr/local/Caskroom/miniforge/base/etc/profile.d/conda.sh\" "
        "\"/usr/local/bioinfo/Miniforge3-24.9.2-0/etc/profile.d/conda.sh\"; do "
        "[ -f \"$c\" ] && . \"$c\" && break; done; "
        f"conda activate {ENV_NAME} 2>/dev/null || true; "
    )


def run_snakemake_step(snakemake_args: list[str], env_vars: dict[str, str], runtime: str) -> int:
    """
    Executes Snakemake across Windows (via WSL bridge), macOS, and Linux.
    Streams stdout and stderr in real time.
    """
    quoted_args = " ".join(f'"{a}"' if " " in a else a for a in snakemake_args)
    env_exports = " ".join(f'export {k}="{v}";' for k, v in env_vars.items())
    conda_init = get_conda_shell_init()
    pipeline_posix = to_engine_path(str(PIPELINE_DIR), runtime)

    # On Apple Silicon macOS (arm64), default CONDA_SUBDIR=osx-64 if not set so Bioconda packages resolve
    if runtime == "macos" and platform.machine() == "arm64":
        env_exports = f'export CONDA_SUBDIR="${{CONDA_SUBDIR:-osx-64}}"; {env_exports}'

    bash_cmd = f'cd "{pipeline_posix}" && {env_exports} {conda_init} snakemake {quoted_args}'

    if runtime == "windows_wsl":
        cmd = ["wsl", "bash", "-c", bash_cmd]
    else:
        cmd = ["bash", "-c", bash_cmd]

    proc = subprocess.Popen(
        cmd,
        cwd=str(ROOT_DIR),
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        bufsize=1,
    )

    for line in iter(proc.stdout.readline, ""):
        print(line, end="", flush=True)
    proc.stdout.close()
    return proc.wait()


# =============================================================================
# 5. CLI PARSER & MAIN ENTRY POINT
# =============================================================================
def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="DeepVir.py",
        description="DeepVir: Viral Metagenomics & Dark Matter Discovery (Universal Python Executor)",
    )
    parser.add_argument("--input", default=str(ROOT_DIR / "data"), help="Directory containing raw reads or contigs")
    parser.add_argument("--output", default=str(ROOT_DIR / "results"), help="Directory where results will be saved")
    parser.add_argument("--sra", default="", help="Text file containing SRA Accession IDs (or direct comma-separated IDs)")
    parser.add_argument("--sra-list", default="", help="Direct comma/space-separated SRA Accession IDs (e.g. SRR10677983,SRR25008973)")
    parser.add_argument("--jobs", type=int, default=8, help="Maximum concurrent jobs/threads")
    parser.add_argument("--profile", default="local", help="Snakemake profile (local or profile_slurm)")
    parser.add_argument("--diamond_db", default="", help="External DIAMOND database (.dmnd)")
    parser.add_argument("--kraken2_db", default="", help="External Kraken2 database folder")
    parser.add_argument("--fastp-length", "--fastp-len", "--length-required", type=int, default=None)
    parser.add_argument("--fastp-quality", "--fastp-qual", "--quality-phred", type=int, default=None)
    parser.add_argument("--temp-dir", default="", help="Scratch directory for temporary files and workdir")
    parser.add_argument("--workdir", default="", help="Override Snakemake working directory")
    parser.add_argument("--conda-prefix", default="", help="Persistent directory for Snakemake Conda environments")

    # Module toggles
    parser.add_argument(
        "--host-filter",
        action="store_true",
        help="Filter out host reads via BWA against one or more NCBI host genome accessions",
    )
    parser.add_argument(
        "--host-accession",
        "--host-accessions",
        default="GCA_022682495.2,GCF_000819615.1",
        help="One or more comma-separated NCBI Genome Assembly Accessions for host/control filtering (default: GCA_022682495.2,GCF_000819615.1 for Desmodus rotundus + PhiX174)",
    )
    parser.add_argument("--assembly", action="store_true", help="Enable De Novo Assembly")
    parser.add_argument("--assembler", "--assemblers", "--denovo", default="", help="Comma-separated assemblers")
    parser.add_argument("--kraken2", action="store_true", help="Enable Kraken2 Contigs Taxonomy")
    parser.add_argument("--diamond", action="store_true", help="Enable DIAMOND Contigs Taxonomy")
    parser.add_argument("--darkmatter", action="store_true", help="Enable PalmAnnot / Dark Matter module")
    parser.add_argument(
        "--rdrp-minscore",
        "--minscore",
        dest="rdrp_minscore",
        type=float,
        default=None,
        help="Minimum RdRp score threshold (0-100) for PalmAnnot dark matter discovery (default: 75)",
    )
    parser.add_argument("--rvdb", action="store_true", help="Enable RVDB validation on dark matter ORFs")
    parser.add_argument("--download-only", action="store_true", help="Only download SRA data")
    parser.add_argument("--remove-download", action="store_true", help="Remove downloaded SRA files after run")
    parser.add_argument("--reads-kraken", action="store_true", help="Enable Kraken2 on reads")
    parser.add_argument("--reads-diamond", action="store_true", help="Enable DIAMOND + MEGAN LCA on reads")

    # Utility flags
    parser.add_argument("-n", "--dry-run", action="store_true", help="Perform Snakemake dry-run only")
    parser.add_argument("-g", "--gui", action="store_true", help="Launch Desktop GUI")
    parser.add_argument("--update-env", action="store_true", help="Update DeepVir Conda environment")
    parser.add_argument("-v", "--version", action="version", version="DeepVir v2.0 (Universal Python Executor) | 10.2026")
    return parser


def launch_gui() -> int:
    """
    Launches the Tkinter Desktop GUI (`pipeline/gui/deepvir_gui.py`).
    - If running on Linux/WSL/macOS outside the 'DeepVir' Conda environment (or if
      `tkinter` is missing from the system Python), automatically initializes Conda
      and activates 'DeepVir' before launching the GUI.
    - On Linux/WSL, if Conda installed the headless `tk-noxft` build (which disables
      TrueType/Xft font scaling and locks Xwayland to a single 6x13 fixed bitmap font),
      automatically upgrades `tk` to `conda-forge::tk=*=xft_*` once.
    """
    gui_script = PIPELINE_DIR / "gui" / "deepvir_gui.py"
    runtime = detect_runtime_platform()

    tk_xft_guard = (
        'if [ -n "$CONDA_PREFIX" ] && ls "$CONDA_PREFIX"/conda-meta/tk-*noxft*.json >/dev/null 2>&1; then '
        'echo "[INFO] Updating Conda Tk from headless noxft to xft for scalable GUI fonts (one-time setup)..."; '
        'conda install -y -c conda-forge "tk=*=xft_*" fonts-conda-ecosystem; '
        'fi; '
    )

    try:
        import tkinter  # noqa: F401
        has_tk = True
    except ImportError:
        has_tk = False

    # On Linux/WSL, check if the currently active Conda environment has tk-noxft
    conda_prefix = os.environ.get("CONDA_PREFIX", "")
    has_noxft = bool(conda_prefix and list((Path(conda_prefix) / "conda-meta").glob("tk-*noxft*.json")))

    if runtime != "windows_wsl" and (not has_tk or has_noxft or os.environ.get("CONDA_DEFAULT_ENV") != ENV_NAME):
        conda_init = get_conda_shell_init()
        gui_posix = to_engine_path(str(gui_script), runtime)
        return subprocess.call(["bash", "-c", f'{conda_init} {tk_xft_guard} python3 "{gui_posix}"'])

    if runtime == "windows_wsl" and not has_tk:
        conda_init = get_conda_shell_init()
        gui_posix = to_engine_path(str(gui_script), runtime)
        return subprocess.call(["wsl", "bash", "-c", f'{conda_init} {tk_xft_guard} python3 "{gui_posix}"'])

    return subprocess.call([sys.executable, str(gui_script)])


def main() -> int:
    parser = build_arg_parser()
    # If launched with no arguments at all, open the GUI automatically for 1-click convenience
    if len(sys.argv) == 1:
        return launch_gui()

    args = parser.parse_args()

    if args.gui:
        return launch_gui()

    runtime = detect_runtime_platform()
    LOGS_DIR.mkdir(parents=True, exist_ok=True)

    print(f"[INFO] Checking Conda environment '{ENV_NAME}'...", flush=True)
    print(
        f"       > Python Executor active (OS: {platform.system()} {platform.machine()} | Mode: {runtime} | Python {platform.python_version()})",
        flush=True,
    )

    # 1. Setup external databases if provided (and sanitize broken symlinks)
    setup_resources(args.diamond_db, args.kraken2_db, runtime)

    # 2. Generate sample list in pure Python
    sample_list_path = generate_sample_list(args.input, args.output, args.sra, args.sra_list)

    print("[INFO] Initializing DeepVir Workflow...", flush=True)

    # 3. Locate main config inside pipeline/config/config.yaml
    main_config = PIPELINE_DIR / "config" / "config.yaml"
    if not main_config.is_file():
        print(f"[ERROR] Could not find {main_config}", flush=True)
        return 1

    # 4. Resolve username, temp_dir, and snakemake_workdir
    user = os.environ.get("USER") or getpass.getuser() or "deepvir_user"
    temp_dir = args.temp_dir.strip() if args.temp_dir else discover_profile_tmpdir(args.profile, user)
    temp_dir = to_engine_path(temp_dir, runtime)
    if user not in temp_dir:
        temp_dir = f"{temp_dir.rstrip('/')}/{user}"

    snakemake_workdir = to_engine_path(args.workdir.strip(), runtime) if args.workdir else temp_dir

    if temp_dir.startswith("/tmp/") and ("slurm" in args.profile or "cluster" in args.profile):
        print(f"[WARNING] You are using a cluster profile ({args.profile}) but temp_dir is inside /tmp/.", flush=True)
        print("[WARNING] Compute nodes usually have isolated /tmp partitions. Consider passing --temp-dir on scratch.", flush=True)
        time.sleep(2)

    # Normalize paths for Snakemake config
    pipeline_posix = to_engine_path(str(PIPELINE_DIR), runtime)
    input_posix = to_engine_path(str(to_host_path(args.input).resolve()), runtime)
    output_posix = to_engine_path(str(to_host_path(args.output).resolve()), runtime)
    logs_posix = to_engine_path(str(LOGS_DIR.resolve()), runtime)
    sample_list_posix = to_engine_path(str(sample_list_path.resolve()), runtime)
    main_config_posix = f"{pipeline_posix}/config/config.yaml"
    snakefile_posix = f"{pipeline_posix}/workflow/Snakefile"

    profile_dir = resolve_profile_dir(args.profile)
    profile_posix = to_engine_path(str(profile_dir), runtime)

    # 5. Generate dynamic run_overrides.yaml inside logs/
    run_overrides_path = generate_run_overrides(
        args=args,
        snakemake_workdir=snakemake_workdir,
        temp_dir=temp_dir,
        input_posix=input_posix,
        output_posix=output_posix,
        logs_posix=logs_posix,
        sample_list_posix=sample_list_posix,
        runtime=runtime,
    )
    run_overrides_posix = to_engine_path(str(run_overrides_path), runtime)

    job_id = os.environ.get("SLURM_JOB_ID", "local")
    shadow_dir = f"{temp_dir}/deepvir_shadow/{job_id}"

    # Persistent Conda prefix: never store inside volatile /tmp unless explicitly requested,
    # because /tmp is wiped on every WSL shutdown or Linux reboot!
    if args.conda_prefix.strip():
        conda_dir = to_engine_path(args.conda_prefix.strip(), runtime)
    elif temp_dir.startswith("/tmp"):
        if runtime == "windows_wsl":
            conda_dir = f"/home/{user}/.deepvir_conda"
        else:
            home_posix = to_engine_path(str(Path.home()), runtime)
            conda_dir = f"{home_posix}/.deepvir_conda"
    else:
        conda_dir = f"{temp_dir}/deepvir_conda"

    if runtime == "windows_wsl":
        subprocess.run(["wsl", "mkdir", "-p", snakemake_workdir, shadow_dir, conda_dir, logs_posix], check=False)
    else:
        Path(snakemake_workdir).mkdir(parents=True, exist_ok=True)
        Path(shadow_dir).mkdir(parents=True, exist_ok=True)
        Path(conda_dir).mkdir(parents=True, exist_ok=True)

    print(f"[INFO] Pipeline Backend   set to: {pipeline_posix}", flush=True)
    print(f"[INFO] Centralized Logs   set to: {logs_posix}", flush=True)
    print(f"[INFO] Snakemake workdir  set to: {snakemake_workdir}", flush=True)
    print(f"[INFO] Global TMPDIR      set to: {temp_dir}", flush=True)
    print(f"[INFO] Shadow directory   set to: {shadow_dir}", flush=True)
    print(f"[INFO] Conda prefix       set to: {conda_dir}", flush=True)
    print(f"[INFO] Profile            set to: {profile_posix}", flush=True)
    print(f"[INFO] Snakefile          set to: {snakefile_posix}", flush=True)

    env_vars = {
        "TMPDIR": temp_dir,
        "TEMP": temp_dir,
        "TMP": temp_dir,
        "PYTHONUNBUFFERED": "1",
    }

    # Step 1: Unlock working directory
    print("\n> Snakemake: Unlocking working directory...", flush=True)
    rc = run_snakemake_step(
        [
            "--profile", profile_posix,
            "--snakefile", snakefile_posix,
            "--directory", snakemake_workdir,
            "--configfile", main_config_posix, run_overrides_posix,
            "--unlock", "--quiet",
        ],
        env_vars,
        runtime,
    )
    if rc != 0:
        return rc

    # Step 2: Perform dry-run
    print("\n> Snakemake: Performing a dry-run...", flush=True)
    rc = run_snakemake_step(
        [
            "--profile", profile_posix,
            "--snakefile", snakefile_posix,
            "--directory", snakemake_workdir,
            "--jobs", str(args.jobs),
            "--use-conda",
            "--conda-prefix", conda_dir,
            "--configfile", main_config_posix, run_overrides_posix,
            "--dry-run",
        ],
        env_vars,
        runtime,
    )
    print("---------------------------------------------------", flush=True)
    if rc != 0 or args.dry_run:
        if rc == 0 and args.dry_run:
            print("\n[OK] Snakemake dry-run completed successfully.", flush=True)
        return rc

    # Step 3: Main workflow execution
    print("\n> Snakemake: Starting main execution...", flush=True)
    rc = run_snakemake_step(
        [
            "--profile", profile_posix,
            "--snakefile", snakefile_posix,
            "--directory", snakemake_workdir,
            "--jobs", str(args.jobs),
            "--use-conda",
            "--conda-prefix", conda_dir,
            "--configfile", main_config_posix, run_overrides_posix,
            "--shadow-prefix", shadow_dir,
            "--keep-going",
            "--scheduler", "greedy",
        ],
        env_vars,
        runtime,
    )
    if rc != 0:
        return rc

    # Step 4: Generate HTML workflow report
    print("\n> Snakemake: Creating DAG & Report...", flush=True)
    run_snakemake_step(
        [
            "--report",
            "--profile", profile_posix,
            "--snakefile", snakefile_posix,
            "--directory", snakemake_workdir,
            "--jobs", str(args.jobs),
            "--use-conda",
            "--conda-prefix", conda_dir,
            "--configfile", main_config_posix, run_overrides_posix,
            "--keep-going",
        ],
        env_vars,
        runtime,
    )

    print("\n[OK] Done. Thank you for using DeepVir!", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
