# DeepVir: First Time Local Usage

This document guides you through the initial setup required to run DeepVir on your local machine (Linux, macOS, or Windows with WSL2), including installing Conda, preparing testing databases, and executing the pipeline using the Python orchestrator (`DeepVir.py`) or the 1-click Desktop GUI.

## 1. Getting Conda & Python

DeepVir relies on **Python 3.8+** to run the universal orchestrator (`DeepVir.py`) and **Conda** for tool dependency management, ensuring reproducibility across environments. Snakemake's native `--use-conda` flag automatically creates isolated environments for each bioinformatics tool from `pipeline/workflow/envs/`.

If you don't have Conda installed, we recommend installing **Miniforge** (which uses `conda-forge` as the default channel and is free for all usage):
1. Visit the [Miniforge GitHub repository](https://github.com/conda-forge/miniforge).
2. Download the installer for your operating system (on Windows, install Miniforge inside your WSL2 Linux distribution).
3. Follow the installation instructions provided on the repository.

Once installed, ensure `conda` and `python3` are available in your terminal:
```bash
conda --version
python3 --version
```

Create the base `DeepVir` Conda environment (which includes Snakemake) once before your first run:
```bash
conda env create -f pipeline/workflow/envs/DeepVir.yaml
```

---

## 2. Quick Start: Test Databases

If you are testing the pipeline and do not want to download massive production databases right away, you can build a lightweight *E. coli* K-12 proteome database for DIAMOND and download the RefSeq Viral database for Kraken2.

### Step 2.1: E. coli K-12 DIAMOND Database

```bash
# Create directory
mkdir -p db/diamond && cd db/diamond

# Download the E. coli K-12 proteome from UniProt
wget "https://rest.uniprot.org/uniprotkb/stream?format=fasta&query=(organism_id:83333)" -O ecoli_proteome.fasta

# Format the FASTA header to keep only the protein accession ID
awk -F'|' '/^>/ {print ">"$2} !/^>/ {print $0}' ecoli_proteome.fasta > ecoli_formatted.fasta

# Build the DIAMOND database
diamond makedb --in ecoli_formatted.fasta -d ecoli_db

# Remove temporary FASTA files
rm -rf ecoli_proteome.fasta ecoli_formatted.fasta
cd ../../
```

### Step 2.2: Kraken2 RefSeq Viral Database

For Kraken2, you can use the lightweight RefSeq Viral database provided by the [Kraken2 AWS Indexes](https://benlangmead.github.io/aws-indexes/k2).

```bash
mkdir -p db/kraken2_viral && cd db/kraken2_viral

# Download the RefSeq Viral database (approx. 0.5 GB)
wget https://genome-idx.s3.amazonaws.com/kraken/k2_viral_20240112.tar.gz

# Extract the database
tar -zxvf k2_viral_20240112.tar.gz

# Clean up
rm k2_viral_20240112.tar.gz
cd ../../
```

---

## 3. Running DeepVir Locally with Python

Once your databases are ready, you can run DeepVir either via the **Interactive Desktop GUI** or directly from the **Python CLI**.

### Step 3.1: Prepare Your Input Directory

Place your raw sequencing reads (`.fastq.gz` files) or pre-assembled contigs (`.fasta`) inside an input directory (for example, `data/`):

```bash
mkdir -p data
cp sample1_R1.fastq.gz data/
cp sample1_R2.fastq.gz data/
```

### Step 3.2 (Option A): Launch via the 1-Click Desktop GUI

You can configure all paths, databases, assemblers, and modules visually:
- **Windows**: Double-click `DeepVir_gui.bat`
- **macOS**: Double-click `DeepVir_gui.command`
- **Linux**: Run `./DeepVir_gui.sh`
- **Terminal (any OS)**:
  ```bash
  python DeepVir.py --gui
  ```

See the **[Desktop GUI Guide](desktop_gui.md)** for a complete description of every panel in the interface.

### Step 3.2 (Option B): Execute via the Python CLI

Run `DeepVir.py` from the repository root, pointing it to your input directory, output directory, local profile, and test databases:

```bash
# From the DeepVir repository root directory
python DeepVir.py \
    --input data \
    --output results \
    --profile local \
    --diamond_db db/diamond/ecoli_db.dmnd \
    --kraken2_db db/kraken2_viral \
    --reads-kraken \
    --assembly \
    --assembler spades \
    --diamond \
    --darkmatter \
    --jobs 8
```

*(Tip: On Windows, you can run `python DeepVir.py ...` directly from PowerShell or Command Prompt; `DeepVir.py` automatically bridges the execution into WSL2 and translates Windows paths like `D:\DeepVir\data` to `/mnt/d/DeepVir/data`.)*

### Explanation of Arguments:

- `--input`: Path to the directory containing your `.fastq.gz` or `.fasta` files.
- `--output`: Path where the scientific pipeline results will be saved (`results/`).
- `--profile local`: Uses the local Snakemake profile (`pipeline/profiles/local`).
- `--diamond_db`: Path to the compiled `.dmnd` database (`ecoli_db.dmnd`).
- `--kraken2_db`: Directory containing the extracted Kraken2 database (`hash.k2d`, `opts.k2d`, `taxo.k2d`).
- `--reads-kraken`: Enables read-level taxonomic classification with Kraken2.
- `--assembly --assembler spades`: Enables *de novo* assembly of reads into contigs using SPAdes.
- `--diamond`: Enables contig-level DIAMOND classification with full NCBI lineage resolution.
- `--darkmatter`: Enables the RdRp 'dark matter' viral discovery module (PalmAnnot).
- `--jobs 8`: Maximum number of parallel cores/threads Snakemake is allowed to use.
- `--dry-run` *(optional)*: Preview the Snakemake DAG and planned jobs without running the tools.

### Step 3.3: Reviewing the Results and Logs

When the pipeline finishes:
- **Scientific Results (`results/`)**:
  - `results/multiqc_report.html`: Interactive MultiQC summary report across all samples.
  - `results/<sample>/darkmatter_spades_kauto/<sample>_Report_Diversity.html`: Interactive HTML report summarizing viral diversity and RdRp candidates.
  - `results/<sample>/diamond_spades_kauto/<sample>_spades_kauto_hits_with_lineage.tsv`: Full DIAMOND output augmented with complete NCBI lineages.
  - `results/<sample>/krona_spades_kauto/`: Interactive Krona charts for exploring taxonomic composition.
- **Execution Logs (`logs/`)**:
  - `logs/<sample>/`: Individual tool execution logs for every Snakemake rule.
  - `logs/run_overrides.yaml`: Snapshot of the exact parameters used in the run.
