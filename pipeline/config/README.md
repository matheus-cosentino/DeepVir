# Configuration Guide — DeepVir

DeepVir is a modular Snakemake pipeline for viral discovery and metatranscriptomic analysis.
All parameters are defined in `config/config.yaml`.

---

## Quick Start

```bash
snakemake --configfile config/config.yaml --sdm conda --cores 8
```

---

## Parameters

### Paths

| Parameter | Default | Description |
|---|---|---|
| `output_dir` | `results` | Directory where all output files will be written |
| `data_dir` | `data` | Directory containing input FASTQ files (local samples) |
| `workdir` | `DeepVir` | Snakemake working directory. On HPC, point to a fast scratch path (e.g. `/scratch/$USER`) |
| `sample_list` | `accessions.txt` | Path to a plain text file with one sample name or SRA accession per line |

> **Tip**: SRA accessions (e.g. `SRR123456`) are automatically downloaded by DeepVir.
> Local samples must already be in `data_dir` following the pattern `{sample}_R1.fastq.gz` / `{sample}_R2.fastq.gz`.

---

### Analysis Modules

Activate or deactivate analysis modules under the `modules:` key.

| Module | Default | Description |
|---|---|---|
| `keep_download` | `true` | Keep downloaded SRA FASTQ files after run (set `false` to delete them on success) |
| `download_only` | `false` | Only download data from SRA, skip all analyses |
| `assembly` | `true` | Perform de novo assembly |
| `kraken2` | `true` | Taxonomic classification with Kraken2 (contig-level) |
| `diamond` | `true` | Taxonomic annotation with DIAMOND BLAST (contig-level) |
| `darkmatter` | `true` | Dark matter search using palmDB / PalmAnnot |
| `rvdb` | `false` | Optional RVDB validation on dark matter ORFs |
| `reads_kraken2` | `true` | Read-level taxonomic classification with Kraken2 |
| `reads_diamond` | `true` | Read-level taxonomic classification with DIAMOND + MEGAN LCA |

---

### De novo Assembly Tool

Choose one or more assemblers under `tool.denovo`:

| Option | Recommended for |
|---|---|
| `spades` | Illumina short reads (RNA-seq, metatranscriptome) — **default** |
| `megahit` | Illumina short reads (faster, less memory) |
| `flye` | Oxford Nanopore long reads |
| `raven` | Oxford Nanopore long reads (faster) |
| `medaka_flye` | ONP Flye assembly + Medaka polishing |
| `medaka_raven` | ONP Raven assembly + Medaka polishing |

```yaml
tool:
  denovo:
    - 'spades'
```

---

### Tool Parameters

#### `fastp` — Read Quality Control

| Parameter | Default | Description |
|---|---|---|
| `length_required` | `50` | Minimum read length after trimming (use `500` for Nanopore) |
| `qualified_quality_phred` | `30` | Minimum Phred quality score (use `10` for long reads) |

#### `spades` — SPAdes Assembler (Illumina)

| Parameter | Default | Description |
|---|---|---|
| `algorithm` | `--rnaviral` | SPAdes mode. Options: `--meta`, `--rna`, `--rnaviral`, `--metaviral`, `--plasmid`, etc. |
| `kmer` | `auto` | K-mer sizes (odd integers ≤ 128). Use `auto` for SPAdes auto-detection. |

#### `flye` — Flye Assembler (Nanopore / PacBio)

| Parameter | Default | Description |
|---|---|---|
| `type` | `nano-corr` | Read type. Options: `nano-raw`, `nano-corr`, `nano-hq`, `pacbio-raw`, `pacbio-corr`, `pacbio-hifi` |

#### `diamond` — DIAMOND BLAST

| Parameter | Default | Description |
|---|---|---|
| `min_contig_len` | `600` | Minimum contig length (bp) to submit to DIAMOND |
| `outfmt` | `6` | DIAMOND output format (6 = BLAST tabular) |
| `max_target_seqs` | `1` | Number of top hits to report per query |
| `evalue` | `0.001` | Maximum e-value threshold |
| `sensitivity` | `--sensitive` | Search sensitivity: `--sensitive`, `--more-sensitive`, `--very-sensitive`, `--ultra-sensitive` |
| `block_size` | `6` | RAM per block in GB (~6–8 GB per unit). Adjust based on available memory. |
| `index_chunks` | `1` | Number of index chunks. `1` = full index in RAM (fastest) |

#### `kraken2` — Kraken2 Classifier

| Parameter | Default | Description |
|---|---|---|
| `confidence` | `0.0` | Confidence score threshold (0.0–1.0). Higher = more conservative |

#### `kraken_biom` — Kraken2 BIOM Export

| Parameter | Default | Description |
|---|---|---|
| `format` | `hdf5` | BIOM output format: `hdf5`, `json`, or `tsv` |
| `max` | `D` | Highest taxonomic rank to include (D=Domain, P=Phylum, ...) |
| `min` | `F` | Lowest taxonomic rank to include (F=Family, G=Genus, S=Species, ...) |

#### `palm_annot` — Dark Matter / PalmAnnot

| Parameter | Default | Description |
|---|---|---|
| `seqtype` | `aa` | Sequence type: `aa` (amino acid) |
| `minscore` | `75` | Minimum palmDB score |
| `minpssmscore` | `10` | Minimum PSSM score |

#### `hmmscan` — HMMER (dark matter validation)

| Parameter | Default | Description |
|---|---|---|
| `evalue` | `1e-5` | Maximum e-value for HMMER hits |

#### `medaka_model` — Medaka polishing (Nanopore only)

| Parameter | Default | Description |
|---|---|---|
| `medaka_model` | `r941_min_sup_g507` | Basecalling model. Use R9.4.1 (`r941_*`) or R10 (`r1041_*`) models matching your flowcell |

---

### Resource Paths (Databases)

All database directories are defined under `resources:`. Paths are relative to the repository root.

| Key | Default path | Description |
|---|---|---|
| `diamond` | `resources/diamond/` | DIAMOND database (NR or custom protein DB) |
| `kraken2` | `resources/kraken2/` | Kraken2 database directory |
| `taxonomy` | `resources/taxonomy/` | NCBI taxonomy files (`prot.accession2taxid`, `nodes.dmp`, `names.dmp`) |
| `krona` | `resources/krona/taxonomy` | Krona taxonomy database |
| `megan` | `resources/megan/` | MEGAN mapping files for LCA |
| `palm_annot_dir` | `resources/palm_annot/` | PalmDB files for dark matter search |
| `rvdb` | `resources/rvdb/` | RVDB HMM profiles (optional) |
| `basta_db` | `resources/basta_db/` | BASTA LCA database (optional) |
| `logo_dirs` | `resources/logo/` | Logo assets for HTML reports |

---

## Sample List Format

The `accessions.txt` file (or any file pointed to by `sample_list`) must contain one entry per line:

```
SRR12345678       # SRA accession — downloaded automatically
MySample01        # Local sample — must exist as data/MySample01_R1.fastq.gz
```

---

## HPC / SLURM Usage

For cluster execution, use a Snakemake profile. A SLURM profile is provided under `profiles/`:

```bash
snakemake --configfile config/config.yaml \
          --sdm conda \
          --profile profiles/slurm \
          --jobs 50
```

See `slurm_scripts/` for example submission scripts.
