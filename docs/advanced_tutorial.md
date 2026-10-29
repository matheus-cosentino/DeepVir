# Advanced Configuration Tutorial

This page gathers the main options for tuning preprocessing with `fastp`, selecting the contig assembler, and customizing downstream parameters in DeepVir.

## 1. Where to Edit Configuration

DeepVir offers two ways to customize pipeline parameters:
1. **Via the Python CLI (`DeepVir.py`) or Desktop GUI (`--gui`)**: Frequently changed parameters such as `fastp` thresholds (`--fastp-length`, `--fastp-quality`) and the *de novo* assembler engine (`--assembler`) can be passed directly on the command line or selected in the GUI, which writes them into `logs/run_overrides.yaml`.
2. **Via the Static Configuration File**: Global technical parameters (such as SPAdes K-mers, Flye read types, DIAMOND E-values, and PalmAnnot score cutoffs) are stored in [pipeline/config/config.yaml](../pipeline/config/config.yaml).

After editing [pipeline/config/config.yaml](../pipeline/config/config.yaml), save the file and run `python DeepVir.py` again—the workflow will automatically read the updated values.

---

## 2. Tuning `fastp` Quality Control

`fastp` filters and cleans raw reads before assembly and taxonomic analysis.

### Setting via Python CLI

You can override `fastp` thresholds directly when calling `DeepVir.py`:

```bash
python DeepVir.py \
    --input data \
    --output results \
    --fastp-length 50 \
    --fastp-quality 30 \
    --assembly --diamond
```

### Setting in `pipeline/config/config.yaml`

In the `fastp` section of [pipeline/config/config.yaml](../pipeline/config/config.yaml), you can adjust the default values:

- `length_required`: minimum read length after trimming.
- `qualified_quality_phred`: minimum Phred quality threshold.

Example:

```yaml
fastp:
  length_required:
    - 50
  qualified_quality_phred:
    - 30
```

### Practical Recommendations

- **Standard Illumina data**: `length_required: 30`–`50` and `qualified_quality_phred: 30`.
- **Lower-quality or degraded RNA samples**: try `20` or `25` for Phred quality and `30` for minimum length.
- **Nanopore / PacBio long-read data**: use a larger minimum length (`500`) and a lower Phred threshold (`10`). *(In the Desktop GUI, click the **Nanopore / Long-read (L: 500, Q: 10)** preset button.)*

---

## 3. Choosing the Contig Assembler

You can select the assembler directly from the CLI with `--assembler`:

```bash
python DeepVir.py \
    --input data \
    --output results \
    --assembly \
    --assembler megahit \
    --diamond --darkmatter
```

Available `--assembler` choices:
- `spades` *(default)*: Best option for Illumina metagenomic and metatranscriptomic reads (`--rnaviral` / `--meta`).
- `megahit`: Faster and lower memory footprint when you prioritize speed on large short-read datasets.
- `flye`: Recommended for long-read data (Nanopore / PacBio).
- `raven`: Fast long-read assembler for Nanopore datasets.
- `medaka_flye`: Flye assembly followed by Medaka consensus polishing for Oxford Nanopore reads.
- `medaka_raven`: Raven assembly followed by Medaka consensus polishing for Oxford Nanopore reads.

You can also change the default assembler list in the `tool.denovo` section of [pipeline/config/config.yaml](../pipeline/config/config.yaml):

```yaml
tool:
  denovo:
    - 'spades'
    # - 'megahit'
    # - 'flye'
```

---

## 4. Tuning SPAdes

If you use SPAdes, two settings in [pipeline/config/config.yaml](../pipeline/config/config.yaml) control its behavior:

### 4.1 Algorithm Mode

In the `spades.algorithm` section:

```yaml
spades:
  algorithm:
    - '--rnaviral'
```

- `--rnaviral` or `--meta`: recommended for viral metatranscriptomic and metagenomic datasets.
- `--isolate`: suitable for high-coverage viral isolate sequencing.

### 4.2 K-mer Sizes

The `spades.kmer` section defines the k-mer sizes used by SPAdes:

```yaml
spades:
  kmer:
    - 'auto'
```

- `auto`: automatic k-mer selection based on read length (recommended starting point).
- Explicit values like `21`, `33`, `55` or combinations such as `21_33_55` can be configured for targeted evaluations.

---

## 5. Long-Read Assembly Settings

When working with long-read data (`--assembler flye`, `raven`, `medaka_flye`, or `medaka_raven`), you can configure the Flye read type and Medaka model in [pipeline/config/config.yaml](../pipeline/config/config.yaml):

```yaml
flye:
  type: "nano-corr"

medaka:
  model: "r941_min_hac_g507"
```

Supported `flye.type` values:
- `pacbio-raw`
- `pacbio-corr`
- `pacbio-hifi`
- `nano-raw`
- `nano-corr`
- `nano-hq`

---

## 6. Adjusting Downstream Parameters

Other key parameters in [pipeline/config/config.yaml](../pipeline/config/config.yaml) include:

- `diamond.min_contig_len`: minimum contig length (bp) retained before DIAMOND annotation (default: `250`).
- `diamond.evalue`: E-value threshold for DIAMOND alignments (default: `1e-3`).
- `kraken2.confidence`: classification confidence threshold for Kraken2 (default: `0.2`).
- `palm_annot.minscore`: minimum catalytic motif score threshold for RdRp candidate detection in the Dark Matter module (default: `50`).

---

## 7. Example Run After Editing Configuration

```bash
python DeepVir.py \
    --input data \
    --output results \
    --profile local \
    --assembly \
    --assembler spades \
    --diamond \
    --darkmatter \
    --jobs 8
```

If you are running on an HPC cluster, change `--profile local` to `--profile profile_slurm`.
