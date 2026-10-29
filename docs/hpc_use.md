# DeepVir HPC (Slurm) Usage Guide

This document explains how to execute the DeepVir pipeline in a High-Performance Computing (HPC) environment managed by the **Slurm** workload manager using the Python executor (`DeepVir.py`).

## Resource Management

When running on an HPC cluster, it is crucial to use a Snakemake profile so that each pipeline rule receives the appropriate CPU, memory, and walltime limits. DeepVir provides a pre-configured Slurm profile located at `pipeline/profiles/profile_slurm`.

**The `--profile profile_slurm` flag is key** to running DeepVir on an HPC cluster. `DeepVir.py` automatically resolves profile names under `pipeline/profiles/` and delegates job submission to Slurm, scheduling partition, CPU, and memory requirements per rule.

## Basic Execution (Assembly + DIAMOND + Dark Matter)

A typical invocation for the complete viral discovery pipeline on an HPC login or head node is:

```bash
python3 DeepVir.py \
    --input /path/to/raw_reads \
    --output /path/to/results \
    --profile profile_slurm \
    --assembly \
    --assembler spades \
    --diamond \
    --darkmatter \
    --diamond_db /path/to/db/diamond/nr.dmnd \
    --kraken2_db /path/to/db/kraken2 \
    --jobs 20
```

### Explanation of Arguments:
- `--profile profile_slurm`: Instructs Snakemake to use the Slurm executor defined in `pipeline/profiles/profile_slurm`.
- `--jobs 20`: Maximum number of parallel jobs that Snakemake will submit to the Slurm queue simultaneously.
- `--assembly --assembler spades`: Enables *de novo* assembly with SPAdes.
- `--diamond` and `--darkmatter`: Enables contig taxonomic annotation and RdRp viral discovery.
- `--diamond_db` and `--kraken2_db`: Paths to the shared cluster databases.

## Optional RVDB Validation

To include optional RVDB HMM validation on top of the Dark Matter analysis, append the `--rvdb` flag:

```bash
python3 DeepVir.py \
    --input /path/to/raw_reads \
    --output /path/to/results \
    --profile profile_slurm \
    --assembly \
    --diamond \
    --darkmatter \
    --rvdb \
    --diamond_db /path/to/db/diamond/nr.dmnd \
    --kraken2_db /path/to/db/kraken2 \
    --jobs 20
```

## Quick Read-Level Profiling

If you only want a rapid initial evaluation of your samples directly on the quality-filtered reads without running resource-intensive assemblies:

```bash
python3 DeepVir.py \
    --input /path/to/raw_reads \
    --output /path/to/results \
    --profile profile_slurm \
    --reads-kraken \
    --reads-diamond \
    --kraken2_db /path/to/db/kraken2 \
    --diamond_db /path/to/db/diamond/nr.dmnd \
    --jobs 20
```

## Previewing the DAG Before Submitting (`--dry-run`)

Before submitting jobs to the Slurm queue, you can validate your inputs and preview the planned rule graph with `--dry-run`:

```bash
python3 DeepVir.py \
    --input /path/to/raw_reads \
    --output /path/to/results \
    --profile profile_slurm \
    --assembly --diamond --darkmatter \
    --dry-run
```

For more advanced configuration of execution parameters (e.g., `fastp` filters, K-mers, and long-read assemblers), see the [Advanced Configuration Tutorial](advanced_tutorial.md).
