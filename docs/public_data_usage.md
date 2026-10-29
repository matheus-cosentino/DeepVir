# Working with Public Data (NCBI SRA)

One of the key advantages of the DeepVir pipeline is the ability to seamlessly incorporate publicly available sequencing datasets from the [NCBI Sequence Read Archive (SRA)](https://www.ncbi.nlm.nih.gov/sra/) directly into your viral discovery workflow.

## 1. Prepare the Accessions File

Create a plain-text file (for example, `data/accessions.txt`) where each line corresponds to an SRA Run identifier (`SRR...`, `ERR...`, or `DRR...`):

```text
SRR25008973
SRR25008974
SRR25008975
SRR25008976
SRR25008977
```

## 2. Download and Analyze SRA Samples with Python

Pass the accessions file to `DeepVir.py` using the `--sra` flag along with your desired analysis modules:

```bash
python DeepVir.py \
    --input data \
    --output results \
    --sra data/accessions.txt \
    --profile local \
    --assembly \
    --assembler spades \
    --diamond \
    --kraken2 \
    --darkmatter \
    --jobs 8
```

## 3. Download Only Mode (`--download-only`)

If you want to download and compress the SRA FASTQ files into your `--input` directory first and run the analysis later, add the `--download-only` flag:

```bash
python DeepVir.py \
    --input data \
    --output results \
    --sra data/accessions.txt \
    --download-only \
    --jobs 4
```

## 4. Automatic Cleanup of Raw SRA Files (`--remove-download`)

By default, downloaded `.fastq.gz` files are kept inside `<input>/sra/`. If disk space is limited and you want Snakemake to automatically delete the raw SRA reads once the workflow finishes, add `--remove-download`:

```bash
python DeepVir.py \
    --input data \
    --output results \
    --sra data/accessions.txt \
    --assembly --diamond --darkmatter \
    --remove-download
```
