import sys
import os
import gzip
from Bio import SeqIO

def main():
    try:
        reads_input = snakemake.input.reads
        reads_file = str(reads_input[0]) if isinstance(reads_input, list) else str(reads_input)

        blast_input = snakemake.input.diamond
        blast_file = str(blast_input[0]) if isinstance(blast_input, list) else str(blast_input)

        output_file = snakemake.output.nohits
        sample_name = snakemake.wildcards.sample
        log_file = snakemake.log[0]
        min_len = int(getattr(snakemake.params, "min_len", 500))
    except NameError:
        print("ERROR: Designed to be run via Snakemake.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"ERROR accessing Snakemake variables: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"Starting sample {sample_name}: extracting Nanopore dark reads (no DIAMOND hit, len >= {min_len})", file=sys.stderr)
    print(f"  Input Reads: {reads_file}", file=sys.stderr)
    print(f"  DIAMOND Hits: {blast_file}", file=sys.stderr)
    print(f"  Output FASTA: {output_file}", file=sys.stderr)

    ids_to_exclude = set()
    if os.path.exists(blast_file) and os.path.getsize(blast_file) > 0:
        try:
            with open(blast_file, "r", encoding="utf-8", errors="ignore") as f:
                for line in f:
                    line = line.strip()
                    if not line:
                        continue
                    parts = line.split("\t")
                    if len(parts) > 1 and "NOT_FOUND" not in parts[1]:
                        # Query ID can contain read name before space/slash
                        q_id = parts[0].split()[0]
                        ids_to_exclude.add(q_id)
        except Exception as e:
            print(f"WARNING: Error reading DIAMOND file {blast_file}: {e}", file=sys.stderr)

    print(f"Found {len(ids_to_exclude)} read IDs with taxonomic hits to exclude.", file=sys.stderr)

    count_kept = 0
    count_total = 0

    if not os.path.exists(reads_file) or os.path.getsize(reads_file) == 0:
        print(f"[INFO] Empty input reads file. Creating empty {output_file}", file=sys.stderr)
        open(output_file, "w").close()
        return

    # Determine file format and compression
    is_gz = reads_file.endswith(".gz")
    is_fastq = any(ext in reads_file.lower() for ext in [".fastq", ".fq"])
    file_format = "fastq" if is_fastq else "fasta"
    open_func = gzip.open if is_gz else open

    with open_func(reads_file, "rt", encoding="utf-8", errors="ignore") as in_handle, open(output_file, "w", encoding="utf-8") as out_handle:
        for record in SeqIO.parse(in_handle, file_format):
            count_total += 1
            rec_id = record.id.split()[0]
            if rec_id not in ids_to_exclude and len(record.seq) >= min_len:
                count_kept += 1
                # Write as clean FASTA
                out_handle.write(f">{record.id}\n{str(record.seq)}\n")

    print(f"Processed {count_total} total reads. Kept {count_kept} dark reads ({count_kept/max(1,count_total)*100:.1f}%).", file=sys.stderr)

if __name__ == "__main__":
    main()

