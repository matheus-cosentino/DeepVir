#!/usr/bin/env python3
"""
Generate assembly statistics for MultiQC
Creates a TSV file with contig statistics
"""

from Bio import SeqIO
import json
import sys
import gzip

def main():
    try:
        fasta_file = snakemake.input.fasta
        output_json = snakemake.output.json
        output_tsv = snakemake.output.tsv
        sample_name = snakemake.wildcards.sample
        tool_name = snakemake.wildcards.tool
        log_file = snakemake.log[0]
    except NameError:
        print("ERROR: This script is designed to be run via Snakemake.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"ERROR: Could not access Snakemake variables: {e}", file=sys.stderr)
        sys.exit(1)

    # Configure logging
    log = open(log_file, 'w')
    
    try:
        # Parse FASTA and calculate statistics
        sequences = []
        for record in SeqIO.parse(fasta_file, "fasta"):
            sequences.append(len(record.seq))
        
        if not sequences:
            log.write(f"WARNING: No sequences found in {fasta_file}\n")
            sequences = [0]
        
        sequences.sort(reverse=True)
        num_contigs = len(sequences)
        total_bp = sum(sequences)
        mean_len = total_bp / num_contigs if num_contigs > 0 else 0
        median_len = sequences[num_contigs // 2] if num_contigs > 0 else 0
        n50 = calculate_n50(sequences)
        max_len = sequences[0] if sequences else 0
        min_len = sequences[-1] if sequences else 0
        
        # Create JSON output for MultiQC
        stats_dict = {
            "sample": sample_name,
            "tool": tool_name,
            "num_contigs": num_contigs,
            "total_bp": int(total_bp),
            "mean_length": round(mean_len, 2),
            "median_length": int(median_len),
            "n50": int(n50),
            "max_length": int(max_len),
            "min_length": int(min_len)
        }
        
        # Write JSON
        with open(output_json, 'w') as f:
            json.dump(stats_dict, f, indent=2)
        
        # Write MultiQC-friendly TSV table
        with open(output_tsv, 'w') as f:
            f.write("sample\ttool\tnum_contigs\ttotal_bp\tmean_length\tmedian_length\tn50\tmax_length\tmin_length\n")
            f.write(f"{sample_name}\t{tool_name}\t{num_contigs}\t{int(total_bp)}\t{mean_len:.2f}\t{int(median_len)}\t{int(n50)}\t{int(max_len)}\t{int(min_len)}\n")
        
        log.write(f"Successfully generated assembly statistics for {sample_name} ({tool_name})\n")
        log.write(f"Contigs: {num_contigs}, Total BP: {total_bp}, N50: {n50}\n")
        
    except Exception as e:
        log.write(f"ERROR: {str(e)}\n")
        log.close()
        sys.exit(1)
    
    log.close()

def calculate_n50(seq_lengths):
    """Calculate N50 from a list of sequence lengths"""
    if not seq_lengths:
        return 0
    sorted_lengths = sorted(seq_lengths, reverse=True)
    total = sum(sorted_lengths)
    target = total / 2
    cumsum = 0
    for length in sorted_lengths:
        cumsum += length
        if cumsum >= target:
            return length
    return 0

if __name__ == "__main__":
    main()
