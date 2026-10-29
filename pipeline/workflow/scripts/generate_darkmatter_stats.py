#!/usr/bin/env python3
"""
Generate Dark Matter statistics for MultiQC
Captures no-hit contigs and RdRp identification
"""

from Bio import SeqIO
import json
import sys
import os
import gzip

def main():
    try:
        nohit_fasta = snakemake.input.nohits
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

    log = open(log_file, 'w')
    
    try:
        # Parse no-hit FASTA
        num_nohits = 0
        total_bp_nohits = 0
        lengths = []
        
        if os.path.exists(nohit_fasta) and os.path.getsize(nohit_fasta) > 0:
            for record in SeqIO.parse(nohit_fasta, "fasta"):
                num_nohits += 1
                seq_len = len(record.seq)
                total_bp_nohits += seq_len
                lengths.append(seq_len)
        
        # Calculate additional statistics
        if lengths:
            lengths.sort(reverse=True)
            mean_len = total_bp_nohits / len(lengths)
            median_len = lengths[len(lengths)//2]
            n50 = calculate_n50(lengths)
            max_len = lengths[0]
            min_len = lengths[-1]
        else:
            mean_len = 0
            median_len = 0
            n50 = 0
            max_len = 0
            min_len = 0
        
        # Create statistics dictionary
        stats_dict = {
            "sample": sample_name,
            "tool": tool_name,
            "no_hit_contigs": int(num_nohits),
            "no_hit_bp": int(total_bp_nohits),
            "mean_no_hit_length": round(mean_len, 2),
            "median_no_hit_length": int(median_len),
            "no_hit_n50": int(n50),
            "max_no_hit_length": int(max_len),
            "min_no_hit_length": int(min_len)
        }
        
        # Write JSON
        with open(output_json, 'w') as f:
            json.dump(stats_dict, f, indent=2)
        
        # Write MultiQC-friendly TSV table
        with open(output_tsv, 'w') as f:
            f.write("sample\ttool\tno_hit_contigs\tno_hit_bp\tmean_no_hit_length\tmedian_no_hit_length\tno_hit_n50\tmax_no_hit_length\tmin_no_hit_length\n")
            f.write(f"{sample_name}\t{tool_name}\t{int(num_nohits)}\t{int(total_bp_nohits)}\t{mean_len:.2f}\t{int(median_len)}\t{int(n50)}\t{int(max_len)}\t{int(min_len)}\n")
        
        log.write(f"Successfully generated Dark Matter statistics for {sample_name} ({tool_name})\n")
        log.write(f"No-hit contigs: {num_nohits}, Total BP: {total_bp_nohits}, N50: {n50}\n")
        
    except Exception as e:
        log.write(f"ERROR: {str(e)}\n")
        import traceback
        log.write(traceback.format_exc())
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
