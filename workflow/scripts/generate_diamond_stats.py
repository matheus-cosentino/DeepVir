#!/usr/bin/env python3
"""
Generate Diamond hit statistics for MultiQC
Extracts hit counts and coverage information
"""

import pandas as pd
import json
import sys

def main():
    try:
        diamond_file = snakemake.input.diamond
        output_json = snakemake.output.json
        output_tsv = snakemake.output.tsv
        sample_name = snakemake.wildcards.sample
        tool_name = getattr(snakemake.wildcards, 'tool', 'reads')
        log_file = snakemake.log[0]
    except NameError:
        print("ERROR: This script is designed to be run via Snakemake.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"ERROR: Could not access Snakemake variables: {e}", file=sys.stderr)
        sys.exit(1)

    log = open(log_file, 'w')
    
    try:
        wildcards = snakemake.wildcards
        sample_name = wildcards.sample
        tool_name = getattr(wildcards, 'tool', 'reads')

        # Read Diamond report in a streaming manner to avoid loading
        # massive read alignment files into memory with pandas.
        total_hits = 0
        sum_identity = 0.0
        sum_evalue = 0.0
        sum_bitscore = 0.0
        unique_queries = set()

        try:
            with open(diamond_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue
                    parts = line.split('\t')
                    if len(parts) >= 12:
                        total_hits += 1
                        unique_queries.add(parts[0])
                        try:
                            sum_identity += float(parts[2])
                            sum_evalue += float(parts[10])
                            sum_bitscore += float(parts[11])
                        except ValueError:
                            pass
        except FileNotFoundError:
            log.write(f"ERROR: Diamond report file not found: {diamond_file}\n")
            sys.exit(1)

        num_queries_with_hits = len(unique_queries)
        if total_hits > 0:
            mean_identity = sum_identity / total_hits
            mean_evalue = sum_evalue / total_hits
            mean_bitscore = sum_bitscore / total_hits
        else:
            mean_identity = 0.0
            mean_evalue = 0.0
            mean_bitscore = 0.0
            log.write(f"WARNING: Diamond report is empty (no hits) for {diamond_file}. Writing zero-hit stats.\n")
        
        # Create statistics dictionary
        stats_dict = {
            "sample": sample_name,
            "tool": tool_name,
            "total_hits": int(total_hits),
            "unique_queries": int(num_queries_with_hits),
            "mean_identity": round(mean_identity, 2),
            "mean_evalue": float(mean_evalue),
            "mean_bitscore": round(mean_bitscore, 2)
        }
        
        # Write JSON
        with open(output_json, 'w') as f:
            json.dump(stats_dict, f, indent=2)
        
        # Write MultiQC-friendly TSV table
        with open(output_tsv, 'w') as f:
            f.write("sample\ttool\ttotal_hits\tunique_queries\tmean_identity\tmean_evalue\tmean_bitscore\n")
            f.write(f"{sample_name}\t{tool_name}\t{int(total_hits)}\t{int(num_queries_with_hits)}\t{mean_identity:.2f}\t{mean_evalue:.2e}\t{mean_bitscore:.2f}\n")
        
        log.write(f"Successfully generated Diamond statistics for {sample_name} ({tool_name})\n")
        log.write(f"Total hits: {total_hits}, Unique queries: {num_queries_with_hits}\n")
        
    except Exception as e:
        log.write(f"ERROR: {str(e)}\n")
        import traceback
        log.write(traceback.format_exc())
        log.close()
        sys.exit(1)
    
    log.close()

if __name__ == "__main__":
    main()
