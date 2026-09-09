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

        # Read Diamond report
        # Format: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore [extra...]
        # outfmt 6 has NO header — we must specify header=None to avoid the first
        # data row being silently consumed as column names.
        BASE_COLS = ["qseqid","sseqid","pident","length","mismatch",
                     "gapopen","qstart","qend","sstart","send","evalue","bitscore"]
        try:
            df = pd.read_csv(diamond_file, sep='\t', comment='#', header=None)
            if len(df.columns) >= len(BASE_COLS):
                # Assign base names; any extra columns get generic names
                extra = [f"col{i}" for i in range(len(BASE_COLS), len(df.columns))]
                df.columns = BASE_COLS + extra
            else:
                df.columns = BASE_COLS[:len(df.columns)]
        except pd.errors.EmptyDataError:
            # No hits in the diamond report — file is empty or only contains comment lines
            log.write(f"WARNING: Diamond report is empty (no hits) for {diamond_file}. Writing zero-hit stats.\n")
            df = pd.DataFrame(columns=BASE_COLS)
        
        # Count unique queries with hits
        num_queries_with_hits = len(df['qseqid'].unique()) if 'qseqid' in df.columns and len(df) > 0 else 0
        total_hits = len(df)
        
        # Calculate statistics
        if len(df) > 0:
            mean_identity = df['pident'].mean() if 'pident' in df.columns else 0
            mean_evalue = df['evalue'].mean() if 'evalue' in df.columns else 0
            mean_bitscore = df['bitscore'].mean() if 'bitscore' in df.columns else 0
        else:
            mean_identity = 0
            mean_evalue = 0
            mean_bitscore = 0
        
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
