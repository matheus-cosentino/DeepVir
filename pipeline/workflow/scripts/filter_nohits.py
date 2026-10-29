from Bio import SeqIO
import sys
import gzip

def main():
    # Access snakemake variables
    try:
        # Garante tratamento correto se input for lista ou string
        blast_input = snakemake.input.diamond
        blast_file = str(blast_input[0]) if isinstance(blast_input, list) else str(blast_input)
        
        fasta_input = snakemake.input.fasta
        fasta_file = str(fasta_input[0]) if isinstance(fasta_input, list) else str(fasta_input)
        
        output_file = snakemake.output.nohits
        sample_name = snakemake.wildcards.sample
        log_file = snakemake.log[0]
    except NameError:
        print("ERROR: This script is designed to be run via Snakemake.", file=sys.stderr)
        print("It expects 'snakemake.input', 'snakemake.output', etc.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"ERROR: Could not access Snakemake variables: {e}", file=sys.stderr)
        sys.exit(1)

    # Configure logging to stderr (which snakemake redirects to the log file)
    print(f"Starting sample {sample_name}: filtering sequences with taxonomic hits", file=sys.stderr)
    print(f"  Log file: {log_file}", file=sys.stderr)
    print(f"  Input FASTA: {fasta_file}", file=sys.stderr)
    print(f"  DIAMOND results: {blast_file}", file=sys.stderr)
    print(f"  Output FASTA: {output_file}", file=sys.stderr)
 
    # Identify sequence IDs that have a 'proper' taxonomic hit and should be excluded.
    ids_to_exclude = set()
    try:
        with open(blast_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                if line.strip():
                    parts = line.split('\t')
                    # Exclude if there's a hit and it's not a 'NOT_FOUND' hit
                    if len(parts) > 1 and "NOT_FOUND" not in parts[1]:
                        ids_to_exclude.add(parts[0])
    except FileNotFoundError:
        print(f"ERROR: DIAMOND hits file not found at {blast_file}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"ERROR reading DIAMOND file {blast_file}: {e}", file=sys.stderr)
        sys.exit(1)
    
    print(f"Found {len(ids_to_exclude)} unique sequences with taxonomic hits to be excluded.", file=sys.stderr)
    
    # Process FASTA file: keep sequences that are not in the exclusion list.
    # This includes sequences with no hits at all, or hits to 'NOT_FOUND'.
    count_kept = 0
    count_total = 0
    
    try:
        # Handle gzipped input if needed
        open_func = open
        read_mode = 'rt'
        if fasta_file.endswith('.gz'):
            open_func = gzip.open
            print("Detected gzipped FASTA input", file=sys.stderr)
            
        with open_func(fasta_file, read_mode) as fasta_handle, open(output_file, 'w') as output_handle:
            for record in SeqIO.parse(fasta_handle, "fasta"):
                count_total += 1
                if record.id not in ids_to_exclude:
                    count_kept += 1
                    SeqIO.write(record, output_handle, "fasta")
                    
    except FileNotFoundError:
        print(f"ERROR: FASTA file not found at {fasta_file}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"ERROR processing FASTA file {fasta_file}: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Summary statistics
    print(f"Processed {count_total} total sequences", file=sys.stderr)
    percent_kept = (count_kept / count_total * 100) if count_total > 0 else 0
    print(f"Saved {count_kept} sequences with no taxonomic hit ({percent_kept:.1f}%)", file=sys.stderr)
    print(f"Output written to: {output_file}", file=sys.stderr)
    print("Filtering complete.", file=sys.stderr)

if __name__ == "__main__":
    main()