###################################################################################
#                                Snakefile                                        #
#                         MSc. Matheus Cosentino                                  # 
###################################################################################
#                                                                                 #
#         ████  █████ █████ ████  █   █ ███ ████                                  #
#         █░░░█ █░░░░░█░░░░░█░░░█ █░  █░ █░░█░░░█                                 #
#         █░░░█░████░░████░░████░░█░░ █░░█░░████░░                                #
#         █░░ █░█░░░░ █░░░░ █░░░░ ░█░█ ░░█░░█░░█░ ░                               #
#         ████ ░█████░█████░█░░░░░  █ ░ ███░█░░░█░                                #
#          ░░░░ ░░░░░░ ░░░░░ ░░      ░ ░ ░░░ ░░  ░                                 #
#           ░░░░  ░░░░░ ░░░░░ ░       ░   ░░░ ░   ░                                #
#                                                                                 #
###################################################################################
#                         version: 2.0  |  10.2026                              #
###################################################################################

###########################################
# --- 1. Import libraryes to be used --- #
##########################################

import os, re, glob, time, sys, subprocess, platform, yaml
import urllib.request
from snakemake.io import expand
from collections import defaultdict
import gzip

# Global Object
PRE_ASSEMBLED_LABEL = "pre_assembled"
TOOL_OUTPUT_MAP = {
    "spades": "contigs.fasta", 
    "megahit": "final.contigs.fa", 
    "flye": "assembly.fasta",
    "raven": "assembly.fasta",
    "medaka_flye": "consensus.fasta", 
    "medaka_raven": "consensus.fasta"
}

# Global Lists to hold sorted SRA samples (Used by get_data.smk constraints)
PAIRED_SRA = []
SINGLE_SRA = []

##########################################################
# --- 2. Function to obtain final outputs per module --- #
##########################################################
## list to generate all output
# workflow/scripts/functions.py

def get_final_outputs(include_multiqc=False):
  final_outputs = []

  # 1. Keep Download
  if MODULES["download_only"]:
    for sample, meta in SAMPLE_META.items():         
      if meta['mode'] == 'SRA':
        final_outputs.extend(meta['files'])
    return final_outputs

  if MODULES["keep_download"]: 
    for sample, meta in SAMPLE_META.items():         
        if meta['mode'] == 'SRA':
            final_outputs.extend(meta['files'])

  # 1b. Host Genome Filtering (BWA Depletion)
  if MODULES.get("host_filter", False):
    for sample, meta in SAMPLE_META.items():
        if meta['mode'] != 'CONTIGS':
            is_paired = (meta['mode'] == 'PAIRED') or (meta['mode'] == 'SRA' and len(meta['files']) == 2)
            hf_rule = "host_filter_paired" if is_paired else "host_filter_unpaired"
            final_outputs.append(f"{OUT_DIR}/{sample}/host_filter/{sample}_{hf_rule}_summary.txt")
            final_outputs.append(f"{OUT_DIR}/{sample}/stats/{sample}_{hf_rule}_stats.tsv")
            if is_paired:
                final_outputs.append(f"{OUT_DIR}/{sample}/host_filter/{sample}_host_filter_paired_1.fastq.gz")
                final_outputs.append(f"{OUT_DIR}/{sample}/host_filter/{sample}_host_filter_paired_2.fastq.gz")
            else:
                final_outputs.append(f"{OUT_DIR}/{sample}/host_filter/{sample}_host_filter_unpaired.fastq.gz")
  
  # 2. Assembly (Gera os arquivos físicos padronizados {sample}_{rule})
  if MODULES["assembly"]:  
    tools_list = MAPPER if isinstance(MAPPER, list) else [MAPPER]
    for tool_name in tools_list:
        if tool_name not in TOOL_OUTPUT_MAP:
             raise ValueError(f"Output filename not defined for tool: {tool_name}")
        
        if tool_name == "spades":
            kmer_list = config["spades"]["kmer"] 
            final_outputs.extend(expand(
                "{output_dir}/{sample}/spades/kmer_{kmer_val}/{sample}_spades_k{kmer_val}.fasta", 
                output_dir=OUT_DIR, sample=SAMPLE, kmer_val=kmer_list
            ))
        elif tool_name.startswith("medaka_"):
            assembler_base = tool_name.split("medaka_", 1)[1]
            final_outputs.extend(expand(
                "{output_dir}/{sample}/{tool}/{sample}_medaka_polish_{assembler}.fasta",
                output_dir=OUT_DIR, sample=SAMPLE, tool=tool_name, assembler=assembler_base
            ))
        else:
            final_outputs.extend(expand(
                "{output_dir}/{sample}/{tool}/{sample}_{tool}.fasta", 
                output_dir=OUT_DIR, sample=SAMPLE, tool=tool_name
            ))

  assembler_list = MAPPER if isinstance(MAPPER, list) else [MAPPER]
  
  # 3. Darkmatter
  if MODULES["darkmatter"]:
    for sample, meta in SAMPLE_META.items():
        # Lógica de expansão para esta amostra
        current_tools = []
        if meta['mode'] == 'CONTIGS':
             current_tools = [PRE_ASSEMBLED_LABEL]
        else:
             for tool in assembler_list:
                if tool == "spades":
                    # O CORREÇÃO CRÍTICA: Expande spades para os k-mers no Darkmatter também
                    kmer_list = config["spades"]["kmer"]
                    for k in kmer_list:
                        current_tools.append(f"spades_k{k}")
                else:
                    current_tools.append(tool)

        final_outputs.extend(expand("{out_dir}/{sample}/darkmatter_{tool}/{sample}_report_summarize_{tool}.html", 
                out_dir=OUT_DIR, sample=sample, tool=current_tools))
        final_outputs.extend(expand("{out_dir}/{sample}/darkmatter_{tool}/{sample}_dm_validate_{tool}_orfs.fasta", 
                out_dir=OUT_DIR, sample=sample, tool=current_tools))
        final_outputs.extend(expand("{out_dir}/{sample}/darkmatter_{tool}/{sample}_dm_validate_{tool}_contigs.fasta", 
                out_dir=OUT_DIR, sample=sample, tool=current_tools))
        if MODULES.get("rvdb", False):
            final_outputs.extend(expand("{out_dir}/{sample}/rvdb_{tool}/{sample}_rvdb_summarize_{tool}.tsv", 
                    out_dir=OUT_DIR, sample=sample, tool=current_tools))
            final_outputs.extend(expand("{out_dir}/{sample}/rvdb_{tool}/{sample}_rvdb_validate_structural_{tool}_orfs.fasta", 
                    out_dir=OUT_DIR, sample=sample, tool=current_tools))
            final_outputs.extend(expand("{out_dir}/{sample}/rvdb_{tool}/{sample}_rvdb_validate_structural_{tool}_contigs.fasta", 
                    out_dir=OUT_DIR, sample=sample, tool=current_tools))
            final_outputs.extend(expand("{out_dir}/{sample}/rvdb_{tool}/{sample}_rvdb_validate_pol_{tool}_orfs.fasta", 
                    out_dir=OUT_DIR, sample=sample, tool=current_tools))
            final_outputs.extend(expand("{out_dir}/{sample}/rvdb_{tool}/{sample}_rvdb_validate_pol_{tool}_contigs.fasta", 
                    out_dir=OUT_DIR, sample=sample, tool=current_tools))

  # 4. Reads (Diamond) + MEGAN LCA
  if MODULES["reads_diamond"]:
    final_outputs.extend(expand("{out_dir}/{sample}/diamond_reads/{sample}_append_lineage_reads.tsv", out_dir=OUT_DIR, sample=SAMPLE))
    final_outputs.extend(expand("{out_dir}/{sample}/krona_reads/{sample}_krona_diamond_reads.html", out_dir=OUT_DIR, sample=SAMPLE))
    final_outputs.extend(expand("{out_dir}/{sample}/megan_reads/{sample}_reads_daa2info.megan", out_dir=OUT_DIR, sample=SAMPLE))

  # 6. Reads (Kraken2)
  if MODULES["reads_kraken2"]:
    for sample in SAMPLE:
        meta = SAMPLE_META.get(sample)
        is_paired = (meta['mode'] == 'PAIRED') or (meta['mode'] == 'SRA' and len(meta['files']) == 2)
        label = "paired" if is_paired else "unpaired"
        final_outputs.append(f"{OUT_DIR}/{sample}/kraken2_reads/{sample}_kraken_biom_reads_{label}.txt")
        final_outputs.append(f"{OUT_DIR}/{sample}/krona_reads/{sample}_krona_reads_kraken_{label}.html")
        final_outputs.append(f"{OUT_DIR}/kraken2_all/all_samples_kraken_rarefaction_plot.pdf")

  # 6. Contigs (Kraken2) & 7. Contigs (Diamond)
  if MODULES["kraken2"] or MODULES["diamond"]:
      for sample, meta in SAMPLE_META.items():
          # Reutiliza a mesma lógica de expansão usada no Darkmatter
          current_tools = []
          if meta['mode'] == 'CONTIGS':
              current_tools = [PRE_ASSEMBLED_LABEL]
          else:
              for tool in assembler_list:
                  if tool == "spades":
                      kmer_list = config["spades"]["kmer"]
                      for k in kmer_list:
                          current_tools.append(f"spades_k{k}")
                  else:
                      current_tools.append(tool)

          # 6. Contigs (Kraken2)
          if MODULES["kraken2"]:
              final_outputs.extend(expand(
                  "{out_dir}/{sample}/kraken2_{tool}/{sample}_kraken_biom_contig_{tool}.txt",
                  out_dir=OUT_DIR, sample=sample, tool=current_tools
              ))
              final_outputs.extend(expand(
                  "{out_dir}/{sample}/krona_{tool}/{sample}_krona_kraken2_{tool}.html",
                  out_dir=OUT_DIR, sample=sample, tool=current_tools
              ))
              # Add Kraken2 statistics
              final_outputs.extend(expand(
                  "{out_dir}/{sample}/stats/{sample}_generate_kraken2_stats_{tool}.tsv",
                  out_dir=OUT_DIR, sample=sample, tool=current_tools
              ))

          # 7. Contigs (Diamond)
          if MODULES["diamond"]:
              final_outputs.extend(expand("{out_dir}/{sample}/diamond_{tool}/{sample}_append_lineage_{tool}.tsv", out_dir=OUT_DIR, sample=sample, tool=current_tools))
              final_outputs.extend(expand("{out_dir}/{sample}/krona_{tool}/{sample}_krona_diamond_{tool}.html", out_dir=OUT_DIR, sample=sample, tool=current_tools))
              # Add Diamond contig statistics
              final_outputs.extend(expand(
                  "{out_dir}/{sample}/stats/{sample}_generate_diamond_stats_{tool}.tsv",
                  out_dir=OUT_DIR, sample=sample, tool=current_tools
              ))

          # Assembly Statistics (for de novo assemblies)
          if MODULES["assembly"] and meta['mode'] != 'CONTIGS':
              final_outputs.extend(expand(
                  "{out_dir}/{sample}/stats/{sample}_generate_assembly_stats_{tool}.tsv",
                  out_dir=OUT_DIR, sample=sample, tool=current_tools
              ))

          # Dark Matter Statistics
          if MODULES["darkmatter"]:
              final_outputs.extend(expand(
                  "{out_dir}/{sample}/stats/{sample}_generate_darkmatter_stats_{tool}.tsv",
                  out_dir=OUT_DIR, sample=sample, tool=current_tools
              ))

  # Reads Diamond statistics
  if MODULES["reads_diamond"]:
    final_outputs.extend(expand("{out_dir}/{sample}/stats/{sample}_generate_diamond_reads_stats.tsv", out_dir=OUT_DIR, sample=SAMPLE))

  # Reads Kraken2 statistics  
  if MODULES["reads_kraken2"]:
    final_outputs.extend(expand("{out_dir}/{sample}/stats/{sample}_generate_kraken2_reads_stats.tsv", out_dir=OUT_DIR, sample=SAMPLE))

  # MultiQC Reports
  if include_multiqc:
      final_outputs.append(f"{OUT_DIR}/all_samples_multiqc_aggregate.html")
      final_outputs.extend(expand("{out_dir}/{sample}/multiqc/{sample}_multiqc_sample.html", out_dir=OUT_DIR, sample=SAMPLE))

  return final_outputs



## in process
# list of all expected inputs of multiqc
# In workflow/scripts/functions.py

def get_multiqc_inputs(wildcards=None, sample=None):
    if wildcards is not None:
        sample_id = wildcards.sample
    elif sample is not None:
        sample_id = sample
    else:
        raise ValueError("[ERROR] get_multiqc_inputs requires either 'wildcards' or 'sample' argument.")

    mqc_inputs = []
    meta = SAMPLE_META.get(sample_id)
    if not meta: return []

    is_paired = (meta['mode'] == 'PAIRED') or (meta['mode'] == 'SRA' and len(meta['files']) == 2)
    label_reads = "paired" if is_paired else "unpaired"
    fastp_rule = "fastp_paired" if is_paired else "fastp_unpaired"
    hf_rule = "host_filter_paired" if is_paired else "host_filter_unpaired"

    # Reads QC and Statistics (Skip for CONTIGS mode)
    if meta['mode'] != 'CONTIGS':
        mqc_inputs.append(os.path.join(OUT_DIR, sample_id, "fastp", f"{sample_id}_{fastp_rule}.json"))
        
        if MODULES.get("host_filter", False):
            mqc_inputs.append(os.path.join(OUT_DIR, sample_id, "stats", f"{sample_id}_{hf_rule}_stats.tsv"))
            mqc_inputs.append(os.path.join(OUT_DIR, sample_id, "host_filter", f"{sample_id}_{hf_rule}_summary.txt"))
        
        # Reads Diamond & Kraken2 Statistics
        if MODULES.get("reads_diamond", False):
            # Add Diamond read statistics summary
            mqc_inputs.append(os.path.join(OUT_DIR, sample_id, "stats", f"{sample_id}_generate_diamond_reads_stats.tsv"))
        
        if MODULES.get("reads_kraken2", False):
            mqc_inputs.append(os.path.join(OUT_DIR, sample_id, "kraken2_reads", f"{sample_id}_kraken2_reads_{label_reads}_report.txt"))
            # Add Kraken2 read statistics summary
            mqc_inputs.append(os.path.join(OUT_DIR, sample_id, "stats", f"{sample_id}_generate_kraken2_reads_stats.tsv"))

    # Contigs QC (Expandido por K-mer)
    assembler_list = MAPPER if isinstance(MAPPER, list) else [MAPPER]
    tools_expanded = []
    
    if meta['mode'] == 'CONTIGS':
        tools_expanded = [PRE_ASSEMBLED_LABEL]
    else:
        for tool in assembler_list:
            if tool == "spades":
                # Expande para k-mers
                kmer_list = config["spades"]["kmer"]
                for k in kmer_list:
                    tools_expanded.append(f"spades_k{k}")
            else:
                tools_expanded.append(tool)

    # Assembly Statistics (Number of contigs, N50, etc.)
    if MODULES.get("assembly", False) and meta['mode'] != 'CONTIGS':
        mqc_inputs.extend(expand("{out_dir}/{sample}/stats/{sample}_generate_assembly_stats_{tool}.tsv",
            out_dir=OUT_DIR, sample=sample_id, tool=tools_expanded))

    # Diamond Contig Statistics
    if MODULES.get("diamond", False):
        # Add Diamond contig statistics summary
        mqc_inputs.extend(expand("{out_dir}/{sample}/stats/{sample}_generate_diamond_stats_{tool}.tsv",
            out_dir=OUT_DIR, sample=sample_id, tool=tools_expanded))

    # Kraken2 Contig Statistics
    if MODULES.get("kraken2", False):
        mqc_inputs.extend(expand("{out_dir}/{sample}/kraken2_{tool}/{sample}_kraken2_contigs_{tool}_report.txt",
            out_dir=OUT_DIR, sample=sample_id, tool=tools_expanded))
        # Add Kraken2 contig statistics summary
        mqc_inputs.extend(expand("{out_dir}/{sample}/stats/{sample}_generate_kraken2_stats_{tool}.tsv",
            out_dir=OUT_DIR, sample=sample_id, tool=tools_expanded))

    # Dark Matter Statistics (Putative dark matter contigs)
    if MODULES.get("darkmatter", False):
        mqc_inputs.extend(expand("{out_dir}/{sample}/stats/{sample}_generate_darkmatter_stats_{tool}.tsv",
            out_dir=OUT_DIR, sample=sample_id, tool=tools_expanded))

    return mqc_inputs

#############################################################
# --- 3. Clean Up SRA Downloads after total processing --- #
############################################################
## ready to test
def cleanup_downloaded_fastqs(config):
    """
    Delete Fastq Data downloaded from SRA after workflown conclusion,
    unless 'keep_download' is TRUE.
    """
    options = config.get("modules", {})
    download_only = options.get("keep_download", False)

    # Verify if there is any SRA within metadata
    sra_samples_exist = any(meta['mode'] == 'SRA' for meta in SAMPLE_META.values())

    # 1: Not possible to clean (No SRA)
    if not sra_samples_exist:
        return

    # 2: User delanded to keep download
    if download_only:
        print(f"\n{blue}[INFO]{nc} Workflow (download_only) complete. *Keeping* SRA downloaded FASTQ files.")
        return

    # 3: Cleaning Space After Workflow completion
    print(f"\n{blue}[INFO]{nc} --- INITIATING POST-SUCCESS CLEANUP ---")
    print("Workflow complete. Removing temporary SRA FASTQ files...")
    
    files_removed = 0
    
    # Iteration over mmeta global in SAMPLE_META
    for sample, meta in SAMPLE_META.items():
        # Clean only SRA flag
        if meta['mode'] == 'SRA':
            # The list possesses paths
            for file_path in meta['files']:
                if os.path.exists(file_path):
                    try:
                        os.remove(file_path)
                        print(f"Removed: {file_path}")
                        files_removed += 1
                    except OSError as e:
                        print(f"{red}[WARNING]{nc} Could not remove {file_path}: {e}")
            
    print(f"{blue}[INFO]{nc} Cleanup complete. {files_removed} files removed.")
    print("---------------------------------------\n")

###########################################################
# --- 4. Obtain Novel Fasta with no Hits Iddentified --- #
##########################################################
## ready to test
def filter_nohits():
    from Bio import SeqIO
    # Access snakemake variables
    try:
        blast_file = snakemake.input.diamond
        fasta_file = snakemake.input.fasta
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
    print(f"Starting sample {sample_name}: filtering no-hit sequences", file=sys.stderr)
    print(f"  Log file: {log_file}", file=sys.stderr)
    print(f"  Input FASTA: {fasta_file}", file=sys.stderr)
    print(f"  DIAMOND results: {blast_file}", file=sys.stderr)
    print(f"  Output FASTA: {output_file}", file=sys.stderr)
 
    # Read BLAST hits
    ids_with_hits = set()
    try:
        with open(blast_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                if line.strip():
                    parts = line.split('\t')
                    if parts:  # Ensure line has content
                        ids_with_hits.add(parts[0])
    except FileNotFoundError:
        print(f"ERROR: DIAMOND hits file not found at {blast_file}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"ERROR reading DIAMOND file {blast_file}: {e}", file=sys.stderr)
        sys.exit(1)
    
    print(f"Found {len(ids_with_hits)} unique sequences with hits", file=sys.stderr)
    
    # Process and write sequences in one pass
    count_no_hits = 0
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
                if record.id not in ids_with_hits:
                    count_no_hits += 1
                    SeqIO.write(record, output_handle, "fasta")
                    
    except FileNotFoundError:
        print(f"ERROR: FASTA file not found at {fasta_file}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"ERROR processing FASTA file {fasta_file}: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Summary statistics
    print(f"Processed {count_total} total sequences", file=sys.stderr)
    print(f"Saved {count_no_hits} sequences with no hits ({count_no_hits/max(1,count_total)*100:.1f}%)", file=sys.stderr)
    print(f"Output written to: {output_file}", file=sys.stderr)
    print("Filtering complete.", file=sys.stderr)

####################################################
# --- 5. Obtain Metadata to Library Processing --- #
####################################################
## ready to test
def identify_data_type(sample_list, data_dir):
    """
    Verify if each sample within the file is a existing contig file, paired fastq or unpaired fastq.
    If none is found, signals for download (SRA).
    
    Returns:
        sample_meta (dict): Dictionary with metadada for sample.
                            Structure: {
                                'sample_A': {'mode': 'PAIRED', 'files': [path1, path2]},
                                'sample_B': {'mode': 'SRA', 'files': [path1_future, path2_future]}
                            }
    """
    sample_meta = {}
    
    for sample in sample_list:
        # Guard: if the sample name is a subdirectory inside data_dir, skip it.
        # This prevents folders (e.g. per-sample BCL2FASTQ dirs) from being
        # mistakenly queued as SRA downloads.
        if os.path.isdir(os.path.join(data_dir, sample)):
            print(f"[WARNING] '{sample}' is a directory inside '{data_dir}' — skipping (not treated as SRA).")
            continue

        # 1. Definição de caminhos esperados (Prioridade de checagem)
        
        # A. Contigs (.fasta or .fa or .fas)
        path_fasta = os.path.join(data_dir, f"{sample}.fasta")
        path_fa = os.path.join(data_dir, f"{sample}.fa")
        path_fas = os.path.join(data_dir, f"{sample}.fas")
        
        # B. Paired FastQ (_R1/_R2 ou _1/_2)
        p_r1 = os.path.join(data_dir, f"{sample}_R1.fastq.gz")
        p_r2 = os.path.join(data_dir, f"{sample}_R2.fastq.gz")
        p_1  = os.path.join(data_dir, f"{sample}_1.fastq.gz")
        p_2  = os.path.join(data_dir, f"{sample}_2.fastq.gz")
        p_01  = os.path.join(data_dir, f"{sample}_R1_001.fastq.gz")
        p_02  = os.path.join(data_dir, f"{sample}_R2_001.fastq.gz")
        
        # C. Unpaired FastQ
        p_unpR1 = os.path.join(data_dir, f"{sample}_R1.fastq.gz")
        p_unp1  = os.path.join(data_dir, f"{sample}_1.fastq.gz")
        p_unp  = os.path.join(data_dir, f"{sample}.fastq.gz")
        p_unp01  = os.path.join(data_dir, f"{sample}_R1_001.fastq.gz")

        # 1. Contig?
        if os.path.exists(path_fasta):
            sample_meta[sample] = {'mode': 'CONTIGS', 'files': [path_fasta]}
        elif os.path.exists(path_fa):
            sample_meta[sample] = {'mode': 'CONTIGS', 'files': [path_fa]}
        elif os.path.exists(path_fas):
            sample_meta[sample] = {'mode': 'CONTIGS', 'files': [path_fas]}
            
        # 2. Paired?
        elif os.path.exists(p_r1) and os.path.exists(p_r2):
            sample_meta[sample] = {'mode': 'PAIRED', 'files': [p_r1, p_r2]}
        elif os.path.exists(p_1) and os.path.exists(p_2):
            sample_meta[sample] = {'mode': 'PAIRED', 'files': [p_1, p_2]}
        elif os.path.exists(p_01) and os.path.exists(p_02):
            sample_meta[sample] = {'mode': 'PAIRED', 'files': [p_01, p_02]}

        # 3. Unpaired?
        elif os.path.exists(p_unpR1):
            sample_meta[sample] = {'mode': 'UNPAIRED', 'files': [p_unpR1]}
        elif os.path.exists(p_unp1):
            sample_meta[sample] = {'mode': 'UNPAIRED', 'files': [p_unp1]}
        elif os.path.exists(p_unp):
            sample_meta[sample] = {'mode': 'UNPAIRED', 'files': [p_unp]}
        elif os.path.exists(p_unp01):
            sample_meta[sample] = {'mode': 'UNPAIRED', 'files': [p_unp01]}
            
        # 4. None? Download (SRA)
        else:
            # SRA downloads always land in the parent data_dir, not a subdirectory
            future_r1 = os.path.join(data_dir, f"{sample}_1.fastq.gz")
            future_r2 = os.path.join(data_dir, f"{sample}_2.fastq.gz")
            
            # API CHECK: Is this Single or Paired?
            layout = get_sra_layout(sample)
            
            if layout == 'SINGLE':
                #print(f"[INFO] {sample} identified as SRA SINGLE END.")
                # We use SRA mode but only list 1 file. 
                # This will trigger fastp_unpaired.
                sample_meta[sample] = {'mode': 'SRA', 'files': [future_r1]}
                SINGLE_SRA.append(sample)
            else:
                #print(f"[INFO] {sample} identified as SRA PAIRED END.")
                sample_meta[sample] = {'mode': 'SRA', 'files': [future_r1, future_r2]}
                PAIRED_SRA.append(sample)

    return sample_meta

####################################################################
# --- 6. Obtain path to Library Reads availablle or downloaded --- #
####################################################################
## ready to test
def get_input_r1(wildcards):
    """Returns R1 local or future"""
    meta = SAMPLE_META[wildcards.sample]
    if meta['mode'] in ['PAIRED', 'SRA']:
        return meta['files'][0]
    return []

def get_input_r2(wildcards):
    """Returns R2 local or future"""
    meta = SAMPLE_META[wildcards.sample]
    if meta['mode'] in ['PAIRED', 'SRA']:
        return meta['files'][1]
    return []

def get_input_unp(wildcards):
    """Returns R1 local or future for unpaired data"""
    meta = SAMPLE_META[wildcards.sample]
    if meta['mode'] in ['UNPAIRED', 'SRA']:
        return meta['files'][0]
    return []

########################################################################
# --- 7. Function to Get Contigs paths according to assembly tool --- #
########################################################################
def get_contigs_path(wildcards):
    """
    Determines input file based on the wildcard label.
    Handles virtual tools like 'spades_k21'.
    """
    tool_name = wildcards.tool
    sample = wildcards.sample
    meta = SAMPLE_META.get(sample)
    
    if not meta:
        raise ValueError(f"Metadata missing for {sample}")

    # 1. Caso Pre-assembled
    if tool_name == PRE_ASSEMBLED_LABEL:
        if meta['mode'] != 'CONTIGS':
            raise ValueError(f"Sample '{sample}' is marked as {meta['mode']}, but '{PRE_ASSEMBLED_LABEL}' was requested.")
        return meta['files'][0]

    # 2. Caso Virtual Tool: SPADES com K-mer (ex: spades_k33)
    if tool_name.startswith("spades_k"):
        try:
            kmer_val = tool_name.split("_k")[1]
            return os.path.join(OUT_DIR, sample, "spades", f"kmer_{kmer_val}", f"{sample}_spades_k{kmer_val}.fasta")
        except IndexError:
             raise ValueError(f"Could not parse kmer value from tool name: {tool_name}")

    # 3. Caso Medaka Polish (medaka_flye / medaka_raven)
    elif tool_name.startswith("medaka_"):
        assembler_base = tool_name.split("medaka_", 1)[1]
        return os.path.join(OUT_DIR, sample, tool_name, f"{sample}_medaka_polish_{assembler_base}.fasta")

    # 4. Caso Genérico (Megahit, Flye, Raven)
    else:
        if tool_name not in TOOL_OUTPUT_MAP:
             raise ValueError(f"Tool '{tool_name}' not recognized in TOOL_OUTPUT_MAP.")
        return os.path.join(OUT_DIR, sample, tool_name, f"{sample}_{tool_name}.fasta")

#####################################################
# --- 8. Helper Functions of Get De novo input --- #
####################################################
def _get_reads_stage_dir():
    """Return 'host_filter' when host_filter module is active, otherwise 'fastp'."""
    if isinstance(MODULES, dict) and MODULES.get("host_filter", False):
        return "host_filter"
    return "fastp"

def get_denovo_r1(wildcards):
    """Get R1 Only for paired-end(Paired/SRA)."""
    meta = SAMPLE_META.get(wildcards.sample)
    stage_dir = _get_reads_stage_dir()
    rule_prefix = f"{stage_dir}_paired"
    if meta and meta['mode'] in ['PAIRED', 'SRA'] and len(meta['files']) > 1:
        return os.path.join(OUT_DIR, wildcards.sample, stage_dir, f"{wildcards.sample}_{rule_prefix}_1.fastq.gz")
    return []

def get_denovo_r2(wildcards):
    """Get R2 Only for paired-end(Paired/SRA)."""
    meta = SAMPLE_META.get(wildcards.sample)
    stage_dir = _get_reads_stage_dir()
    rule_prefix = f"{stage_dir}_paired"
    if meta and meta['mode'] in ['PAIRED', 'SRA'] and len(meta['files']) > 1:
        return os.path.join(OUT_DIR, wildcards.sample, stage_dir, f"{wildcards.sample}_{rule_prefix}_2.fastq.gz")   
    return []

def get_denovo_unpaired(wildcards):
    """
    Get files for the flag -s (single/unpaired).
    - PAIRED: Return 'orphans' paired reads.
    - UNPAIRED or SRA-SINGLE: Return 'unpaired' main reads.
    """
    meta = SAMPLE_META.get(wildcards.sample)
    stage_dir = _get_reads_stage_dir()

    if not meta:
        return []
    
    # 1. True Paired (SRA or Local) -> Return Orphans
    if meta['mode'] in ['PAIRED', 'SRA'] and len(meta['files']) > 1:
        return os.path.join(OUT_DIR, wildcards.sample, stage_dir, f"{wildcards.sample}_{stage_dir}_paired_orphans.fastq.gz")
    
    # 2. True Single (SRA or Local) -> Return Unpaired Main Reads
    is_sra_single = (meta['mode'] == 'SRA' and len(meta['files']) == 1)
    
    if meta['mode'] == 'UNPAIRED' or is_sra_single:
        return os.path.join(OUT_DIR, wildcards.sample, stage_dir, f"{wildcards.sample}_{stage_dir}_unpaired.fastq.gz")
        
    return []

def get_ONP_input(wildcards):
    """
    Specific input function for Flye.
    1. Checks if the sample is UNPAIRED (Nanopore/PacBio).
    2. If PAIRED, returns empty list (skips Flye for Illumina samples).
    """
    meta = SAMPLE_META.get(wildcards.sample)
    stage_dir = _get_reads_stage_dir()

    if not meta:
        return []

    is_sra_single = (meta['mode'] == 'SRA' and len(meta['files']) == 1)

    if meta['mode'] == 'UNPAIRED' or is_sra_single:
        return os.path.join(OUT_DIR, wildcards.sample, stage_dir, f"{wildcards.sample}_{stage_dir}_unpaired.fastq.gz")
    
    return []

#####################################################
# --- 9. Helper Functions of Get De novo params --- #
####################################################
def get_spades_params(wildcards, input):
    """Build args dynamically."""
    cmd = ""
    
    # 1. R1 & R2 (For Paired End)
    # Check if attributes exist and are not empty
    if hasattr(input, 'r1') and hasattr(input, 'r2') and input.r1 and input.r2:
        cmd += f"-1 {input.r1} -2 {input.r2} "
    
    # 2. Extra/Unpaired (For Single End or Orphans)
    # FIX: Use 'input.extra' because that is the name defined in the rule
    if hasattr(input, 'extra') and input.extra:
        cmd += f"-s {input.extra}"
        
    return cmd

def get_megahit_params(wildcards, input):
    """Build args dynamically."""
    cmd = ""
    
    # 1. R1 & R2
    if hasattr(input, 'r1') and hasattr(input, 'r2') and input.r1 and input.r2:
        cmd += f"-1 {input.r1} -2 {input.r2} "
    
    # 2. Extra/Unpaired
    # FIX: Use 'input.extra' here as well
    if hasattr(input, 'extra') and input.extra:
        cmd += f"-r {input.extra}"
        
    return cmd


################################
# --- 10. Check SRA Layout --- #
################################
def get_sra_layout(accession):
    """
    Queries ENA API to check if an SRA accession is PAIRED or SINGLE.
    Defaults to PAIRED on failure to be safe, or SINGLE if specified.
    """
    url = f"https://www.ebi.ac.uk/ena/portal/api/filereport?accession={accession}&result=read_run&fields=library_layout&format=tsv"
    try:
        with urllib.request.urlopen(url) as response:
            data = response.read().decode('utf-8')
            # Response format: header\naccession\tLAYOUT
            if "SINGLE" in data:
                return "SINGLE"
            return "PAIRED"
    except Exception as e:
        print(f"[WARNING] Could not check layout for {accession}: {e}. Defaulting to PAIRED.")
        return "PAIRED"


##########################################
# --- 11. Diamond Database Detection --- #
##########################################

def get_diamond_db_input(wildcards):
    """
    Retorna a lista de arquivos do banco de dados Diamond.
    Lida tanto com arquivo único (.dmnd) quanto com banco BLAST particionado (nr.*).
    Caso nenhum banco externo/local exista em resources/diamond, retorna o caminho do
    banco base de teste (ecoli_db.dmnd) para acionar a regra download_default_diamond_db.
    """
    _script_dir = os.path.dirname(os.path.abspath(__file__))
    _project_root = os.path.normpath(os.path.join(_script_dir, "..", ".."))
    db_dir = os.path.join(_project_root, "resources", "diamond")

    # 1. Tenta achar qualquer arquivo .dmnd existente e válido (ignora symlinks quebrados)
    dmnd_file = [f for f in glob.glob(os.path.join(db_dir, "*.dmnd")) if os.path.exists(f)]
    if dmnd_file:
        return dmnd_file

    # 2. Procura por arquivos de índice do BLAST particionado (ex: .acc, .phr, .psq, .pal)
    blast_idx = [
        f for f in (
            glob.glob(os.path.join(db_dir, "*.acc"))
            + glob.glob(os.path.join(db_dir, "*.phr"))
            + glob.glob(os.path.join(db_dir, "*.psq"))
            + glob.glob(os.path.join(db_dir, "*.pal"))
        )
        if os.path.exists(f)
    ]
    if blast_idx:
        return [
            f for f in glob.glob(os.path.join(db_dir, "*"))
            if not os.path.basename(f).startswith(".") and os.path.exists(f)
        ]

    # 3. Fallback automático de 1ª execução: solicita o banco base E. coli K-12 (ecoli_db.dmnd)
    return [os.path.join(db_dir, "ecoli_db.dmnd")]

def get_diamond_db_name(wildcards):
    """
    Retorna o NOME BASE para o comando do Diamond (-d).
    Caso nenhum banco exista ainda (antes da regra download_default_diamond_db rodar),
    retorna 'ecoli_db.dmnd'.
    """
    _script_dir = os.path.dirname(os.path.abspath(__file__))
    _project_root = os.path.normpath(os.path.join(_script_dir, "..", ".."))
    db_dir = os.path.join(_project_root, "resources", "diamond")

    # Caso 1: Arquivo .dmnd único e válido
    dmnd_files = [f for f in glob.glob(os.path.join(db_dir, "*.dmnd")) if os.path.exists(f)]
    if dmnd_files:
        return os.path.basename(dmnd_files[0])

    # Caso 2: Banco BLAST particionado
    acc_files = [f for f in glob.glob(os.path.join(db_dir, "*.acc")) if os.path.exists(f)]
    if acc_files:
        filename = os.path.basename(acc_files[0])
        prefix = re.sub(r'\.\d+\.acc$', '', filename)
        if prefix == filename:
            prefix = os.path.splitext(filename)[0]
        return prefix

    # Caso 3: Banco base padrão (será criado pela regra download_default_diamond_db)
    return "ecoli_db.dmnd"

########################################################
# --- 12. Check Presence of header in Diamond File --- #
########################################################

def get_header_check(hit_file_path):
    """
    Checks if the DIAMOND hit file has a header (starting with 'qseqid').
    Returns 1 if header is present, 0 otherwise.
    """
    try:
        # Tenta abrir o arquivo. Não precisamos de gzip aqui, 
        # pois o DIAMOND output não deve ser comprimido.
        with open(hit_file_path, 'r') as f:
            first_line = f.readline()
            if first_line.startswith('qseqid'):
                return 1
            return 0
    except Exception as e:
        print(f"[ERROR] Could not read hit file header check: {e}", file=sys.stderr)
        return 0

####################################################
# --- 13. Get Taxonomy reports of all samples --- #
###################################################
def get_all_kraken_reports(wildcards):
    paths = []
    for s in SAMPLE:
        meta = SAMPLE_META.get(s)
        if meta:
             is_paired = (meta['mode'] == 'PAIRED') or (meta['mode'] == 'SRA' and len(meta['files']) == 2)
             label = "paired" if is_paired else "unpaired"
             paths.append(os.path.join(OUT_DIR, s, "kraken2_reads", f"{s}_kraken2_reads_{label}_report.txt"))
    return paths