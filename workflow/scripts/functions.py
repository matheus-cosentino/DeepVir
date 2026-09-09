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
#                         version: 1.0  |  09.2026                              #
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
  
  # 2. Assembly (Mantemos igual, gera os arquivos físicos)
  if MODULES["assembly"]:  
    tools_list = MAPPER if isinstance(MAPPER, list) else [MAPPER]
    for tool_name in tools_list:
        file_name = TOOL_OUTPUT_MAP.get(tool_name)
        if not file_name:
             raise ValueError(f"Output filename not defined for tool: {tool_name}")
        
        if tool_name == "spades":
            kmer_list = config["spades"]["kmer"] 
            final_outputs.extend(expand(
                "{output_dir}/{sample}/spades/kmer_{kmer_val}/{filename}", 
                output_dir=OUT_DIR, sample=SAMPLE, filename=file_name, kmer_val=kmer_list
            ))
        else:
            final_outputs.extend(expand(
                "{output_dir}/{sample}/{tool}/{filename}", 
                output_dir=OUT_DIR, sample=SAMPLE, tool=tool_name, filename=file_name
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

        final_outputs.extend(expand("{out_dir}/{sample}/darkmatter_{tool}/{sample}_Report_Diversity.html", 
                out_dir=OUT_DIR, sample=sample, tool=current_tools))
        final_outputs.extend(expand("{out_dir}/{sample}/darkmatter_{tool}/{sample}_RdRp_Orfs.fasta", 
                out_dir=OUT_DIR, sample=sample, tool=current_tools))
        if MODULES.get("rvdb", False):
            final_outputs.extend(expand("{out_dir}/{sample}/rvdb_{tool}/{sample}_structural_contigs.fasta", 
                    out_dir=OUT_DIR, sample=sample, tool=current_tools))
            final_outputs.extend(expand("{out_dir}/{sample}/rvdb_{tool}/{sample}_pol_contigs.fasta", 
                    out_dir=OUT_DIR, sample=sample, tool=current_tools))

  # 4. Reads (Diamond) + MEGAN LCA
  if MODULES["reads_diamond"]:
    final_outputs.extend(expand("{out_dir}/{sample}/diamond_reads/{sample}_reads_hits_with_lineage.tsv", out_dir=OUT_DIR, sample=SAMPLE))
    final_outputs.extend(expand("{out_dir}/{sample}/krona_reads/{sample}_reads_diamond_krona.html", out_dir=OUT_DIR, sample=SAMPLE))
    final_outputs.extend(expand("{out_dir}/{sample}/megan_reads/{sample}_reads_summary.megan", out_dir=OUT_DIR, sample=SAMPLE))

  # 6. Reads (Kraken2)
  if MODULES["reads_kraken2"]:
    for sample in SAMPLE:
        meta = SAMPLE_META.get(sample)
        is_paired = (meta['mode'] == 'PAIRED') or (meta['mode'] == 'SRA' and len(meta['files']) == 2)
        label = "paired" if is_paired else "unpaired"
        final_outputs.append(f"{OUT_DIR}/{sample}/kraken2_reads/{sample}_{label}_reads_biom.txt")
        final_outputs.append(f"{OUT_DIR}/{sample}/krona_reads/{sample}_{label}_kraken2_krona.html")
        final_outputs.append(f"{OUT_DIR}/kraken2_all/Rarefaction_Curve.pdf")

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
                  "{out_dir}/{sample}/kraken2_{tool}/{sample}_{tool}_contig_biom.txt",
                  out_dir=OUT_DIR, sample=sample, tool=current_tools
              ))
              final_outputs.extend(expand(
                  "{out_dir}/{sample}/krona_{tool}/{sample}_{tool}_kraken2_krona.html",
                  out_dir=OUT_DIR, sample=sample, tool=current_tools
              ))
              # Add Kraken2 statistics
              final_outputs.extend(expand(
                  "{out_dir}/{sample}/stats/{sample}_{tool}_kraken2_stats.tsv",
                  out_dir=OUT_DIR, sample=sample, tool=current_tools
              ))

          # 7. Contigs (Diamond)
          if MODULES["diamond"]:
              final_outputs.extend(expand("{out_dir}/{sample}/diamond_{tool}/{sample}_{tool}_hits_with_lineage.tsv", out_dir=OUT_DIR, sample=sample, tool=current_tools))
              final_outputs.extend(expand("{out_dir}/{sample}/krona_{tool}/{sample}_{tool}_diamond_krona.html",out_dir=OUT_DIR, sample=sample, tool=current_tools))
              # Add Diamond contig statistics
              final_outputs.extend(expand(
                  "{out_dir}/{sample}/stats/{sample}_{tool}_diamond_stats.tsv",
                  out_dir=OUT_DIR, sample=sample, tool=current_tools
              ))

          # Assembly Statistics (for de novo assemblies)
          if MODULES["assembly"] and meta['mode'] != 'CONTIGS':
              final_outputs.extend(expand(
                  "{out_dir}/{sample}/stats/{sample}_{tool}_assembly_stats.tsv",
                  out_dir=OUT_DIR, sample=sample, tool=current_tools
              ))

          # Dark Matter Statistics
          if MODULES["darkmatter"]:
              final_outputs.extend(expand(
                  "{out_dir}/{sample}/stats/{sample}_{tool}_darkmatter_stats.tsv",
                  out_dir=OUT_DIR, sample=sample, tool=current_tools
              ))

  # Reads Diamond statistics
  if MODULES["reads_diamond"]:
    final_outputs.extend(expand("{out_dir}/{sample}/stats/{sample}_reads_diamond_stats.tsv", out_dir=OUT_DIR, sample=SAMPLE))

  # Reads Kraken2 statistics  
  if MODULES["reads_kraken2"]:
    final_outputs.extend(expand("{out_dir}/{sample}/stats/{sample}_reads_kraken2_stats.tsv", out_dir=OUT_DIR, sample=SAMPLE))

  # MultiQC Reports
  if include_multiqc:
      final_outputs.append(f"{OUT_DIR}/multiqc_report.html")
      final_outputs.extend(expand("{out_dir}/{sample}/multiqc/multiqc_report.html", out_dir=OUT_DIR, sample=SAMPLE))

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
    fastp_suffix = "paired" if is_paired else "unp"

    # Reads QC and Statistics (Skip for CONTIGS mode)
    if meta['mode'] != 'CONTIGS':
        mqc_inputs.append(os.path.join(OUT_DIR, sample_id, "fastp", f"{sample_id}_{fastp_suffix}.json"))
        
        # Reads Diamond & Kraken2 Statistics
        if MODULES.get("reads_diamond", False):
            # Add Diamond read statistics summary
            mqc_inputs.append(os.path.join(OUT_DIR, sample_id, "stats", "{}_reads_diamond_stats.tsv".format(sample_id)))
        
        if MODULES.get("reads_kraken2", False):
            mqc_inputs.append(os.path.join(OUT_DIR, sample_id, "kraken2_reads", f"{sample_id}_{label_reads}_reads_report.txt"))
            # Add Kraken2 read statistics summary
            mqc_inputs.append(os.path.join(OUT_DIR, sample_id, "stats", "{}_reads_kraken2_stats.tsv".format(sample_id)))

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
    if meta['mode'] != 'CONTIGS':
        mqc_inputs.extend(expand("{out_dir}/{sample}/stats/{sample}_{tool}_assembly_stats.tsv",
            out_dir=OUT_DIR, sample=sample_id, tool=tools_expanded))

    # Diamond Contig Statistics
    if MODULES.get("diamond", False):
        # Add Diamond contig statistics summary
        mqc_inputs.extend(expand("{out_dir}/{sample}/stats/{sample}_{tool}_diamond_stats.tsv",
            out_dir=OUT_DIR, sample=sample_id, tool=tools_expanded))

    # Kraken2 Contig Statistics
    if MODULES.get("kraken2", False):
        mqc_inputs.extend(expand("{out_dir}/{sample}/kraken2_{tool}/{sample}_{tool}_contig_report.txt",
            out_dir=OUT_DIR, sample=sample_id, tool=tools_expanded))
        # Add Kraken2 contig statistics summary
        mqc_inputs.extend(expand("{out_dir}/{sample}/stats/{sample}_{tool}_kraken2_stats.tsv",
            out_dir=OUT_DIR, sample=sample_id, tool=tools_expanded))

    # Dark Matter Statistics (Putative dark matter contigs)
    if MODULES.get("darkmatter", False):
        mqc_inputs.extend(expand("{out_dir}/{sample}/stats/{sample}_{tool}_darkmatter_stats.tsv",
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
    
    #print(f"[INFO] Checking input availability for {len(sample_list)} samples...")

    for sample in sample_list:
        # If the sample name matches a subdirectory inside data_dir, look for
        # FASTQ/FASTA files *inside* that subdirectory before trying data_dir.
        # This handles BCL2FASTQ-style output where each sample has its own folder.
        sample_dir = os.path.join(data_dir, sample)
        if os.path.isdir(sample_dir):
            print(f"[INFO] '{sample}' is a subdirectory inside '{data_dir}' — searching for reads inside it.")
            data_dir_for_sample = sample_dir
        else:
            data_dir_for_sample = data_dir

        # 1. Definição de caminhos esperados (Prioridade de checagem)
        
        # A. Contigs (.fasta or .fa or .fas)
        path_fasta = os.path.join(data_dir_for_sample, f"{sample}.fasta")
        path_fa = os.path.join(data_dir_for_sample, f"{sample}.fa")
        path_fas = os.path.join(data_dir_for_sample, f"{sample}.fas")
        
        # B. Paired FastQ (_R1/_R2 ou _1/_2)
        p_r1 = os.path.join(data_dir_for_sample, f"{sample}_R1.fastq.gz")
        p_r2 = os.path.join(data_dir_for_sample, f"{sample}_R2.fastq.gz")
        p_1  = os.path.join(data_dir_for_sample, f"{sample}_1.fastq.gz")
        p_2  = os.path.join(data_dir_for_sample, f"{sample}_2.fastq.gz")
        p_01  = os.path.join(data_dir_for_sample, f"{sample}_R1_001.fastq.gz")
        p_02  = os.path.join(data_dir_for_sample, f"{sample}_R2_001.fastq.gz")
        
        # C. Unpaired FastQ
        p_unpR1 = os.path.join(data_dir_for_sample, f"{sample}_R1.fastq.gz")
        p_unp1  = os.path.join(data_dir_for_sample, f"{sample}_1.fastq.gz")
        p_unp  = os.path.join(data_dir_for_sample, f"{sample}.fastq.gz")
        p_unp01  = os.path.join(data_dir_for_sample, f"{sample}_R1_001.fastq.gz")

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
            # Pega o numero depois do "k" (spades_k33 -> 33)
            kmer_val = tool_name.split("_k")[1]
            filename = TOOL_OUTPUT_MAP.get("spades")
            
            return os.path.join(OUT_DIR, sample, "spades", f"kmer_{kmer_val}", filename)
        except IndexError:
             raise ValueError(f"Could not parse kmer value from tool name: {tool_name}")

    # 3. Caso Genérico (Megahit, Flye, Raven, ou o próprio 'spades' se chamado incorretamente)
    else:
        # Se por algum erro a regra chamar 'spades' puro, isso vai falhar lá na frente
        # mas aqui retornamos o caminho padrão.
        filename = TOOL_OUTPUT_MAP.get(tool_name)
        if not filename:
             # Tenta limpar sufixos caso haja algo estranho, ou lança erro
             raise ValueError(f"Tool '{tool_name}' not recognized in TOOL_OUTPUT_MAP.")
            
        return os.path.join(OUT_DIR, sample, tool_name, filename)

#####################################################
# --- 8. Helper Functions of Get De novo input --- #
####################################################
def get_denovo_r1(wildcards):
    """Get R1 Only for paired-end(Paired/SRA)."""
    meta = SAMPLE_META.get(wildcards.sample)
    # FIX: Check if len(files) > 1 to ensure it is actually paired
    if meta and meta['mode'] in ['PAIRED', 'SRA'] and len(meta['files']) > 1:
        return os.path.join(OUT_DIR, wildcards.sample, "fastp", f"{wildcards.sample}_1.fastq.gz")
    return []

def get_denovo_r2(wildcards):
    """Get R2 Only for paired-end(Paired/SRA)."""
    meta = SAMPLE_META.get(wildcards.sample)
    # FIX: Check if len(files) > 1 to ensure it is actually paired
    if meta and meta['mode'] in ['PAIRED', 'SRA'] and len(meta['files']) > 1:
        return os.path.join(OUT_DIR, wildcards.sample, "fastp", f"{wildcards.sample}_2.fastq.gz")   
    return []

def get_denovo_unpaired(wildcards):
    """
    Get files for the flag -s (single/unpaired).
    - PAIRED: Return 'orphans' paired reads.
    - UNPAIRED or SRA-SINGLE: Return 'unp' main reads.
    """
    meta = SAMPLE_META.get(wildcards.sample)

    if not meta:
        return []
    
    # 1. True Paired (SRA or Local) -> Return Orphans
    if meta['mode'] in ['PAIRED', 'SRA'] and len(meta['files']) > 1:
        return os.path.join(OUT_DIR, wildcards.sample, "fastp", f"{wildcards.sample}_orphans.fastq.gz")
    
    # 2. True Single (SRA or Local) -> Return Unpaired Main Reads
    # We check if mode is UNPAIRED OR if mode is SRA with only 1 file
    is_sra_single = (meta['mode'] == 'SRA' and len(meta['files']) == 1)
    
    if meta['mode'] == 'UNPAIRED' or is_sra_single:
        # This forces Snakemake to use fastp_unpaired because it needs the '_unp' file
        return os.path.join(OUT_DIR, wildcards.sample, "fastp", f"{wildcards.sample}_unp.fastq.gz")
        
    return []

def get_ONP_input(wildcards):
    """
    Specific input function for Flye.
    1. Checks if the sample is UNPAIRED (Nanopore/PacBio).
    2. If PAIRED, returns empty list (skips Flye for Illumina samples).
    """
    meta = SAMPLE_META.get(wildcards.sample)

    if not meta:
        return []

    is_sra_single = (meta['mode'] == 'SRA' and len(meta['files']) == 1)

    if meta['mode'] == 'UNPAIRED' or is_sra_single:
        # For single end, we use the fastp 'unp' output
        return os.path.join(OUT_DIR, wildcards.sample, "fastp", f"{wildcards.sample}_unp.fastq.gz")
    
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
    """
    # Usa caminho absoluto derivado da localização deste script para evitar
    # problemas com --directory do Snakemake (que muda o CWD do processo).
    _script_dir = os.path.dirname(os.path.abspath(__file__))
    _project_root = os.path.normpath(os.path.join(_script_dir, "..", ".."))
    db_dir = os.path.join(_project_root, "resources", "diamond")
    
    # Tenta achar um .dmnd clássico
    dmnd_file = glob.glob(os.path.join(db_dir, "*.dmnd"))
    if dmnd_file:
        return dmnd_file
        
    # Se não, procura por arquivos de índice do BLAST (ex: .acc, .phr, .psq)
    # Pegamos tudo que estiver na pasta para garantir que o Snakemake monitore
    blast_files = glob.glob(os.path.join(db_dir, "*"))
    
    if not blast_files:
        raise ValueError(f"Nenhum banco de dados Diamond encontrado em {db_dir}")
        
    return blast_files

def get_diamond_db_name(wildcards):
    """
    Retorna o NOME BASE para o comando do Diamond (-d).
    """
    # Usa caminho absoluto derivado da localização deste script para evitar
    # problemas com --directory do Snakemake (que muda o CWD do processo).
    _script_dir = os.path.dirname(os.path.abspath(__file__))
    _project_root = os.path.normpath(os.path.join(_script_dir, "..", ".."))
    db_dir = os.path.join(_project_root, "resources", "diamond")
    
    # Caso 1: Arquivo .dmnd único
    dmnd_files = glob.glob(os.path.join(db_dir, "*.dmnd"))
    if dmnd_files:
        return os.path.basename(dmnd_files[0]) # Ex: database.dmnd
        
    # Caso 2: Banco BLAST particionado
    # Procura arquivos comuns do índice para deduzir o prefixo (ex: nr.00.acc -> nr)
    acc_files = glob.glob(os.path.join(db_dir, "*.acc"))
    if acc_files:
        filename = os.path.basename(acc_files[0])
        # Remove sufixos .número.acc para pegar o prefixo limpo (ex: nr.78.acc -> nr)
        # Ajuste a regex conforme o padrão exato dos seus arquivos
        prefix = re.sub(r'\.\d+\.acc$', '', filename) 
        # Fallback simples se não tiver número
        if prefix == filename:
            prefix = os.path.splitext(filename)[0]
        return prefix
        
    # Fallback genérico: pega o nome do primeiro arquivo sem extensão
    all_files = glob.glob(os.path.join(db_dir, "*"))
    if not all_files:
        print(f"[WARNING] Diamond DB directory '{db_dir}' is empty. Returning placeholder for dry-run.", file=sys.stderr)
        return "missing_db"
    first_file = os.path.basename(all_files[0])
    return first_file.split('.')[0]

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
             paths.append(os.path.join(OUT_DIR, s, "kraken2_reads", f"{s}_{label}_reads_report.txt"))
    return paths