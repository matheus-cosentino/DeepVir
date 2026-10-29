# Expected Outputs & Execution Logs

DeepVir separates scientific deliverables from execution logs across two dedicated directories:
1. **`results/` (`--output`)**: All scientific outputs, assembled contigs, taxonomic tables, viral dark matter candidates, and interactive HTML dashboards.
2. **`logs/`**: All per-sample Snakemake rule logs, workflow execution logs, GUI logs, and run configuration snapshots.

---

## 1. Scientific Outputs (`results/`)

Once executed via `python DeepVir.py`, your results are organized inside `--output` (e.g., `results/`):

```text
results/
├── kraken2_all/                   # Merged Kraken2 outputs and Alpha-Diversity
│   ├── all_samples.biom           # BIOM format table encompassing all evaluated samples
│   ├── Rarefaction_Curve.pdf      # Visual plot depicting species richness estimates
│   └── OTU_table.tab              # Tabular format of operational taxonomic units
├── SRR10677983/
│   ├── fastp/                     # QC & Trimmed adapter reports
│   │   ├── SRR10677983_1.fastq.gz                        # Trimmed forward reads
│   │   ├── SRR10677983_2.fastq.gz                        # Trimmed reverse reads
│   │   ├── SRR10677983_orphans.fastq.gz                  # Unpaired survival reads
│   │   ├── SRR10677983.html                              # fastp interactive HTML report
│   │   └── SRR10677983.json                              # fastp metrics in JSON format
│   ├── kraken2_reads/             # Fast Kraken2 taxonomic classification on reads
│   │   ├── SRR10677983_paired_reads_report.txt           # Standard Kraken2 taxonomic report for paired reads
│   │   ├── SRR10677983_paired_reads_output.txt           # Detailed per-read classification output
│   │   └── SRR10677983_paired_reads_biom.txt             # Kraken2 results in BIOM format
│   ├── diamond_reads/             # DIAMOND reads classification outputs
│   │   ├── SRR10677983_reads_report.daa                  # DIAMOND alignment in DAA format
│   │   ├── SRR10677983_reads_report.txt                  # Tabular alignment results (outfmt 6)
│   │   └── SRR10677983_reads_hits_with_lineage.tsv       # DIAMOND hits with full NCBI lineage
│   ├── megan_reads/               # MEGAN Last Common Ancestor taxonomic classification
│   │   └── SRR10677983_reads_summary.megan               # Extracted MEGAN summary table
│   ├── spades/                    # Assembled FASTA contigs
│   │   └── kmer_auto/             # De Novo assembly produced by SPAdes auto k-mer definition
│   │       └── contigs.fasta                             # Final assembled contigs in FASTA format
│   ├── diamond_spades_kauto/      # Annotated DIAMOND outputs vs NR
│   │   ├── SRR10677983_spades_kauto_report.txt           # Standard alignment report against NR database
│   │   ├── SRR10677983_spades_kauto_hits_with_header.tsv # Alignment hits with informative tab-separated headers
│   │   └── SRR10677983_spades_kauto_hits_with_lineage.tsv # Final comprehensive DIAMOND alignment results with NCBI lineages
│   ├── kraken2_spades_kauto/      # Kraken2 contig classification
│   │   ├── SRR10677983_spades_kauto_contig_report.txt    # Standard Kraken2 report for assembled contigs
│   │   ├── SRR10677983_spades_kauto_contig_output.txt    # Detailed per-contig classification output
│   │   └── SRR10677983_spades_kauto_contig_biom.txt      # Contig taxonomic distribution in BIOM format
│   ├── rvdb_spades_kauto/         # Optional RVDB HMM screening summaries
│   │   ├── SRR10677983_RVDB_results.tsv                  # Primary tabular results evaluated against RVDB
│   │   ├── SRR10677983_RVDB_Summary.csv                  # Aggregated CSV summary report of RVDB profile alignments
│   │   ├── SRR10677983_structural_Orfs.fasta             # ORFs translated from hits to viral structural proteins
│   │   ├── SRR10677983_structural_contigs.fasta          # Original nucleotide contigs mapped to structural proteins
│   │   ├── SRR10677983_pol_Orfs.fasta                    # ORFs translated from hits to viral polymerase proteins
│   │   └── SRR10677983_pol_contigs.fasta                 # Original nucleotide contigs mapped to polymerase proteins
│   ├── darkmatter_spades_kauto/   # ORFs & Filtered RdRp viral candidates
│   │   ├── SRR10677983_ORFs.fasta                        # Translated ORFs from unmapped contigs
│   │   ├── SRR10677983_RdRp.fev                          # Palmprint raw structural feature evaluation hits
│   │   ├── SRR10677983_RdRp.fasta                        # Sequence data of putative RdRp alignments
│   │   ├── SRR10677983_RdRp.tsv                          # Tabular summary of RdRp predictions
│   │   ├── SRR10677983_Report_Diversity.html             # Comprehensive interactive RMarkdown HTML report
│   │   ├── SRR10677983_contigs_summary.tsv               # Assembly length metrics (Mean, Quartiles, etc.)
│   │   ├── SRR10677983_contigs_summary_per_Kingdom.tsv   # Contig taxonomic breakdown at Kingdom level
│   │   ├── SRR10677983_KINGDOM_Classification_Log10.png  # Barplot comparing Kingdom-level hits
│   │   ├── SRR10677983_Viral_Diamond.tsv                 # Subset of DIAMOND hits assigned to Viruses
│   │   ├── SRR10677983_Virus_Classification_Log10.png    # Barplot of Viral Family diversity
│   │   ├── SRR10677983_viral_contigs_summary.tsv         # Assembly metrics exclusively for viral families
│   │   ├── SRR10677983_RdRp_Orfs.fasta                   # Filtered FASTA containing only target RdRp ORFs
│   │   └── SRR10677983_RdRp_contigs.fasta                # Original nucleotide contigs embedding the RdRp hit
│   └── krona_spades_kauto/        # Interactive HTML plots of relative abundances
│       ├── SRR10677983_spades_kauto_kraken2_krona.html   # Interactive Krona chart for Kraken2 taxonomy
│       └── SRR10677983_spades_kauto_diamond_krona.html   # Interactive Krona chart for DIAMOND taxonomy
├── multiqc_report.html            # Aggregated sequencing QC metrics dashboard
├── multiqc_data/                  # Exported MultiQC data directory
└── stats/                         # Per-sample and per-tool summary statistics (assembly, diamond, kraken2, darkmatter)
```

---

## 2. Execution Logs (`logs/`)

All execution logs and dynamic configuration overrides are centralized in the root `logs/` directory:

```text
logs/
├── run_overrides.yaml                   # Active run configuration generated by DeepVir.py
├── run_overrides_<YYYYMMDD_HHMMSS>.yaml # Timestamped archive of each execution's parameters
├── deepvir_desktop_gui.log              # Desktop GUI launcher and bridge diagnostics log
├── workflow/                            # Global workflow logs (e.g., MultiQC aggregation, Krona taxonomy, RVDB setup)
└── SRR10677983/                         # Per-sample rule execution logs
    ├── sra_download.log                 # fasterq-dump / pigz download log
    ├── fastp_paired.log                 # fastp trimming log
    ├── spades_kauto.log                 # SPAdes assembly log
    ├── kraken2_spades_kauto.log         # Kraken2 contig classification log
    ├── diamond_blastx_spades_kauto.log  # DIAMOND alignment log
    ├── palm_annot_spades_kauto.log      # PalmAnnot RdRp search log
    └── report_summarize_spades_kauto.log # RMarkdown diversity report generation log
```
