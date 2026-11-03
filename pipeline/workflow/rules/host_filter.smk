# Host Genome Download, BWA Indexing, and Multi-Host Read Depletion Rules

import re as _re

_RAW_HOST_ACC = config.get("host_filter", {}).get("accession", "GCA_022682495.2")
if isinstance(_RAW_HOST_ACC, (list, tuple)):
    _RAW_HOST_STR = ",".join(str(x) for x in _RAW_HOST_ACC)
else:
    _RAW_HOST_STR = str(_RAW_HOST_ACC or "GCA_022682495.2")

_HOST_ACC_LIST = []
for _tok in _re.split(r"[,;\s]+", _RAW_HOST_STR.strip()):
    _clean = _tok.strip()
    if _clean and _clean not in _HOST_ACC_LIST:
        _HOST_ACC_LIST.append(_clean)
if not _HOST_ACC_LIST:
    _HOST_ACC_LIST = ["GCA_022682495.2"]

_HOST_SET_ID = "__".join(_HOST_ACC_LIST)
_HOST_LABEL = ", ".join(_HOST_ACC_LIST)
_HOST_INDEX_DIR = os.path.join(RESOURCES_DIR, "host", "index", _HOST_SET_ID)
_HOST_COMBINED_FASTA = os.path.join(_HOST_INDEX_DIR, "host_reference.fna")
_HOST_SUMMARY_SCRIPT = _os.path.join(_PROJECT_ROOT, "workflow", "scripts", "summarize_host_filter.py")


rule download_host_genome:
    wildcard_constraints:
        host_acc = r"[A-Za-z0-9_.\-]+"
    output:
        genome = os.path.join(RESOURCES_DIR, "host", "genomes", "{host_acc}", "{host_acc}_download_host_genome.fna")
    params:
        host_dir = os.path.join(RESOURCES_DIR, "host", "genomes", "{host_acc}"),
        zip_file = os.path.join(RESOURCES_DIR, "host", "genomes", "{host_acc}", "{host_acc}_dataset.zip")
    conda:
        HOST_FILTER
    resources:
        threads = config.get("simple_Task_threads", 4)
    log:
        os.path.join(LOG_DIR, "workflow", "{host_acc}_download_host_genome.log")
    shell:
        """
        mkdir -p "{params.host_dir}" "$(dirname "{log}")"
        rm -rf "{params.host_dir}/ncbi_dataset" "{params.zip_file}"
        datasets download genome accession "{wildcards.host_acc}" --include genome --filename "{params.zip_file}" > "{log}" 2>&1
        unzip -o "{params.zip_file}" -d "{params.host_dir}" >> "{log}" 2>&1
        sed "s/^>/>{wildcards.host_acc}|/" "{params.host_dir}"/ncbi_dataset/data/*/*_genomic.fna > "{output.genome}" 2>> "{log}"
        rm -rf "{params.host_dir}/ncbi_dataset" "{params.zip_file}" "{params.host_dir}/README.md" "{params.host_dir}/md5sum.txt"
        """


rule bwa_index_host:
    input:
        genomes = expand(
            os.path.join(RESOURCES_DIR, "host", "genomes", "{host_acc}", "{host_acc}_download_host_genome.fna"),
            host_acc=_HOST_ACC_LIST
        )
    output:
        combined = _HOST_COMBINED_FASTA,
        idx = multiext(_HOST_COMBINED_FASTA, ".amb", ".ann", ".bwt", ".pac", ".sa")
    params:
        index_dir = _HOST_INDEX_DIR
    conda:
        HOST_FILTER
    resources:
        threads = config.get("simple_Task_threads", 4)
    log:
        os.path.join(LOG_DIR, "workflow", f"{_HOST_SET_ID}_bwa_index_host.log")
    shell:
        """
        mkdir -p "{params.index_dir}" "$(dirname "{log}")"
        cat {input.genomes} > "{output.combined}" 2> "{log}"
        bwa index "{output.combined}" >> "{log}" 2>&1
        """


rule host_filter_paired:
    input:
        r1 = os.path.join(OUT_DIR, "{sample}", "fastp", "{sample}_fastp_paired_1.fastq.gz"),
        r2 = os.path.join(OUT_DIR, "{sample}", "fastp", "{sample}_fastp_paired_2.fastq.gz"),
        orphans = os.path.join(OUT_DIR, "{sample}", "fastp", "{sample}_fastp_paired_orphans.fastq.gz"),
        genome = _HOST_COMBINED_FASTA,
        idx = multiext(_HOST_COMBINED_FASTA, ".amb", ".ann", ".bwt", ".pac", ".sa")
    output:
        r1 = os.path.join(OUT_DIR, "{sample}", "host_filter", "{sample}_host_filter_paired_1.fastq.gz"),
        r2 = os.path.join(OUT_DIR, "{sample}", "host_filter", "{sample}_host_filter_paired_2.fastq.gz"),
        orphans = os.path.join(OUT_DIR, "{sample}", "host_filter", "{sample}_host_filter_paired_orphans.fastq.gz"),
        report = os.path.join(OUT_DIR, "{sample}", "host_filter", "{sample}_host_filter_paired_summary.txt"),
        tsv = os.path.join(OUT_DIR, "{sample}", "stats", "{sample}_host_filter_paired_stats.tsv"),
        json = os.path.join(OUT_DIR, "{sample}", "stats", "{sample}_host_filter_paired_stats.json")
    params:
        accession = _HOST_LABEL,
        out_dir = os.path.join(OUT_DIR, "{sample}", "host_filter"),
        stats_dir = os.path.join(OUT_DIR, "{sample}", "stats"),
        summary_script = _HOST_SUMMARY_SCRIPT
    conda:
        HOST_FILTER
    resources:
        threads = config.get("simple_Task_threads", 4)
    log:
        os.path.join(LOG_DIR, "{sample}", "{sample}_host_filter_paired.log")
    benchmark:
        os.path.join(OUT_DIR, "benchmarks", "host_filter_paired", "{sample}_host_filter_paired.tsv")
    shell:
        """
        mkdir -p "{params.out_dir}" "{params.stats_dir}" "$(dirname "{log}")"

        PE_BAM="{params.out_dir}/{wildcards.sample}_pe.bam"
        ORPH_BAM="{params.out_dir}/{wildcards.sample}_orphans.bam"
        PE_COUNTS="{params.out_dir}/{wildcards.sample}_pe_host_counts.tmp"
        ORPH_COUNTS="{params.out_dir}/{wildcards.sample}_orph_host_counts.tmp"

        bwa mem -t {resources.threads} "{input.genome}" "{input.r1}" "{input.r2}" 2> "{log}" \
            | samtools view -@ {resources.threads} -b -F 2304 -o "$PE_BAM" - 2>> "{log}"

        PE_TOTAL=$(samtools view -@ {resources.threads} -c "$PE_BAM")
        PE_HOST=$(samtools view -@ {resources.threads} -c -F 4 "$PE_BAM")
        samtools view -@ {resources.threads} -F 4 "$PE_BAM" | cut -f 3 | cut -d'|' -f1 | sort | uniq -c > "$PE_COUNTS" || true

        samtools fastq -@ {resources.threads} -f 12 -F 2304 \
            -1 "{output.r1}" -2 "{output.r2}" -0 /dev/null -s /dev/null -n "$PE_BAM" 2>> "{log}"

        PE_UNMAPPED=$(samtools view -@ {resources.threads} -c -f 12 "$PE_BAM")
        rm -f "$PE_BAM"

        ORPH_TOTAL=0
        ORPH_HOST=0
        ORPH_UNMAPPED=0
        : > "$ORPH_COUNTS"
        if [ -s "{input.orphans}" ]; then
            bwa mem -t {resources.threads} "{input.genome}" "{input.orphans}" 2>> "{log}" \
                | samtools view -@ {resources.threads} -b -F 2304 -o "$ORPH_BAM" - 2>> "{log}"
            ORPH_TOTAL=$(samtools view -@ {resources.threads} -c "$ORPH_BAM")
            ORPH_HOST=$(samtools view -@ {resources.threads} -c -F 4 "$ORPH_BAM")
            ORPH_UNMAPPED=$(samtools view -@ {resources.threads} -c -f 4 "$ORPH_BAM")
            samtools view -@ {resources.threads} -F 4 "$ORPH_BAM" | cut -f 3 | cut -d'|' -f1 | sort | uniq -c > "$ORPH_COUNTS" || true
            samtools fastq -@ {resources.threads} -f 4 -F 2304 "$ORPH_BAM" | gzip -c > "{output.orphans}" 2>> "{log}"
            rm -f "$ORPH_BAM"
        else
            gzip -c /dev/null > "{output.orphans}"
        fi

        python3 "{params.summary_script}" paired \
            "{wildcards.sample}" "{params.accession}" \
            "$PE_TOTAL" "$PE_HOST" "$PE_UNMAPPED" \
            "$ORPH_TOTAL" "$ORPH_HOST" "$ORPH_UNMAPPED" \
            "{output.report}" "{output.tsv}" "{output.json}" \
            "$PE_COUNTS" "$ORPH_COUNTS" >> "{log}" 2>&1

        rm -f "$PE_COUNTS" "$ORPH_COUNTS"
        """


rule host_filter_unpaired:
    input:
        r1 = os.path.join(OUT_DIR, "{sample}", "fastp", "{sample}_fastp_unpaired.fastq.gz"),
        genome = _HOST_COMBINED_FASTA,
        idx = multiext(_HOST_COMBINED_FASTA, ".amb", ".ann", ".bwt", ".pac", ".sa")
    output:
        r1 = os.path.join(OUT_DIR, "{sample}", "host_filter", "{sample}_host_filter_unpaired.fastq.gz"),
        report = os.path.join(OUT_DIR, "{sample}", "host_filter", "{sample}_host_filter_unpaired_summary.txt"),
        tsv = os.path.join(OUT_DIR, "{sample}", "stats", "{sample}_host_filter_unpaired_stats.tsv"),
        json = os.path.join(OUT_DIR, "{sample}", "stats", "{sample}_host_filter_unpaired_stats.json")
    params:
        accession = _HOST_LABEL,
        out_dir = os.path.join(OUT_DIR, "{sample}", "host_filter"),
        stats_dir = os.path.join(OUT_DIR, "{sample}", "stats"),
        summary_script = _HOST_SUMMARY_SCRIPT
    conda:
        HOST_FILTER
    resources:
        threads = config.get("simple_Task_threads", 4)
    log:
        os.path.join(LOG_DIR, "{sample}", "{sample}_host_filter_unpaired.log")
    benchmark:
        os.path.join(OUT_DIR, "benchmarks", "host_filter_unpaired", "{sample}_host_filter_unpaired.tsv")
    shell:
        """
        mkdir -p "{params.out_dir}" "{params.stats_dir}" "$(dirname "{log}")"

        SE_BAM="{params.out_dir}/{wildcards.sample}_se.bam"
        SE_COUNTS="{params.out_dir}/{wildcards.sample}_se_host_counts.tmp"

        bwa mem -t {resources.threads} "{input.genome}" "{input.r1}" 2> "{log}" \
            | samtools view -@ {resources.threads} -b -F 2304 -o "$SE_BAM" - 2>> "{log}"

        TOTAL_READS=$(samtools view -@ {resources.threads} -c "$SE_BAM")
        HOST_READS=$(samtools view -@ {resources.threads} -c -F 4 "$SE_BAM")
        UNMAPPED_READS=$(samtools view -@ {resources.threads} -c -f 4 "$SE_BAM")
        samtools view -@ {resources.threads} -F 4 "$SE_BAM" | cut -f 3 | cut -d'|' -f1 | sort | uniq -c > "$SE_COUNTS" || true

        samtools fastq -@ {resources.threads} -f 4 -F 2304 "$SE_BAM" | gzip -c > "{output.r1}" 2>> "{log}"
        rm -f "$SE_BAM"

        python3 "{params.summary_script}" unpaired \
            "{wildcards.sample}" "{params.accession}" \
            "$TOTAL_READS" "$HOST_READS" "$UNMAPPED_READS" \
            "{output.report}" "{output.tsv}" "{output.json}" \
            "$SE_COUNTS" >> "{log}" 2>&1

        rm -f "$SE_COUNTS"
        """
