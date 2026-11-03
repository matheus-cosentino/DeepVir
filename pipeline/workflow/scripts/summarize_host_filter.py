#!/usr/bin/env python3
"""
summarize_host_filter.py
Generates the human-readable summary report (.txt), MultiQC custom table (.tsv),
and JSON metadata (.json) for the DeepVir host genome filtering (BWA) module.
Supports single or multiple target host genome accessions with per-genome breakdown.
"""

import json
import os
import re
import sys


def parse_per_host_counts(counts_files: list[str], acc_list: list[str]) -> dict[str, int]:
    """
    Parses `uniq -c` output files where each line is `<count> <host_acc>`
    and returns a dictionary mapping each target host accession to its mapped read count.
    """
    counts = {acc: 0 for acc in acc_list}
    for path in counts_files:
        if not path or not os.path.isfile(path):
            continue
        with open(path, "r", encoding="utf-8", errors="replace") as f:
            for line in f:
                parts = line.strip().split()
                if len(parts) >= 2 and parts[0].isdigit():
                    cnt = int(parts[0])
                    host_id = parts[1].strip()
                    counts[host_id] = counts.get(host_id, 0) + cnt
    return counts


def main() -> int:
    mode = sys.argv[1]  # "paired" or "unpaired"
    sample = sys.argv[2]
    accession_str = sys.argv[3]
    acc_list = [a.strip() for a in re.split(r"[,;\s]+", accession_str) if a.strip()]
    accession_label = ", ".join(acc_list) if acc_list else accession_str

    if mode == "paired":
        pe_total = int(sys.argv[4])
        pe_host = int(sys.argv[5])
        pe_unmapped_reads = int(sys.argv[6])
        orph_total = int(sys.argv[7])
        orph_host = int(sys.argv[8])
        orph_unmapped = int(sys.argv[9])
        report_path = sys.argv[10]
        tsv_path = sys.argv[11]
        json_path = sys.argv[12]
        counts_files = sys.argv[13:]

        total_reads = pe_total + orph_total
        non_host_reads = pe_unmapped_reads + orph_unmapped
        host_reads = max(0, total_reads - non_host_reads)
    else:
        total_reads = int(sys.argv[4])
        host_reads = int(sys.argv[5])
        non_host_reads = int(sys.argv[6])
        report_path = sys.argv[7]
        tsv_path = sys.argv[8]
        json_path = sys.argv[9]
        counts_files = sys.argv[10:]
        pe_total = total_reads
        pe_unmapped_reads = non_host_reads
        orph_total = 0
        orph_unmapped = 0

    host_pct = (host_reads / total_reads * 100.0) if total_reads > 0 else 0.0
    non_host_pct = (non_host_reads / total_reads * 100.0) if total_reads > 0 else 0.0

    per_host_counts = parse_per_host_counts(counts_files, acc_list)
    per_host_summary = {}
    for acc, cnt in per_host_counts.items():
        pct = (cnt / total_reads * 100.0) if total_reads > 0 else 0.0
        per_host_summary[acc] = {"reads": cnt, "percent": round(pct, 2)}

    with open(report_path, "w", encoding="utf-8") as rf:
        rf.write("========================================================\n")
        rf.write(f"DeepVir Host Genome Filtering Summary ({sample})\n")
        rf.write("========================================================\n")
        rf.write(f"Target Host Genome(s) : {accession_label}\n")
        rf.write(f"Total Input Reads     : {total_reads}\n")
        rf.write(f"Host-Associated Reads : {host_reads} ({host_pct:.2f}%)\n")
        rf.write(f"Retained Non-Host     : {non_host_reads} ({non_host_pct:.2f}%)\n")
        if len(per_host_summary) > 1 or any(v["reads"] > 0 for v in per_host_summary.values()):
            rf.write("--------------------------------------------------------\n")
            rf.write("Per-Host Genome Breakdown (Primary Mapped Reads):\n")
            for acc, stats in per_host_summary.items():
                rf.write(f"  - {acc:<20}: {stats['reads']} reads ({stats['percent']:.2f}%)\n")
        if mode == "paired":
            rf.write("--------------------------------------------------------\n")
            rf.write(f"Paired Reads (Total)  : {pe_total}\n")
            rf.write(f"Paired Reads Retained : {pe_unmapped_reads}\n")
            rf.write(f"Orphan Reads (Total)  : {orph_total}\n")
            rf.write(f"Orphan Reads Retained : {orph_unmapped}\n")
        rf.write("========================================================\n")

    with open(tsv_path, "w", encoding="utf-8") as tf:
        tf.write(
            "sample\thost_accession\ttotal_reads\thost_reads\tnon_host_reads\thost_percent\tnon_host_percent\n"
        )
        tf.write(
            f"{sample}\t{accession_label}\t{total_reads}\t{host_reads}\t{non_host_reads}\t{host_pct:.2f}\t{non_host_pct:.2f}\n"
        )

    with open(json_path, "w", encoding="utf-8") as jf:
        json.dump(
            {
                "sample": sample,
                "host_accession": accession_label,
                "host_accessions": acc_list,
                "total_reads": total_reads,
                "host_reads": host_reads,
                "non_host_reads": non_host_reads,
                "host_percent": round(host_pct, 2),
                "non_host_percent": round(non_host_pct, 2),
                "per_host_breakdown": per_host_summary,
            },
            jf,
            indent=4,
        )

    return 0


if __name__ == "__main__":
    sys.exit(main())
