#!/usr/bin/env python3
import sys
import gzip
import os

try:
    hit_file = snakemake.input.hit_file
    taxid_map = snakemake.input.taxid_map
    out_file = snakemake.output.ids
    log_file = snakemake.log[0]
except NameError:
    print("Must be run via Snakemake", file=sys.stderr)
    sys.exit(1)

log = open(log_file, 'w')
log.write(f"Reading hits from {hit_file}\n")

accessions_to_find = set()
hits = []
with open(hit_file, 'r') as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        parts = line.split('\t')
        if len(parts) >= 2:
            accessions_to_find.add(parts[1])
        hits.append(parts)

log.write(f"Found {len(accessions_to_find)} unique accessions.\n")
log.write(f"Scanning {taxid_map}...\n")

acc_to_taxid = {}
open_func = gzip.open if taxid_map.endswith('.gz') else open

try:
    with open_func(taxid_map, 'rt') as f:
        for line in f:
            if not accessions_to_find:
                break
            parts = line.split('\t')
            if len(parts) >= 3:
                acc = parts[0]
                acc_ver = parts[1]
                taxid = parts[2]
                
                if acc_ver in accessions_to_find:
                    acc_to_taxid[acc_ver] = taxid
                    accessions_to_find.remove(acc_ver)
                elif acc in accessions_to_find:
                    acc_to_taxid[acc] = taxid
                    accessions_to_find.remove(acc)
except Exception as e:
    log.write(f"Error reading taxonomy map: {e}\n")

log.write(f"Found {len(acc_to_taxid)} taxids.\n")

with open(out_file, 'w') as out:
    for parts in hits:
        if len(parts) >= 2:
            sseqid = parts[1]
            taxid = acc_to_taxid.get(sseqid, "NOT_FOUND")
            out.write('\t'.join(parts) + f'\t{taxid}\n')

log.write("Done.\n")
log.close()
