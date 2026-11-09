#!/usr/bin/env python3
import sys
import re
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

DB_TAGS = {
    "gb", "ref", "emb", "dbj", "sp", "tr", "gi", "pdb", "prf", "pir",
    "tpg", "tpe", "tpd", "gnl", "lcl", "pat", "bbs", "genbank", "refseq",
}


def _variants(token):
    """Accession with version, and without version (e.g. WP_000001.1 -> WP_000001)."""
    token = token.strip().strip('"\'')
    if not token:
        return []
    out = [token]
    if '.' in token:
        out.append(token.rsplit('.', 1)[0])
    return out


def candidate_accessions(sseqid):
    """
    Returns lookup keys for a DIAMOND subject ID (column 2), in priority order.
    Supports any database ID layout:
      - RVDB:        acc|p_bank|p_acc|n_bank|n_acc|descr[sp]  -> p_acc first, then acc
      - NCBI style:  ref|WP_000001.1|, gb|AAA123.1|, gi|123|ref|NP_1.1|
      - UniProt:     sp|P12345|NAME_ECOLI, tr|A0A001|...
      - Plain:       WP_000001.1  /  P12345  (protein-only custom DBs)
      - Any '|', ':', ';', ',' separated list of IDs (each token tried)
    The original ID is always tried (first for plain IDs, last for composite IDs).
    """
    sseqid = sseqid.strip()
    if not any(sep in sseqid for sep in '|:;, '):
        return _variants(sseqid)

    fields = [f for f in re.split(r'[|:;,\s]+', sseqid) if f]
    cands = []

    # RVDB layout: >acc|p_bank|p_acc|n_bank|n_acc|descr -> protein accession is field 3
    pipe_fields = sseqid.split('|')
    if len(pipe_fields) >= 5 and pipe_fields[1].lower() in DB_TAGS | {"", "na"}:
        cands += _variants(pipe_fields[2])
        cands += _variants(pipe_fields[0])
    # Generic: every token that is not a bank tag or a pure number (gi numbers are not in the map)
    for tok in fields:
        if tok.lower() in DB_TAGS or tok.isdigit():
            continue
        cands += _variants(tok)
    cands += _variants(sseqid)

    seen = set()
    return [c for c in cands if c and not (c in seen or seen.add(c))]


accessions_to_find = set()
hits = []
with open(hit_file, 'r') as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        parts = line.split('\t')
        if len(parts) >= 2:
            accessions_to_find.update(candidate_accessions(parts[1]))
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
                
                taxid = taxid.strip()
                if acc_ver in accessions_to_find:
                    acc_to_taxid[acc_ver] = taxid
                    accessions_to_find.remove(acc_ver)
                if acc in accessions_to_find:
                    acc_to_taxid[acc] = taxid
                    accessions_to_find.remove(acc)
except Exception as e:
    log.write(f"Error reading taxonomy map: {e}\n")

log.write(f"Found {len(acc_to_taxid)} taxids.\n")

with open(out_file, 'w') as out:
    for parts in hits:
        if len(parts) >= 2:
            sseqid = parts[1]
            taxid = next(
                (acc_to_taxid[c] for c in candidate_accessions(sseqid) if c in acc_to_taxid),
                "NOT_FOUND",
            )
            out.write('\t'.join(parts) + f'\t{taxid}\n')

log.write("Done.\n")
log.close()
