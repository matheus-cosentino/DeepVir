import sys
import os
from Bio import SeqIO
from Bio.Seq import Seq

def translate_frames(dna_seq_str, min_aa=60, stitch=True):
    """
    Translates a nucleotide sequence in all 6 frames.
    Finds open reading frames (bounded by stop codons or ends).
    If stitch=True, stitches adjacent ORFs on the same strand that exhibit frameshifts (1-2 nt indel gaps).
    """
    results = []
    seq_len = len(dna_seq_str)
    if seq_len < 180:
        return results

    bio_seq = Seq(dna_seq_str)
    rev_seq = bio_seq.reverse_complement()

    # Track segments: list of dicts with strand, frame_offset, start_aa, end_aa, aa_seq
    segments = []

    for strand, s_obj, is_rev in [("+", bio_seq, False), ("-", rev_seq, True)]:
        for f_offset in (0, 1, 2):
            sub_len = (len(s_obj) - f_offset) // 3 * 3
            if sub_len < 180:
                continue
            sub_seq = s_obj[f_offset:f_offset + sub_len]
            try:
                aa_str = str(sub_seq.translate(table=1))
            except Exception:
                continue

            pos = 0
            for chunk in aa_str.split("*"):
                c_len = len(chunk)
                if c_len >= min_aa:
                    segments.append({
                        "strand": strand,
                        "is_rev": is_rev,
                        "frame": f_offset,
                        "start_aa": pos,
                        "end_aa": pos + c_len,
                        "aa": chunk,
                        "len": c_len
                    })
                pos += c_len + 1

    # Add all unstitched single ORFs
    seen_seqs = set()
    for s in segments:
        if s["aa"] not in seen_seqs:
            seen_seqs.add(s["aa"])
            results.append((f"orf_{s['strand']}_f{s['frame']}_{s['start_aa']}_{s['end_aa']}", s["aa"]))

    # Frameshift Stitching (Strategy 2)
    if stitch:
        for strand in ("+", "-"):
            strand_segs = [s for s in segments if s["strand"] == strand]
            strand_segs.sort(key=lambda x: x["start_aa"])

            for i in range(len(strand_segs)):
                curr = strand_segs[i]
                stitched_aa = curr["aa"]
                last_end = curr["end_aa"]
                stitched_frames = [str(curr["frame"])]
                last_frame = curr["frame"]

                for j in range(i + 1, len(strand_segs)):
                    nxt = strand_segs[j]
                    gap = nxt["start_aa"] - last_end
                    # Frameshift signature: adjacent on same strand in DIFFERENT reading frames
                    # Gap between -5 aa (overlap) and 25 aa
                    if -5 <= gap <= 25 and nxt["frame"] != last_frame:
                        stitched_aa += "X" + nxt["aa"]
                        last_end = nxt["end_aa"]
                        last_frame = nxt["frame"]
                        stitched_frames.append(str(nxt["frame"]))

                if len(stitched_frames) > 1 and len(stitched_aa) >= min_aa:
                    if stitched_aa not in seen_seqs:
                        seen_seqs.add(stitched_aa)
                        tag = f"stitched_{strand}_frames_{'_'.join(stitched_frames)}"
                        results.append((tag, stitched_aa))

    return results

def main():
    try:
        fasta_file = str(snakemake.input.fasta)
        output_file = str(snakemake.output.orfs)
        sample_name = snakemake.wildcards.sample
        frameshift_mode = bool(snakemake.params.frameshift_mode)
        min_aa = int(getattr(snakemake.params, "min_aa", 60))
        log_file = snakemake.log[0]
    except NameError:
        print("ERROR: Designed to be run via Snakemake.", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"ERROR accessing Snakemake variables: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"Starting frameshift-tolerant ORF prediction for sample {sample_name}", file=sys.stderr)
    print(f"  Input FASTA: {fasta_file}", file=sys.stderr)
    print(f"  Output ORFs: {output_file}", file=sys.stderr)
    print(f"  Frameshift Stitching: {frameshift_mode}, Min AA: {min_aa}", file=sys.stderr)

    if not os.path.exists(fasta_file) or os.path.getsize(fasta_file) == 0:
        print("[INFO] Input FASTA is empty. Writing empty ORF output.", file=sys.stderr)
        open(output_file, "w").close()
        return

    count_in = 0
    count_orfs = 0

    with open(output_file, "w", encoding="utf-8") as out_h:
        for record in SeqIO.parse(fasta_file, "fasta"):
            count_in += 1
            seq_str = str(record.seq).upper()
            orfs = translate_frames(seq_str, min_aa=min_aa, stitch=frameshift_mode)
            for idx, (tag, aa_seq) in enumerate(orfs, 1):
                count_orfs += 1
                header = f">{record.id}_orf{idx}_{tag}_len{len(aa_seq)}\n"
                out_h.write(header)
                out_h.write(f"{aa_seq}\n")

    print(f"Processed {count_in} nucleotide sequences, generated {count_orfs} ORFs/stitched peptides.", file=sys.stderr)

if __name__ == "__main__":
    main()

