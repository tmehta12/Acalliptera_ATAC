#!/usr/bin/env python3
"""Filter a tabular (outfmt 6) BLASTn of the mitochondrial genome against the assembly.

A hit is kept if it has >= 93% identity, e-value <= 1e-10 and its alignment length covers >= 75% of the
mitochondrial genome. The scaffolds in column 2 of the output are the ones removed from the BAM files.

Usage: python filter_mito_blast_hits.py genome_mt.blast mitochondrial_genome.fasta
Output: <input minus '.blast'>.filtered.blast ('no hits' if nothing passes)
"""
import os
import sys

PIDENT_MIN = 93
EVALUE_MAX = 1e-10
COVERAGE_MIN = 75  # % of the mitochondrial genome covered by the alignment


def read_fasta(path):
    seqs, name = {}, None
    with open(path) as handle:
        for line in handle:
            line = line.rstrip()
            if line.startswith(">"):
                name = line[1:]
                seqs[name] = ""
            else:
                seqs[name] += line
    return seqs


def main(blast_path, fasta_path):
    out_path = blast_path[:-len(".blast")] + ".filtered.blast"
    if os.path.getsize(blast_path) == 0:
        print("BLAST output is empty: no scaffolds match the mitochondrial genome")
        with open(out_path, "w") as out:
            out.write("no hits")
        return

    seqs = read_fasta(fasta_path)
    kept = []
    with open(blast_path) as handle:
        for line in handle:
            hit = line.rstrip("\n").split("\t")
            if float(hit[2]) < PIDENT_MIN or float(hit[10]) > EVALUE_MAX:
                continue
            if hit[0] not in seqs:
                continue
            mt_len = len(seqs[hit[0]])
            coverage = round(int(hit[3]) / mt_len * 100)
            if coverage >= COVERAGE_MIN:
                kept.append(hit + [mt_len, coverage])

    with open(out_path, "w") as out:
        if kept:
            out.write("\n".join("\t".join(str(x) for x in hit) for hit in kept) + "\n")
        else:
            out.write("no hits")
    print(f"{len(kept)} hit(s) written to {out_path}")


if __name__ == "__main__":
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])
