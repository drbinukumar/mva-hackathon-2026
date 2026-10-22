#!/usr/bin/env python3
"""Evidence panel for the two BUB1B alleles (no patient data needed).

1) Nonsense-mediated decay (NMD) prediction for the stop-gained allele using the
   50-55 nt rule: a premature stop more than ~50 nt upstream of the last exon-exon
   junction is expected to trigger NMD; a stop in the last exon or within ~50 nt of
   the last junction usually escapes it.
2) Which annotated protein domain contains each affected residue.

Inputs
  --exons    BED of the exons of ONE transcript (the one used for HGVS, NM_001211).
  --strand   + or - for the gene.
  --stop-pos genomic position (1-based) of the stop-gained variant.
  --domains  TSV: name, start_aa, end_aa (from UniProt; see docs/bub1b_domains.tsv).
  --residues amino-acid positions to look up, e.g. 737 1002.
Check the printed exon table against your transcript before trusting the result.
"""
import argparse

def load_exons(path, strand):
    ex = []
    for line in open(path):
        if not line.strip() or line.startswith(("#", "track", "browser")):
            continue
        f = line.split()
        ex.append((int(f[1]) + 1, int(f[2])))          # 1-based inclusive
    ex.sort(reverse=(strand == "-"))                    # transcript order
    return ex

def nmd(exons, strand, pos):
    idx = next((i for i, (s, e) in enumerate(exons) if s <= pos <= e), None)
    if idx is None:
        return None, "position is not inside any exon in the BED"
    last = len(exons) - 1
    s, e = exons[idx]
    to_exon_end = (e - pos) if strand == "+" else (pos - s)
    dist = to_exon_end + sum(b - a + 1 for a, b in exons[idx + 1:last])
    if idx == last:
        return 0, "in the LAST exon: predicted to ESCAPE NMD (truncated protein possible)"
    if idx == last - 1 and to_exon_end < 50:
        return to_exon_end, "within 50 nt of the last junction: predicted to ESCAPE NMD"
    return dist, "%d nt upstream of the last exon-exon junction (> 50-55 nt): predicted to TRIGGER NMD" % dist

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--exons", required=True)
    ap.add_argument("--strand", choices=["+", "-"], required=True)
    ap.add_argument("--stop-pos", type=int, required=True)
    ap.add_argument("--domains")
    ap.add_argument("--residues", type=int, nargs="*", default=[])
    a = ap.parse_args()
    ex = load_exons(a.exons, a.strand)
    print("Exons in transcript order (1-based): %d exons" % len(ex))
    for i, (s, e) in enumerate(ex, 1):
        print("  exon %2d  %d-%d  (%d bp)" % (i, s, e, e - s + 1))
    d, msg = nmd(ex, a.strand, a.stop_pos)
    print("\nStop-gained variant at %d: %s" % (a.stop_pos, msg))
    if a.domains:
        doms = [l.rstrip("\n").split("\t") for l in open(a.domains) if l.strip() and not l.startswith("#")]
        for r in a.residues:
            hit = [n for n, s, e in doms if int(s) <= r <= int(e)]
            print("Residue %d: %s" % (r, ", ".join(hit) if hit else "no annotated domain in the supplied table"))

if __name__ == "__main__":
    main()
