#!/usr/bin/env python3
"""Transparent, reproducible EPCR for the Track 1 predictions file.

For every candidate (single variant or compound-het pair):
  1. Sum the ACMG/AMP evidence points of each variant (candidates.json).
  2. Convert points to a posterior probability of pathogenicity with the
     Bayesian point framework of Tavtigian et al. 2018:
         odds_path = 350 ** (points / 8);  P = odds*prior / ((odds-1)*prior + 1)
  3. For a pair use the WEAKER allele (min), because both alleles are needed.
  4. For primary candidates multiply by the Exomiser gene-level combined score
     (0-1, phenotype + variant), read from the Exomiser genes.tsv. Secondary
     (incidental) findings are not phenotype-driven, so they use step 3 only.
  5. Clamp to [0.001, 1] because the form requires epcr in (0, 1].

Only the RANKING matters for scoring; with a single case these values are
transparent heuristics, not calibrated probabilities.

Usage:
  python compute_epcr.py --genes-tsv exomiser_result.genes.tsv \
      --candidates candidates.json --out predictions.csv
"""
import argparse, csv, json, sys

COLS = ["proband_id", "chrom_1", "pos_1", "ref_1", "alt_1", "chrom_2", "pos_2",
        "ref_2", "alt_2", "epcr", "finding_type", "notes"]

def posterior(points, prior):
    odds = 350.0 ** (points / 8.0)
    return odds * prior / ((odds - 1.0) * prior + 1.0)

def gene_scores(path):
    scores = {}
    with open(path) as fh:
        rd = csv.DictReader(fh, delimiter="\t")
        rd.fieldnames = [c.lstrip("#") for c in rd.fieldnames]
        for r in rd:
            g, s = r["GENE_SYMBOL"], float(r["EXOMISER_GENE_COMBINED_SCORE"])
            scores[g] = max(s, scores.get(g, 0.0))
    return scores

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--genes-tsv", required=True)
    ap.add_argument("--candidates", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    cfg = json.load(open(a.candidates))
    gs = gene_scores(a.genes_tsv)
    rows = []
    for c in cfg["candidates"]:
        post = [posterior(sum(v["criteria"].values()), cfg["prior"]) for v in c["variants"]]
        epcr = min(post)
        detail = "min allele posterior %.3f" % epcr
        if c["finding_type"] == "primary":
            g = gs.get(c["gene"])
            if g is None:
                sys.exit("Gene %s not found in %s" % (c["gene"], a.genes_tsv))
            epcr *= g
            detail += " x Exomiser gene score %.3f" % g
        epcr = min(1.0, max(0.001, epcr))
        v = c["variants"]
        r = [cfg["proband_id"], v[0]["chrom"], v[0]["pos"], v[0]["ref"], v[0]["alt"]]
        r += [v[1]["chrom"], v[1]["pos"], v[1]["ref"], v[1]["alt"]] if len(v) > 1 else ["", "", "", ""]
        r += [round(epcr, 4), c["finding_type"], c["notes"] + " EPCR = " + detail + "."]
        rows.append(r)
    rows.sort(key=lambda r: (r[10] != "primary", -r[9]))
    with open(a.out, "w", newline="") as fh:
        w = csv.writer(fh, lineterminator="\r\n"); w.writerow(COLS); w.writerows(rows)
    for r in rows:
        print(r[1], r[2], r[10], r[9])

if __name__ == "__main__":
    main()
