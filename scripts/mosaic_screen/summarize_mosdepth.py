#!/usr/bin/env python3
"""Summarise 1 Mb window depths per chromosome (exploratory mosaic aneuploidy screen).

ratio   = chromosome median window depth / median of all autosomal windows
z       = robust z of the ratio against the 22 autosomal ratios
          ((ratio - median) / (1.4826 * MAD))
f_gain  = 2*(ratio-1): cell fraction if the signal were a mosaic trisomy
f_loss  = 2*(1-ratio): cell fraction if the signal were a mosaic monosomy
Caveats: GC/mappability bias shifts some chromosomes (e.g. 16, 19, 22) a little;
sex chromosomes depend on karyotype and are not z-scored; blood WGS may miss
low-level mosaicism. Treat as descriptive, not a validated call.
"""
import gzip, sys
import numpy as np

def norm(c): return c[3:] if c.startswith("chr") else c

def main(path):
    d = {}
    with gzip.open(path, "rt") as fh:
        for l in fh:
            f = l.split("\t")
            v = float(f[3])
            if v > 0:
                d.setdefault(norm(f[0]), []).append(v)
    autos = [str(i) for i in range(1, 23)]
    ref = np.median(np.concatenate([d[c] for c in autos if c in d]))
    ratios = {c: np.median(d[c]) / ref for c in d if c in autos or c in ("X", "Y")}
    ar = np.array([ratios[c] for c in autos if c in ratios])
    med = np.median(ar); mad = 1.4826 * np.median(np.abs(ar - med)) or 1e-9
    print("chrom\tn_windows\tmedian_depth\tratio\tz_vs_autosomes\tf_gain\tf_loss")
    for c in autos + ["X", "Y"]:
        if c not in ratios: continue
        r = ratios[c]
        z = (r - med) / mad if c in autos else float("nan")
        print("%s\t%d\t%.2f\t%.4f\t%.2f\t%.3f\t%.3f" % (c, len(d[c]), np.median(d[c]), r, z, max(0, 2*(r-1)), max(0, 2*(1-r))))

if __name__ == "__main__":
    main(sys.argv[1])
