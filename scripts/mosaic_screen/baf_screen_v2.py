import argparse
import gzip
import statistics as st
from collections import defaultdict

parser = argparse.ArgumentParser(
    description="Summarize autosomal heterozygous-site allele balance from a VCF."
)
parser.add_argument("vcf", help="Input VCF or gzip-compressed VCF")
args = parser.parse_args()
vcf_path = args.vcf


autosomes = {str(i) for i in range(1, 23)}

baf = defaultdict(list)
dp = defaultdict(list)

open_vcf = gzip.open if vcf_path.endswith((".gz", ".bgz", ".bgzf")) else open

with open_vcf(vcf_path, "rt") as f:
    for line in f:
        if line.startswith("#"):
            continue
        x = line.rstrip("\n").split("\t")
        if x[0] not in autosomes or x[6] != "PASS":
            continue
        keys = x[8].split(":")
        vals = x[9].split(":")
        if "GT" not in keys or "AD" not in keys:
            continue
        if vals[keys.index("GT")] not in ("0/1", "0|1", "1|0"):
            continue
        try:
            r, a = map(int, vals[keys.index("AD")].split(","))
        except ValueError:
            continue
        tot = r + a
        if tot < 10:
            continue
        baf[x[0]].append(a / tot)
        dp[x[0]].append(tot)

gmed = st.median([d for v in dp.values() for d in v])

print("chrom\tn_hets\tdepth_ratio\tmean_baf\tfolded_baf\tfrac_below_0.4\tfrac_above_0.6")
for c in sorted(baf, key=int):
    b = baf[c]
    n = len(b)
    print("\t".join([
        c, str(n),
        f"{st.median(dp[c]) / gmed:.3f}",
        f"{sum(b) / n:.4f}",
        f"{sum(abs(v - 0.5) for v in b) / n:.4f}",
        f"{sum(v < 0.4 for v in b) / n:.3f}",
        f"{sum(v > 0.6 for v in b) / n:.3f}",
    ]))
