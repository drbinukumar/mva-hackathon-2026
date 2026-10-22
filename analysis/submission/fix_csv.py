import csv, sys
proband, team = sys.argv[1], sys.argv[2]
rows = [r for r in csv.reader(open("binukumar_bub1b-compound-het.csv", newline="")) if r]
for r in rows[1:]:
    r[0] = proband
    for i in (1, 5):
        if r[i] and not r[i].startswith("chr"):
            r[i] = "chr" + r[i]
out = f"{team}_bub1b-compound-het.csv"
csv.writer(open(out, "w", newline=""), lineterminator="\n").writerows(rows)
print("wrote", out)
