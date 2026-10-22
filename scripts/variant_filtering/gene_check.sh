ANN=/path/to/project/analysis/annovar_output/SAMPLE_ID_full.hg38_multianno.txt

for g in FANCD2 SERPINA1; do
awk -F'\t' -v g=$g '$12==g && $11=="exonic" {
  known = ($179!="." || $189!="." || $207!=".")
  af=0
  if ($179!="." && $179+0>af) af=$179+0
  if ($189!="." && $189+0>af) af=$189+0
  if ($207!="." && $207+0>af) af=$207+0
  fn = ($14!="synonymous SNV")
  n++
  if (!known) { unk++; if (fn) unkfn++ }
  else if (af<0.01 && fn) rarefn++
} END {
  print g": exonic", n+0, "| no frequency data", unk+0, "(nonsynonymous or indel:", unkfn+0, ") | rare under 1% with frequency data, nonsynonymous or indel:", rarefn+0
}' $ANN
done
