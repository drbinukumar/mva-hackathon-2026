#!/bin/bash
#SBATCH --job-name=var_filter
#SBATCH --error=/path/to/project/logs/var_filter.%J.err
#SBATCH --output=/path/to/project/logs/var_filter.%J.out
#SBATCH --partition=compute
#SBATCH --qos=common
#SBATCH --account=common
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=06:00:00

module load bcftools-1.15

var_prior=/path/to/project/analysis/var_filtering
data_dir=/path/to/project/analysis/annovar_output
gene_panel=/path/to/project/analysis/gene_panel/full_gene_list.txt

cd $var_prior

# Reheader (add sample name column) - header naturally preserved here
{ paste <(head -n1 $data_dir/SAMPLE_ID_full.hg38_multianno.txt | cut -f1-252) <(bcftools query -l $data_dir/SAMPLE_ID_full.hg38_multianno.vcf | paste -sd'\t' -); tail -n +2 $data_dir/SAMPLE_ID_full.hg38_multianno.txt; } > SAMPLE_ID.hg38_multianno.reheaded.txt

# a. Exonic/UTR/splicing + AF <5% filter (header preserved via NR==1 in awk chain)
egrep -w 'Chr|exonic|UTR3|UTR5|splicing|AAChange.refGene' SAMPLE_ID.hg38_multianno.reheaded.txt | awk -F'\t' 'NR==1 || !($14 == "synonymous SNV" || $14 == "unknown")' | awk -F'\t' '(NR==1) || ($179 < 0.05) && ($189 < 0.05) && ($207 < 0.05)' > SAMPLE_ID_0.05.txt

# b. ClinVar pathogenic - header preserved via NR==1
awk -F'\t' 'NR == 1 || tolower($170) ~ /pathogenic/' SAMPLE_ID_0.05.txt | awk -F'\t' 'NR==1 {print $0"\tCategory"; next} {print $0 "\tClinVar"}' > clinvar.txt

# c. SIFT/PolyPhen/CADD (PP3 logic) - header preserved in each intermediate
awk -F'\t' 'NR==1 { print } NR>1 && ($17 == "." || $26 == ".") { $17 = 0; $26 = 0; } ($17 <= 0.05 && $26 >= 0.447) { print }' SAMPLE_ID_0.05.txt > 1.txt
awk -F'\t' 'NR==1 { print } NR>1 && ($17 == "." || $108 == ".") { $17 = 0; $108 = 0; } ($17 <= 0.05 && $108 >= 20) { print }' SAMPLE_ID_0.05.txt > 2.txt
awk -F'\t' 'NR==1 { print } NR>1 && ($108 == "." || $26 == ".") { $108 = 0; $26 = 0; } ($108 >= 20 && $26 >= 0.447) { print }' SAMPLE_ID_0.05.txt > 3.txt
# Combine and dedupe, header preserved by keeping first line separately
head -n1 1.txt > pp3.txt
tail -n +2 -q 1.txt 2.txt 3.txt | sort -k2n -u | awk -F'\t' '{print $0 "\tPP3"}' >> pp3.txt
# Fix header to include Category label
sed -i '1s/$/\tCategory/' pp3.txt

# d. LoF variants - header preserved
egrep -w 'Chr|frameshift deletion|frameshift insertion|stopgain|startloss|splicing' SAMPLE_ID_0.05.txt | awk -F'\t' 'NR==1 {print $0"\tCategory"; next} {print $0 "\tP-pLoF"}' > lof.txt

# e. Intronic variants flagged by regsnp-intron - header preserved
egrep -w 'Chr|intronic' SAMPLE_ID.hg38_multianno.reheaded.txt | awk -F'\t' 'NR==1 {print $0"\tCategory"; next} $237 ~ /disease-causing/ {print $0 "\tIntronicSplice"}' > intronic_regsnp.txt

# Merge exonic evidence, remove benign/conflicting - header explicitly written first
head -n1 lof.txt | cut -f1-252 > PCL_header.txt
echo -e "Prediction" | paste PCL_header.txt - > PCL.txt
tail -n +2 -q lof.txt clinvar.txt pp3.txt | awk -F'\t' 'tolower($170) !~ /benign|conflicting/ {key=""; for(i=1;i<NF;i++) key=key $i OFS; values[key]=(key in values) ? values[key]","$NF:$NF} END {for(k in values) print k values[k]}' OFS='\t' >> PCL.txt

# Novel-in-databases check - header preserved
awk 'BEGIN {FS=OFS="\t"} NR==1 {print $0, "Databases_Status"; next} {status = ($179=="." && $207=="." && $189=="." && $170==".") ? "Novel" : "."; print $0, status}' PCL.txt > SAMPLE_ID_final.txt

# Homozygous variants - header preserved via 'Chr' match in egrep
egrep 'Chr|1/1|1\|1' SAMPLE_ID_final.txt > SAMPLE_ID_homo_var.txt

# Gene panel filtering - header explicitly written first, then matches appended
head -n1 SAMPLE_ID_final.txt > mva_panel_candidates.txt
grep -w -f $gene_panel SAMPLE_ID_final.txt >> mva_panel_candidates.txt

head -n1 intronic_regsnp.txt > mva_panel_intronic_candidates.txt
grep -w -f $gene_panel intronic_regsnp.txt >> mva_panel_intronic_candidates.txt

echo "Filtering completed: $(date)"
echo "Row counts (including header):"
wc -l SAMPLE_ID_0.05.txt clinvar.txt lof.txt pp3.txt PCL.txt SAMPLE_ID_final.txt mva_panel_candidates.txt mva_panel_intronic_candidates.txt
