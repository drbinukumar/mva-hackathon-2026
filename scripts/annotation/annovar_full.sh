#!/bin/bash
#SBATCH --job-name=annovar_full
#SBATCH --error=/path/to/project/logs/annovar_full.%J.err
#SBATCH --output=/path/to/project/logs/annovar_full.%J.out
#SBATCH --partition=compute
#SBATCH --qos=common
#SBATCH --account=common
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=20
#SBATCH --cpus-per-task=1
#SBATCH --time=24:00:00

cd /path/to/project/analysis/annovar_output

perl /path/to/tools/annovar_2025/annovar/table_annovar.pl \
--thread 20 /path/to/project/data/SAMPLE_ID_FLOWCELL_ID.vcf.gz \
/path/to/tools/annovar_2025/annovar/humandb \
-buildver hg38 \
-out SAMPLE_ID_full \
-remove \
-protocol refGeneWithVer,refGene,avsnp151,dbnsfp47a,clinvar_20250409,1000g2015aug_all,exac03,esp6500siv2_all,gnomad41_exome,gnomad41_genome,gme,dbscsnv11,regsnpintron,revel \
-operation g,g,f,f,f,f,f,f,f,f,f,f,f,f \
-otherinfo -nastring . \
-polish -vcfinput

echo "Full annotation completed: $(date)"
