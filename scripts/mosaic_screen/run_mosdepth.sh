#!/bin/bash
#SBATCH --job-name=mosdepth
#SBATCH --error=mosdepth.%J.err
#SBATCH --output=mosdepth.%J.out
#SBATCH --partition=compute
#SBATCH --qos=common
#SBATCH --account=common
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=06:00:00
# Per-chromosome read depth in 1 Mb windows (exploratory mosaic-aneuploidy screen).
# Submit:  sbatch --export=ALL,BAM=/path/to/sample.sorted.bam,OUTDIR=/path/to/out scripts/mosaic_screen/run_mosdepth.sh
set -euo pipefail
: "${BAM:?set BAM}"; : "${OUTDIR:?set OUTDIR}"
HERE="${SLURM_SUBMIT_DIR:-$(pwd)}/scripts/mosaic_screen"
mkdir -p "$OUTDIR"
mosdepth --threads 4 --fast-mode --no-per-base --by 1000000 "$OUTDIR/depth" "$BAM"
python "$HERE/summarize_mosdepth.py" "$OUTDIR/depth.regions.bed.gz" > "$OUTDIR/chrom_depth_summary.tsv"
cat "$OUTDIR/chrom_depth_summary.tsv"
