#!/bin/bash
#SBATCH --job-name=exomiser_run
#SBATCH --error=/path/to/project/logs/exomiser_run.%J.err
#SBATCH --output=/path/to/project/logs/exomiser_run.%J.out
#SBATCH --partition=compute
#SBATCH --qos=common
#SBATCH --account=common
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=06:00:00

cd /path/to/project/tools/exomiser/exomiser-cli-14.0.0

java -Xms8g -Xmx28g -jar exomiser-cli-14.0.0.jar \
  --analysis /path/to/project/analysis/mva_analysis.yml

echo "Exomiser run completed: $(date)"
