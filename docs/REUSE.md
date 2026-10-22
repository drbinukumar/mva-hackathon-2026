# Reusing this pipeline on another ultra-rare case

This repository contains code and small configuration files only. It never contains
patient data (VCF, BAM, FASTQ, annotated tables). Keep those in a work directory
outside the repository (`WORKDIR`, see `config/paths.env.example`).

## 1. Set up
```
conda env create -f environment.yml && conda activate mva_track1
cp config/paths.env.example config/paths.env   # edit; this file is git-ignored
source config/paths.env
```
ANNOVAR and Exomiser are installed separately (see the comments in `environment.yml`).

## 2. Order of the steps
1. `scripts/snakemake/Snakefile` : alignment and variant-calling workflow
2. `scripts/annotation/annovar_full.sh` : ANNOVAR annotation of the VCF
3. `scripts/exomiser/run_exomiser.sh` : phenotype-driven ranking with the case's HPO terms
4. `scripts/variant_filtering/` : rarity / impact filters and the chromosomal-instability gene check
5. `scripts/sv_calling/run_manta.sh` : structural variants (Manta)
6. `scripts/mosaic_screen/` : allele-balance screen and per-chromosome depth (mosdepth)
7. `scripts/evidence/bub1b_allele_evidence.py` : NMD and domain evidence for a candidate gene
8. `scripts/epcr/compute_epcr.py` : turns ACMG evidence plus the Exomiser gene score into the
   predictions CSV (edit `candidates.json` for a new case)

## 3. Adapting to another disease
* Replace the gene list in `docs/` with the genes of the disease family you care about.
* Replace the HPO terms in the Exomiser analysis file with those of the new case.
* Edit `candidates.json` with the candidate variants and the ACMG criteria you judged to apply.
* Cluster users: Slurm headers in the scripts use the `cancerq2` partition; change them for your site.

## 4. Limits
Single-case analysis. EPCR values are transparent heuristics (ACMG points converted with the
Tavtigian Bayesian framework, scaled by the Exomiser gene score), not calibrated probabilities.
Phase of compound-heterozygous variants needs parental or long-read data.
