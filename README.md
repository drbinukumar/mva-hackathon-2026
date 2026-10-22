# MVA Hackathon 2026 - Track 1 (variant prediction)

Code, configuration and small summary outputs for a Track 1 entry to the
"Rare Disease, Real Kid: MVA Hackathon 2026". No patient-level data
(FASTQ, BAM, VCF, annotated tables) is included in this repository.

## Candidate
A compound-heterozygous candidate in **BUB1B** (chr15, GRCh38):
- chr15:40209701 T>G, stop-gained (c.2210T>G, p.Leu737Ter). Predicted to trigger
  nonsense-mediated decay (747 nt upstream of the last exon junction).
- chr15:40220612 T>G, missense (c.3006T>G, p.Asn1002Lys) in the last exon. Classified VUS.
- Phase (trans configuration) could not be confirmed from the proband's data alone.
A pharmacogenomic variant (DPYD rs1801160) is listed as an incidental secondary finding.

## Repository layout
- `submission/`: predictions CSV and the script that builds it
- `scripts/annotation/`: ANNOVAR annotation workflow
- `scripts/exomiser/`: Exomiser analysis launcher
- `scripts/variant_filtering/`: annotation-based filtering and gene checks
- `scripts/snakemake/`: alignment and phasing workflow
- `scripts/sv_calling/`: Manta launcher
- `scripts/mosaic_screen/`: allele-balance and per-chromosome depth screens
- `scripts/evidence/`: BUB1B allele evidence panel (NMD, domains) and AlphaGenome query
- `scripts/epcr/`: EPCR calculation
- `config/`: configuration templates (copy `paths.env.example` to `paths.env`; not committed)
- `docs/`: gene lists, BUB1B exon BED, reuse guide, plain-language summary
- `results_summary/`: small aggregate outputs only
- `environment.yml`: software environment

## How the predictions were produced
1. Organiser VCF annotated with ANNOVAR; phenotype-driven ranking with Exomiser 14.0.0
   using the proband's HPO terms (BUB1B ranked 1 of 377 genes).
2. Rarity and impact filters, plus a chromosomal-instability gene list check.
3. Evidence panel for the two alleles (`scripts/evidence/`).
4. EPCR computed by `scripts/epcr/compute_epcr.py` from hand-assigned ACMG evidence points
   (Bayesian point framework) scaled by the Exomiser gene score. Only ranking matters for
   scoring; values are heuristics, not calibrated probabilities.

## Exploratory analyses (descriptive, not validated calls)
- Per-chromosome read depth (mosdepth, 1 Mb windows): no clear mosaic gain or loss;
  small elevations on GC-rich chromosomes are consistent with GC bias.
- AlphaGenome regulatory prediction for the missense variant (non-commercial terms):
  predicted BUB1B RNA expression essentially unchanged (ALT/REF 0.9994).

## Validation notes
- The variant-filtering script uses fixed ANNOVAR column numbers; verify them against
  your annotation header before interpreting results.
- The Snakemake workflow aligns FASTQ inputs and separately extracts variants from a VCF;
  confirm both belong to the same sample and reference build before phasing.
- The Manta launcher expects an already configured workflow directory.
- Exomiser configuration paths are placeholders and must match the installed version.
- Scripts with placeholder paths are templates, not validated production workflows.

## Limitations
Single case, no parental data (phase unconfirmed), one VUS allele, mosaicism screens have
limited sensitivity, repeat expansions not assessed.

## Reproducing
See `docs/REUSE.md`. Set paths in `config/paths.env`; adapt Slurm headers to your site.
Record exact tool and database versions, parameters and commands for every run, in a
private log that contains no patient information or credentials.

## Privacy and license
Do not add identifiable or sensitive clinical information. Released under CC BY 4.0
(see `LICENSE`, including third-party notes).
