# MVA Hackathon — Track 1

This repository contains reusable scripts and supporting resources for
variant analysis and candidate prioritization for the MVA hackathon.

## Repository layout

- `analysis/submission/`: submission CSV utility script.
- `config/`: example analysis configuration templates.
- `docs/`: gene lists and exon BED resources.
- `scripts/annotation/`: ANNOVAR annotation workflow.
- `scripts/exomiser/`: Exomiser analysis launcher.
- `scripts/mosaic_screen/`: allele-balance screening utility.
- `scripts/snakemake/`: example alignment and phasing workflow.
- `scripts/sv_calling/`: Manta workflow launcher.
- `scripts/variant_filtering/`: annotation-based filtering and gene checks.

## Data and privacy

Patient-derived genomic data, FASTQ/BAM/CRAM/VCF files, generated
annotation outputs, credentials, and private configuration files must
not be committed to this repository.

Use synthetic or otherwise authorized test data when validating scripts.
Do not publish identifiable or confidential clinical information.

## Configuration requirements

Scripts currently contain placeholder paths and sample identifiers.
Replace these with paths and identifiers appropriate to the authorized
analysis environment before running a workflow.

Do not assume that all workflows are connected into one end-to-end
pipeline. Review each script's inputs, outputs, software requirements,
reference genome, and sample provenance before execution.

## Software and reference requirements

The workflows may require tools such as BWA, SAMtools, BCFtools,
WhatsHap, ANNOVAR, Exomiser, Manta, Snakemake, Python, Java, and
associated annotation databases.

Exact versions and database releases must be recorded for each actual
analysis run. Do not infer the versions used in a previous run from
the current shell environment.

## Important validation notes

- The variant-filtering script uses fixed ANNOVAR output column numbers.
  Verify these against the exact annotation output header and protocol
  configuration before interpreting results.
- The Snakemake workflow aligns FASTQ inputs and separately extracts
  variants from a VCF. Verify that these inputs correspond to the same
  sample and compatible reference build before phasing.
- The Manta launcher expects an already configured workflow directory.
  Confirm that the configuration and input BAMs exist before execution.
- The Exomiser configuration template contains placeholder paths and
  must be validated against the installed Exomiser version.
- The BAF utility should be tested with representative plain and
  compressed VCF files, including missing or incomplete genotype fields.

## Reproducibility

For every analysis, record the input identifiers, reference assembly,
software versions, database releases, parameters, commands, and output
locations in an appropriate private run log. Do not put sensitive
patient information or credentials in public logs.

## Status

This repository is being prepared for review. Scripts with placeholder
paths or unverified assumptions are templates and should not be treated
as validated production workflows.
