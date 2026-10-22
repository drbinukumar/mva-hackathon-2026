import os
from alphagenome.data import genome
from alphagenome.models import dna_client

model = dna_client.create(os.environ["ALPHAGENOME_API_KEY"])

chrom, pos, ref, alt = "chr15", 40220612, "T", "G"
GENE_START, GENE_END = 40161069, 40221122   # BUB1B span from the exon BED

start = pos - 2**19
interval = genome.Interval(chromosome=chrom, start=start, end=start + 2**20)
variant = genome.Variant(chromosome=chrom, position=pos,
                         reference_bases=ref, alternate_bases=alt)

out = model.predict_variant(
    interval=interval, variant=variant,
    ontology_terms=["UBERON:0001157"],   # quickstart placeholder tissue
    requested_outputs=[dna_client.OutputType.RNA_SEQ],
)

ref_t, alt_t = out.reference.rna_seq, out.alternate.rna_seq
print("track shape:", ref_t.values.shape, "interval:", ref_t.interval)
i0 = GENE_START - ref_t.interval.start
i1 = GENE_END - ref_t.interval.start
r, a = ref_t.values[i0:i1].sum(), alt_t.values[i0:i1].sum()
print("predicted BUB1B RNA-seq signal  REF=%.3f  ALT=%.3f  ALT/REF=%.4f" % (r, a, a / r))
