# 
# mcv2go.py
#
# Exports mcv2go a dict mapping MCV terms to SO terms and IDs
#
#
# FIXME. This mapping should come out of the database. For now hardcode it. FIXME.
mcv2soData = [
    # ['MCV term', 'SO term', 'SO id']
    ['antisense lncRNA gene', 'antisense_lncRNA_gene', 'SO:0002182'],
    ['bidirectional promoter lncRNA gene', 'bidirectional_promoter_lncRNA_gene', 'SO:0002185'],
    ['gene', 'gene', 'SO:0000704'],
    ['gene segment', 'gene_segment', 'SO:3000000'],
    ['heritable phenotypic marker', 'heritable_phenotypic_marker','SO:0001500'],
    ['lincRNA gene', 'lincRNA_gene', 'SO:0001641'],
    ['lncRNA gene', 'lncRNA_gene', 'SO:0002127'],
    ['miRNA gene', 'miRNA_gene', 'SO:0001265'],
    ['non-coding RNA gene', 'ncRNA_gene', 'SO:0001263'],
    ['polymorphic pseudogene', 'polymorphic_pseudogene', 'SO:0001841'],
    ['protein coding gene', 'protein_coding_gene', 'SO:0001217'],
    ['pseudogene', 'pseudogene', 'SO:0000336'],
    ['pseudogenic gene segment', 'pseudogenic_gene_segment', 'SO:0001741'],
    ['rRNA gene', 'rRNA_gene', 'SO:0001637'],
    ['ribozyme gene', 'ribozyme_gene', 'SO:0002181'],
    ['RNase MRP RNA gene', 'RNase_MRP_RNA_gene', 'SO:0001640'],
    ['RNase P RNA gene', 'RNase_P_RNA_gene', 'SO:0001639'],
    ['scRNA gene', 'scRNA_gene', 'SO:0001266'],
    ['sense intronic lncRNA gene', 'sense_intronic_lncRNA_gene', 'SO:0002184'],
    ['sense overlapping lncRNA gene', 'sense_overlap_lncRNA_gene', 'SO:0002183'],
    ['snRNA gene', 'snRNA_gene', 'SO:0001268'],
    ['snoRNA gene', 'snoRNA_gene', 'SO:0001267'],
    ['SRP RNA gene', 'SRP_RNA_gene', 'SO:0001269'],
    ['tRNA gene', 'tRNA_gene', 'SO:0001272'],
    ['telomerase RNA gene', 'telomerase_RNA_gene', 'SO:0001643'],
    ['unclassified gene', 'gene', 'SO:0000704'],
    ['unclassified non-coding RNA gene', 'ncRNA_gene', 'SO:0001263'],
]
mcv2so = {}
for r in mcv2soData:
  mcv2so[r[0]] = r[1]

