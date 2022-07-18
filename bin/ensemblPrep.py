
import sys
import gff3
from OrderedSet import OrderedSet


EXCLUDE_SOURCES = OrderedSet([
    "NCBI"
])

EXCLUDE_TYPES = OrderedSet([
    "chromosome",
    "biological_region",
    "supercontig",
    "three_prime_UTR",
    "five_prime_UTR"
])

filtFcn = lambda f: f.type not in EXCLUDE_TYPES and f.source not in EXCLUDE_SOURCES
feats = filter(filtFcn, gff3.iterate(sys.stdin))
for m in gff3.models(feats):
    for f in gff3.flattenModel(m):
        if f.attributes.get("ID","").startswith("transcript:"):
            f.Name = f.transcript_id
        f.source = "ENSEMBL"
        if len(f.parents) == 0:
            f.attributes["curie"] = "ENSEMBL:" + f.ID.split(":")[1]
        biotype = f.attributes.get("biotype", None)
        if biotype and len(f.parents) == 0:
            if biotype == "protein_coding":
                biotype = "protein_coding_gene"
            f.attributes["so_term_name"] = biotype
        f.attributes.pop("biotype", None)
        f.attributes.pop("version", None)
        f.attributes.pop("description", None)
        f.attributes.pop("logic_name", None)
        f.attributes.pop("gene_id", None)
        f.attributes.pop("transcript_support_level", None)
        f.attributes.pop("rank", None)
        f.attributes.pop("constitutive", None)
        f.attributes.pop("ensembl_end_phase", None)
        f.attributes.pop("ensembl_phase", None)
        if f.type == "CDS":
            f.Name = f.protein_id
        elif f.type == "lnc_RNA":
            f.type = "lncRNA"
        elif f.type == "miRNA" and len(f.children):
            # Ensembl miRNA models look ike this:
            #     gene -> miRNA -> exon
            # However, the miRNA actually has the coordinates of the immature transcript.
            # Here we change the type. Output models look like:
            #     gene -> pre_miRNA -> exon
            # FIXME: when Ensembl fixes their representations, revise this code.
            f.type = "pre_miRNA"
        sys.stdout.write(str(f))

