
import sys
import gff3
from OrderedSet import OrderedSet
from itertools import ifilter

EXCLUDE_SOURCES = OrderedSet([
    "mirbase"
])
EXCLUDE_TYPES = OrderedSet([
    "chromosome",
    "biological_region",
    "supercontig",
    "three_prime_UTR",
    "five_prime_UTR"
])

filtFcn = lambda f: f.source not in EXCLUDE_SOURCES and f.type not in EXCLUDE_TYPES
feats = ifilter(filtFcn, gff3.iterate(sys.stdin))
for m in gff3.models(feats):
    for f in gff3.flattenModel(m):
	if not f.type in EXCLUDE_TYPES and not f.source in EXCLUDE_SOURCES:
	    f.source = "ENSEMBL"
	    if len(f.parents) == 0:
	        f.attributes["curie"] = "ENSEMBL:" + f.ID.split(":")[1]
	    biotype = f.attributes.get("biotype", None)
	    if biotype and len(f.parents) == 0:
		if biotype == "protein_coding":
		    biotype += "_gene"
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
	    sys.stdout.write(str(f))
