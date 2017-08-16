
import sys
import gff3
from OrderedSet import OrderedSet
from itertools import ifilter

EXCLUDE_TYPES = OrderedSet([
    "chromosome",
    "biological_region"
])

feats = ifilter(lambda f: f.type not in EXCLUDE_TYPES, gff3.iterate(sys.stdin))
for m in gff3.models(feats):
    for f in gff3.flattenModel(m):
	if not f.type in EXCLUDE_TYPES:
	    f.attributes.pop("Name", None)
	    f.attributes.pop("version", None)
	    f.attributes.pop("description", None)
	    f.attributes.pop("logic_name", None)
	    f.attributes.pop("gene_id", None)
	    if len(f.parents) == 0:
	        f.attributes["curie"] = "ENSEMBL:" + f.ID.split(":")[1]
	    sys.stdout.write(str(f))
