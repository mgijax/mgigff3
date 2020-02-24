#
# exome.py
#
# Post process mgi gff file to produce mgi exome file.
#

import sys
import gff3
import types

#
sys.stdout.write(gff3.HEADER)
for feats in gff3.models(sys.stdin, flatten=True):
    #
    m = feats[0]
    exons = [f for f in feats if "exon" in f.type or "match_part" == f.type] # includes pseudogenic exons
    #
    unique = {}
    for e in exons:
        k = (e.start, e.end)
        if k in unique:
            ee = unique[k]
            ee.provider.add(e.source)
        else:
            unique[k] = e
            if m.type == "pseudogene":
                e.type = "pseudogenic_exon"
            e.attributes = {
                "mgiSymbol": m.Name,
                "mgiId" : m.curie,
                "provider" : set([e.source])
            }
            e.source = "MGI"
    #
    ks = list(unique.keys())
    ks.sort()
    for i,k in enumerate(ks):
        e = unique[k]
        e.ID = "%s_%03d"%(m.curie, i+1)
        e.provider = list(e.provider)
        sys.stdout.write(str(e))
    # end for
#end for
