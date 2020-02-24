#
# canonicalize.py
#
# Turns tree models in canonical form DAG-shaped) models by merging of identical subfeatures.
# 

import gff3
import sys

def mergeExon(ex, f):
    ex.Parent.extend(f.Parent)
    if f.source not in ex.source:
        ex.source += ("," + f.source)

for feats in gff3.models(sys.stdin, flatten=True):
    exons = {}   # (start,end) -> index into ofeats
    ofeats = []  # list of feats in current model. Root is 0th item.
    # merge exons. exons merge if they have the same coordinates.
    for f in feats:
        if f.type == "exon":
            k = (f.start,f.end)
            if k in exons:
                # merge f with a previously seen exon, ex
                # 
                i = exons[k]
                ex = ofeats[i]
                mergeExon(ex, f)
                # 
                # ...and move ex to end of the list
                ofeats[i] = None
                exons[k] = len(ofeats)
                ofeats.append(ex)
            else:
                # first time seeing this exon.
                exons[k] = len(ofeats)
                ofeats.append(gff3.Feature(f))
        else:
            # other feature types just append
            ofeats.append(gff3.Feature(f))

    # merge transcripts. transcripts merge if all their exons merged and their CDSs merged
    for f in feats:
        if f.type == "exon":
           pass

    # output
    for f in ofeats:
        if f:
            sys.stdout.write(str(f))

