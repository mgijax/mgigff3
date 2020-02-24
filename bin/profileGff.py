#
# profileGff.py
#
# Generates a report of feature types and parent/child relationships
# in a gff3 file.
#

import gff3
import sys

roots = {}
leaves = {}
mids = {}
paths = {}
exemplars = {}

def count(f, path, root):
    if len(f.parents) == 0:
        stn = f.attributes.get("so_term_name",None)
        tp = f.type + ("[%s]"%stn if stn else "")
        roots[tp] = roots.get(tp,0) + 1
    if len(f.parents) and len(f.children):
        mids[f.type] = mids.get(f.type,0) + 1
    if len(f.children) == 0:
        leaves[f.type] = leaves.get(f.type,0) + 1
        pth = '|'.join(path)
        paths[pth] = paths.get(pth, 0) + 1

        if "ID" in root.attributes:
            exemplars.setdefault(pth,set()).add(root.ID)
    else:
        for c in f.children:
            count(c, path+[c.type], root)

def pcounts(msg, counts):
    ks = counts.keys()
    ks.sort()
    for k in ks:
        print "\t".join([msg, k, str(counts[k]) ])

 
def main(featureSources):
    for fSource in featureSources:
        for m in gff3.models(fSource):
            stn = m.attributes.get("so_term_name",None)
            count(m, [m.type+("[%s]"%stn if stn else "")], m)
    #
    pcounts("root", roots)
    pcounts("mid", mids)
    pcounts("leaf", leaves)
    for k in exemplars:
        es = list(exemplars[k])
        es.sort()
        if len(es) > 5:
            es = es[::len(es)/5]
        paths[k] = str(paths[k])+ "\t" + "," .join(es)
    pcounts("path", paths)

#
ins = sys.argv[1:]if len(sys.argv) >= 2 else [ sys.stdin ]
main(ins)
