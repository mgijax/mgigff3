#
# countPCrels.py
#
# Shows the counts of distinct pairs (p.type, c.type) where
# p is a feature and c is a child of p.
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
	roots[f.type] = roots.get(f.type,0) + 1
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
	    count(m, [m.type], m)
    #
    pcounts("root", roots)
    pcounts("mid", mids)
    pcounts("leaf", leaves)
    pcounts("path", paths)
    for k in exemplars:
	es = list(exemplars[k])
	es.sort()
	if len(es) > 5:
	    es = es[::len(es)/5]
        exemplars[k] = "," .join(es)
    pcounts("exemplars", exemplars)

#
ins = sys.argv[1:]if len(sys.argv) >= 2 else [ sys.stdin ]
main(ins)
