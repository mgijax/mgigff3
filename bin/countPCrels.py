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
def count(f, path):
    if len(f.parents) == 0:
	roots[f.type] = roots.get(f.type,0) + 1
    if len(f.parents) and len(f.children):
        mids[f.type] = mids.get(f.type,0) + 1
    if len(f.children) == 0:
	leaves[f.type] = leaves.get(f.type,0) + 1
	pth = ' -> '.join(path)
	paths[pth] = paths.get(pth, 0) + 1
    else:
	for c in f.children:
	    count(c, path+[c.type])

def pcounts(msg, counts):
    print msg
    ks = counts.keys()
    ks.sort()
    for k in ks:
	print "\t"+k, counts[k]
 
def main():
    for m in gff3.models(sys.stdin, flatten=True):
	for f in m:
	    if len(f.parents) == 0:
		count(f, [f.type])
    pcounts("Freq distr of root node types:", roots)
    pcounts("Freq distr of midlevel node types:", mids)
    pcounts("Freq distr of leaf node types:", leaves)
    pcounts("Freq distr of root-to-leaf paths:", paths)


main()
