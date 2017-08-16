#
# mergeModels.py
#

import sys
import gff3
import itertools
import heapq

# global counter
N = 0
def nextCount():
    global N
    N += 1
    return N
    
# merges n inputs streams of features into one. 
# Args:
#  modelIterators list of iterators over models
# Returns:
#  iterator over the merged stream
def merge(*modelIters):
    global N
    mis = [ itertools.imap(lambda m:(m.seqid, m.start, nextCount(), m), i) for i in modelIters ]
    for m in heapq.merge(*mis):
        yield m[3]

def printModel(m):
    for f in gff3.flattenModel(m):
	sys.stdout.write(str(f))

# Merges a list of top-level features into one. Feats is a list of top-level features
# to be merged into a single top-level feature. Features 1-n are merged into feature 0.
# Merging means:
#  - children of feature i become children of feature 0
#  - the start/end positions of feature 0 are reset to the min start/max end positions of 
#    all the features 0..n
#  - attribute values are combined
def mergeFeatures(feats):
    m = feats[0]
    for f in feats[1:]:
        for c in f.children:
	    c.Parent=m.ID
	    c.parents = [m]
	    m.children.add(c)

def flush(q):
    ix = {}
    mgiFeats = []
    for f in q:
	try:
	    ix[f.curie] = f
	except:
	    print "#### NO curie id ####", f
	    continue
        if f.source == "MGI":
	    f.ID = f.curie
	    mgiFeats.append(f);
    for f in mgiFeats:
        mergeFeatures([f] + filter(None, [ ix.get(x, None) for x in f.attributes.get('Dbxref',[]) ]))
    q[:] = []
    for m in mgiFeats:
	for f in gff3.flattenModel(m):
	    sys.stdout.write(str(f))
	    
def main():
    q = []	# queue of features
    ix = {}	# index by curie id
    iters = [gff3.models(x) for x in sys.argv[1:]]
    for m in merge(*iters):
	if len(q) == 0 or  m.overlaps(q[0]):
	    q.append(m)
	else:
	    flush(q)
	    q.append(m)
    flush(q)

main()
