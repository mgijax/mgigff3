#
# merge.py
#

import gff3
import sys

class ModelMerger:
    def __init__(self):
	self.idMaker = gff3.IdMaker()
	self.pendingMgi = []
	self.pendingNonMgi = []
	self.windowSize = 200000

    def printModel(self, feats):
	for f in feats:
	    sys.stdout.write(str(f))

    def reassignIDs(self, m):
        feats = gff3.flattenModel(m)
	idmap = {}
	for f in feats:
	    if 'ID' not in f.attributes:
	        f.ID = self.idMaker.next(f.type)
	    elif f.ID in idmap:
	        f.ID = idmap[f.ID]
	    else:
		i = self.idMaker.next(f.type)
		idmap[f.ID] = i
		f.ID = i
	for f in feats:
	    pds = [ p.ID for p in f.parents ]
	    f.Parent = ",".join(pds)

    def mergeModels(self, m, f):
	m.start = min(m.start, f.start)
	m.end   = max(m.end,   f.end)
	f._merged.append(m.curie)
	m._merged.append(f.curie)
	for c in f.children:
	    c.Parent=[m.ID]
	    c.parents = [m]
	    m.children.add(c)

    def log(self, m):
        sys.stderr.write(m)

    # Flushes the current pending queues (pendingMgi and pendingNonMgi)
    # based on the latest feature from the merged stream.
    # If the start position of f is greater then the end position of 
    # an earlier cached feature, m, then neither f nor anything that 
    # comes after f overlaps m.
    #
    def flush(self, f=None):
	flushed = []
	while len(self.pendingMgi):
	    m = self.pendingMgi[0]
	    if f and f.start - m.end < self.windowSize:
		break
	    for xr in m.attributes.get("Dbxref",[]):
	        if not xr in m._merged:
		    self.log("Dangling reference (%s) in %s\n" % (xr, m.curie))
	    flushed.append(m)
	    self.pendingMgi.pop(0)

	while len(self.pendingNonMgi):
	    m = self.pendingNonMgi[0]
	    if f and f.start - m.end < self.windowSize:
		break
	    if len(m.attributes.get("_merged", [])) == 0:
		self.log("Orphan model: ID=%s curie=%s type=%s\n" % (m.ID, m.curie, m.type))
	    self.pendingNonMgi.pop(0)

	for m in flushed:
	    self.reassignIDs(m)
	    self.printModel(gff3.flattenModel(m))

    def addMgi(self, m):
	m._merged = []
	self.pendingMgi.append(m)
	for f in self.pendingNonMgi:
	    if f.curie in m.attributes.get("Dbxref",[]):
		self.mergeModels(m, f)

    def addNonMgi(self, f):
	f._merged = []
	self.pendingNonMgi.append(f)
	for m in self.pendingMgi:
	    if f.curie in m.attributes.get("Dbxref",[]):
		self.mergeModels(m, f)

    def main(self):
	iters = [gff3.models(x) for x in sys.argv[1:]]
	for m in gff3.merge(*iters):
	    self.flush(m)
	    if m.source == "MGI":
		self.addMgi(m)
	    else:
		self.addNonMgi(m)

	self.flush()

ModelMerger().main()
