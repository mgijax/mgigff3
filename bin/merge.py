#
# merge.py
#
# Merges n sorted gff files into one and merges models from other providers under the corresponding MGI feature.
#
#

import gff3
import sys

WINDOWSIZE = 200000

class ModelMerger:
    def __init__(self, wsize=WINDOWSIZE):
	self.idMaker = gff3.IdMaker()
	self.pendingMgi = []
	self.pendingNonMgi = []
	self.windowSize = wsize

    # Prints a model to standard out. 
    #
    def printModel(self, feats):
	for f in feats:
	    sys.stdout.write(str(f))

    def propagatePartIds(self, m):
	if not "gene_id" in m.attributes:
	    return
        gene_id = m.gene_id
	for t in m.children:
	    t.gene_id = gene_id
	    tid = t.transcript_id if "transcript_id" in t.attributes else None
	    for e in t.children:
	        e.gene_id = gene_id
		if tid:
		    e.transcript_id = tid

    # Merges model f into m. f is the root of a gene model from some provider. 
    # m is the MGI feature and is the root of a merged model hierarchy that
    # grows as models from multiple providers are merged.
    #
    # The model hierarchy under f is copied into m, rather than moved, because
    # of complex (non 1-1) associations among models. (I.e., f may be merged
    # into more than one m).
    #
    def mergeModels(self, m, f):
	m.start = min(m.start, f.start)
	m.end   = max(m.end,   f.end)
	feats = gff3.copyModel(gff3.flattenModel(f))
	#
	f._merged.append(m.curie)
	m._merged.append(f.curie)
	for c in feats[0].children:
	    c.Parent=[m.ID]
	    c.parents.clear()
	    c.parents.add(m)
	    m.children.add(c)

    # Message logger
    #
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
	        if xr and xr not in m._merged:
		    self.log("Dangling reference (%s) in %s\n" % (xr, m.curie))
	    flushed.append(m)
	    m.attributes.pop('_merged', None)
	    self.pendingMgi.pop(0)

	while len(self.pendingNonMgi):
	    m = self.pendingNonMgi[0]
	    if f and f.start - m.end < self.windowSize:
		break
	    if len(m.attributes.get("_merged", [])) == 0:
		self.log("Orphan model: ID=%s curie=%s type=%s\n" % (m.ID, m.curie, m.type))
	    m.attributes.pop('_merged', None)
	    self.pendingNonMgi.pop(0)

	for m in flushed:
	    self.propagatePartIds(m)
	    gff3.reassignIDs(gff3.flattenModel(m), self.idMaker)

	return flushed

    # Adds MGI feature m to the pending queue. Merges in any nonMgi models from the
    # other queue that are referenced in Dbxrefs.
    #
    def addMgi(self, m):
	m._merged = []
	self.pendingMgi.append(m)
	for f in self.pendingNonMgi:
	    if f.curie in m.attributes.get("Dbxref",[]):
		self.mergeModels(m, f)

    # Adds nonMgi model f to the pending queue. Merges f into any MGI feature from the
    # other queue that references f in its Dbxrefs.
    #
    def addNonMgi(self, f):
	f._merged = []
	self.pendingNonMgi.append(f)
	for m in self.pendingMgi:
	    if f.curie in m.attributes.get("Dbxref",[]):
		self.mergeModels(m, f)

    # Last check
    def validate(self, m):
	if len(m.children) == 0:
	    self.log("Gene model is not 3 levels: culling:")
	    self.log(str(m))
	    return None
        return m

    # Merges the sorted gff files listed on the command line into a single stream,
    # then merges features based on Dbxrefs. Maintains a kind of moving window in the
    # form of two "pending queues". As items are consumed from the sorted stream, they are
    # appended to the queues. Items are removed from the front of the queue when the 
    # current position (start coordinate of the latest input feature) moves beyond 
    # their endpoints. 
    #
    def merge(self, files):
	iters = [gff3.models(x) for x in files]
	for m in gff3.merge(*iters):
	    for mm in self.flush(m):
		if self.validate(mm):
		    yield mm
	    if m.source == "MGI":
		self.addMgi(m)
	    else:
		self.addNonMgi(m)

	for mm in self.flush():
	    if self.validate(mm):
	        yield mm

    # end class ModelMerger

#####

if __name__ == "__main__":
    for m in ModelMerger().merge(sys.argv[1:]):
        for f in gff3.flattenModel(m):
	    sys.stdout.write(str(f))
	sys.stdout.write("###\n")
