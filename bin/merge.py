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
	feats = gff3.flattenModel(f)
	for f2 in feats:
	    f2.attributes["mgi_id"] = m.curie
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
	    m.attributes.pop('_merged', None)
	    feats = gff3.flattenModel(m)
	    gff3.reassignIDs(feats, self.idMaker)
	    self.printModel(feats)

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

    # Merges the sorted gff files listed on the command line into a single stream,
    # then merges features based on Dbxrefs. Maintains a kind of moving window in the
    # form of two "pending queues". As items are consumed from the sorted stream, they are
    # appended to the queues. Items are removed from the front of the queue when the 
    # current position (start coordinate of the latest input feature) moves beyond 
    # their endpoints. 
    #
    def main(self, files):
	iters = [gff3.models(x) for x in files]
	for m in gff3.merge(*iters):
	    self.flush(m)
	    if m.source == "MGI":
		self.addMgi(m)
	    else:
		self.addNonMgi(m)

	self.flush()

    # end class ModelMerger

#####

if __name__ == "__main__":
    ModelMerger().main(sys.argv[1:])
