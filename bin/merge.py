#
# merge.py
#
# Merges models into the MGI gene.
# 
#	$ python merge.py file1 file2 [file3 ...]
#
# The files must contain data for a single chromosome (and all must be the same chromosome!)
# Example:
#	$ python merge.py ../working/mgi.chr14.gff ../working/ncbi.chr14.gff > test.gff
# Each file must be sorted appropriately. Model features clustered, no forward Parent references, etc,
#
# This script merges at two levels: 
#    1. merges n input files into one output file 
#    2. merges m models of a gene from m providers into a single model.
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

    # Message logger
    #
    def log(self, m):
        sys.stderr.write(m)

    # Prints a model to standard out. 
    #
    def printModel(self, feats):
	for f in feats:
	    sys.stdout.write(str(f))

    # Propagates each gene_id to all features in the gene.
    # Propagates transcript_id to all features in the transcript.
    #
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

    #
    def reassignIDs(self, m):
	idMap = {}
	idCounts = {}
	feats = gff3.flattenModel(m)
        for f in feats:
	    # handle CDSs, where multiple features have same ID.
	    if f.attributes.get("ID", None) in idMap:
	        f.ID = idMap[f.ID]
		continue
	    # generate ID
	    tp = f.type
	    tp = tp.replace("pseudogenic_","")
	    if tp.endswith("RNA") or "transcript" in tp:
	        tp = "transcript"
	    tp = tp[0].lower()
	    count = idCounts[tp] = idCounts.setdefault(tp, 0) + 1
	    newid = m.curie if f is m else "%s.%s%d" % (m.curie,tp,count)
	    # update f, idMap
	    if "ID" in f.attributes: idMap[f.ID] = newid
	    f.ID = newid
	# change Parent refs accordingly
	for f in feats:
	    if "Parent" in f.attributes:
	        f.Parent = [ idMap[pid] for pid in f.Parent ]

    # Flushes the current pending queues (pendingMgi and pendingNonMgi)
    # based on the latest feature from the merged stream.
    # Returns the flushed items from the pendingMgi queue, i.e., the root
    # features of the flushed models.
    #
    # Flushing depends of the fact the file is (roughly) sorted by increasing
    # start position of the features. As each feature f is read from the input,
    # the cache is scanned to see what can be flushed. For a given feature g in
    # the cache, if g.end < f.start, then g is not overlapped by f and cannot
    # be overlapped by anything following f, and so g can be flushed.
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
	    self.reassignIDs(m)

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

    # Hook to do any last checking. 
    # Return m if valid, else None.
    #
    def validate(self, m):
	if len(m.children) == 0:
	    self.log("Warning: Gene has no subfeatures: " + str(m))
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
