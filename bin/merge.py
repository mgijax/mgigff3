#
# merge.py
#
# Merges models into the MGI gene.
# 
#       $ python merge.py file1 file2 [file3 ...]
#
# The files must contain data for a single chromosome (and all must be the same chromosome!)
# Example:
#       $ python merge.py ../working/mgi.chr14.gff ../working/ncbi.chr14.gff > test.gff
# Each file must be sorted appropriately. Model features clustered, no forward Parent references, etc,
#
# This script merges at two levels: 
#    1. merges n input files into one output file 
#    2. merges m models of a gene from m providers into a single model.
#

from lib import gff3
import sys
import types

WINDOWSIZE = 200000

class ModelMerger:
    def __init__(self, wsize=WINDOWSIZE):
        self.idMaker = gff3.IdMaker()
        self.seqid = None
        self.window = []
        self.outputMgi = set()
        self.windowSize = wsize

    # Message logger
    #
    def log(self, m):
        sys.stderr.write(m)

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

    #
    def reassignIDs(self, m):
        self.idMap = {}
        self.idCounts = {}
        feats = gff3.flattenModel(m)
        for f in feats:
            # handle CDSs, where multiple features have same ID.
            if f.attributes.get("ID", None) in self.idMap:
                f.ID = self.idMap[f.ID]
                continue
            # generate ID
            tp = f.type
            tp = tp.replace("pseudogenic_","")
            if tp.endswith("RNA") or "transcript" in tp:
                tp = "transcript"
            tp = tp.lower()
            count = self.idCounts[tp] = self.idCounts.setdefault(tp, 0) + 1
            newid = m.ID if f is m else "%s_%s_%d" % (m.ID.replace(":","_"),tp,count)
            # update f, self.idMap
            if 'ID' in f.attributes: self.idMap[f.ID] = newid
            f.ID = newid
        # change Parent refs accordingly
        for f in feats:
            if "Parent" in f.attributes:
                f.Parent = [ self.idMap[pid] for pid in f.Parent ]
            if "Derives_from" in f.attributes:
                f.Derives_from = self.idMap[f.Derives_from]

    #
    def flush(self, f=None):
        flushed = []
        while len(self.window):
            m = self.window[0]
            if f and f.seqid == m.seqid and f.start - m.end < self.windowSize:
                break
            flushed.append(m)
            self.window.pop(0)

        for m in flushed:
            self.propagatePartIds(m)
            self.reassignIDs(m)
            self.outputMgi.add(m.curie)
            for xr in m.Dbxref:
                if xr not in m._merged:
                    self.log("Dangling xref %s in gene %s\n" % (xr, m.curie))
            m.attributes.pop("_merged", None)

        return flushed

    # Initializes a model for the current window. Injects MGI gene information into 
    # the root nodes of a provider model.
    def initModel (self, m, mgiGene):
        g = mgiGene
        #
        m = gff3.copyModel(gff3.flattenModel(m))[0]
        if m.seqid != mgiGene.seqid:
            self.log("Warning: chromosome mismatch: Gene=%s, MGI chr=%s, %s chr=%s\n" % (g.curie, g.seqid, m.source, m.seqid))
        #
        self.seqid = m.seqid

        # process type
        m.source = "MGI"
        m.type = mgiGene.type
        
        # save orig attributes
        ma = m.attributes
        # replace with mgiGene attributes
        m.attributes = {}
        m.attributes.update(mgiGene.attributes)
        # restore children
        # for each child, fix Parent attribute
        for mc in m.children:
            mc.Parent = [ m.ID ]
        # Keep track of provider IDs that have been merged into this gene
        m._merged = [ ma["curie"] ]
        return m

    # Merges model f into m. f is the root of a gene model from some provider. 
    # m is the MGI feature and is the root of a merged model hierarchy that
    # grows as models from multiple providers are merged. 
    # The root features of the providers' models are merged into the
    # mgi root feature, and the descendant trees from the provider model are grafted onto
    # the growing mgi tree.
    #
    # The model hierarchy under f is copied into m, rather than moved, because
    # of complex (non 1-1) associations among models. (I.e., f may be merged
    # into more than one m).
    #
    # Merging a gene with a pseudogene causes the addition of a biotypeConflict attribute.
    #
    def mergeModel(self, m, f):
        # sanity checks:
        if m.seqid != f.seqid or m.strand != f.strand:
            self.log("Cannot merge because chromosomes or strands disagree.\n%s\n%s\n" % (str(m), str(f)))
            return m
        m.start = min(m.start, f.start)
        m.end   = max(m.end,   f.end)
        feats = gff3.copyModel(gff3.flattenModel(f))
        #
        if ("pseudo" in m.type and "pseudo" not in f.type) \
        or ("pseudo" not in m.type and "pseudo" in f.type): 
            m.biotypeConflict = "true"
        #
        for c in feats[0].children:
            c.Parent=[m.ID]
            c.parents.clear()
            c.parents.add(m)
            m.children.add(c)
        #
        m._merged.append(f.curie)
        return m

    # processNext (nonMgi) model m.
    # One of 4 outcomes:
    #   1. m has no matching MGI gene, and gets kicks out and logged.
    #   2. The matching MGI gene has already been output. This gets kicked out and logged.
    #   3. The matching MGI gene is in the current window. Model m is merged.
    #   4. The matching MGI gene is not in the window and has not been seen. Initialize new model 
    #      in window for MGI gene + model m.
    #
    def processNext(self, m):
        gs  = self.id2mgiFeat.get(m.curie, None)
        if gs is None:
            self.log("Orphan model: %s \n" % str(m.curie))
            return
        for g in gs:
            if g.curie in self.outputMgi:
                self.log("MGI gene has already been output: " + str(g))
                self.log("Detected at: " + str(m))
                return
            for mm in self.window:
                if mm.curie == g.curie:
                    self.mergeModel(mm, m)
                    break
            else:
               self.window.append(self.initModel(m, g))

    #
    def loadMgiData (self, mgiFile) :
        self.mgiFeats = []
        self.id2mgiFeat = {}
        for f in gff3.iterate(mgiFile):
            self.mgiFeats.append(f)
            for xref in f[8].get('Dbxref',[]):
                if xref in self.id2mgiFeat:
                    self.log("WARNING: ID %s associated with multiple genes.\n" % xref)
                    self.id2mgiFeat[xref].append(f)
                else:
                    self.id2mgiFeat[xref] = [f]

    # Merges the sorted gff files listed on the command line into a single stream,
    # then merges features based on Dbxrefs. Maintains a kind of moving window in the
    # form of two "pending queues". As items are consumed from the sorted stream, they are
    # appended to the queues. Items are removed from the front of the queue when the 
    # current position (start coordinate of the latest input feature) moves beyond 
    # their endpoints. 
    #
    def merge(self, mgiFile, providerFiles):
        self.loadMgiData(mgiFile)
        iters = [gff3.models(x) for x in providerFiles]
        for m in gff3.merge(*iters):
            for mm in self.flush(m):
                yield mm
            self.processNext(m)
        #
        for mm in self.flush():
            yield mm
        #
        self.finalCheck()

    #
    def finalCheck (self) :
        chrFeats = filter(lambda f: f.seqid == self.seqid, self.mgiFeats)
        missingChrFeats = filter(lambda f: f.curie not in self.outputMgi, chrFeats)
        for f in missingChrFeats:
            self.log("Missing: " + str(f))
        
    # end class ModelMerger

#####

if __name__ == "__main__":
    merger = ModelMerger()
    for m in merger.merge(sys.argv[1], sys.argv[2:]):
        for f in gff3.flattenModel2(m):
            sys.stdout.write(str(f))
        sys.stdout.write("###\n")
