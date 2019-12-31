#
# mgiComputedMerge.py
#
# Merges the filtered blat results (PSLFILE) with the gene data from MGI (MGISEQIDFILE)
# and generates the MGI blatted model file.
#
# % python mgiComputedMerge.py PSLFILE MGISEQIDFILE > OUTFILE
#
# The file is sorted in the standard way.
#

import sys
import gff3
import psl

class MgiComputedMerger:
    def __init__(self):
	if len(sys.argv) != 3:
	    self.usage()
	self.pslFile = sys.argv[1]
	self.mgiFile = sys.argv[2]
	self.seqid2gene = {}
	self.seqid2type = {}
	self.mgi2feats = {}

    def log(self, m):
        sys.stderr.write(m)

    def usage(self):
	print "USAGE: python %s PSLFILE MGISEQIDFILE > OUTFILE" % sys.argv[0]
	sys.exit(-1)

    def loadSeqidFile(self):
	fd = open(self.mgiFile, 'r')
	for line in fd:
	    seqid, mgiid, symbol, name, mgi_type, so_type, seq_type = line.strip().split("\t")
	    self.seqid2gene[seqid] = mgiid
	    self.seqid2type[seqid] = seq_type
	    feat = self.mgi2feats.get(mgiid, None)
	    if feat is None:
		feat = gff3.Feature()
		feat.ID = mgiid.replace('MGI:', 'MGI_C57BL6J_')
		feat.source = "MGI"
		feat.type = "gene"	# FIXME
		feat.curie = mgiid
		feat.Name = symbol
		feat.description = name
		feat.mgi_type = mgi_type
		feat.so_term_name = so_type
		self.mgi2feats[mgiid] = [ feat ]

    def loadPslFile(self):
	# Each PSL line is parsed and turned into a gff3 feature hierarchy (match->match_part*)
	# Here we iterate over the model roots (matches); the match_parts dangle below (f.children)
        self.counts = {}
	for m in gff3.models(psl.toGff(self.pslFile)):
	    # 
	    seqid = m.qName.split(".")[0]	# get the seqid w/o version number
	    mgiid = self.seqid2gene[seqid]	# lookup the corresponding mgiid
	    mfeats = self.mgi2feats[mgiid]	# list containing the gene followed by its match features
            mfeats.append(m)
            self.counts[seqid] = self.counts.setdefault(seqid, 0) + 1

    def qNameNoVersion (self, f) :
        return f.qName.split(".")[0] 

    def logRejects (self, msg, rejects = None) :
        self.log(msg)
        if rejects and len(rejects) > 0:
            self.log(''.join([str(s) for s in rejects]))
        self.log('\n')
        
    def processAlignments (self) :    
        for mgiid, mfeats in self.mgi2feats.items():
            # the top-level mgi feature
            mf = mfeats[0]
            # throw out mult matches for the same seqid (ie only keep single matches)
            singles = filter(lambda m: self.counts[self.qNameNoVersion(m)] == 1, mfeats[1:])
            multiples = filter(lambda m: self.counts[self.qNameNoVersion(m)] > 1, mfeats[1:])
            if len(multiples) :
                self.logRejects("REJECTING SEQUENCES due to multiple matches\n", multiples)
            # if no single matches, reject the gene
            if len(singles) == 0:
                self.logRejects("REJECTING GENE - Gene(%s) - No unique matches.\n" % mgiid)
                mf._rejected = True
                continue
            # tweak the models
            for s in singles:
                for ss in gff3.flattenModel(s):
		    # chromosome: replace "chr5" for example with just "5"
		    ss.seqid = ss.seqid.replace("chr","")
		    # tag with the MGI#
		    ss.mgi_id = mgiid
                # Attach the match feature to the gene
	        s.Parent = [ mf.ID ]
            # if singles do not all agree on chromosome, reject the lot
            chroms = set([ s.seqid for s in singles ])
            if len(chroms) > 1:
                self.logRejects("REJECTING GENE - Gene(%s) - Sequences match to multiple chromosomes\n" % mgiid, singles)
                mf._rejected = True
                continue
            # Divide into DNA and RNA matches
            rnaSingles = filter(lambda s: self.seqid2type[self.qNameNoVersion(s)] == "RNA", singles)
            dnaSingles = filter(lambda s: self.seqid2type[self.qNameNoVersion(s)] == "DNA", singles)
            # If gene has only DNA matches, then
            #   - all matches must be on same strand, else no model for gene
            #   - gene's strand inferred
            # Otherwise (gene has RNA matches)
            #   - all RNA matches must be on same strand, else no model
            #   - DNA matches on same strand are used. Those on other strand are kicked out/reported
            checkStrands = dnaSingles if len(rnaSingles) == 0 else rnaSingles
            strands = set([ s.strand for s in checkStrands ])
            if len(strands) > 1 :
                self.logRejects("REJECTING GENE - Gene(%s) - sequences match to both strands\n" % mgiid, checkStrands)
                mf._rejected = True
                continue
            # Update the gene-level feature
            mf.seqid = list(chroms)[0]
            mf.strand = list(strands)[0]

            if len(rnaSingles) == 0 :
              mfeats2 = dnaSingles
            else:
              mfeats2 = rnaSingles + filter(lambda s: s.strand == mf.strand, dnaSingles)
              rejects = filter(lambda s: s.strand != mf.strand, dnaSingles)
              if len(rejects) > 0 : 
                  self.logRejects("REJECTING DNA SEQUENCES for Gene(%s) - matching wrong strand.\n" % mgiid, rejects)
            
            mf.start = min([ s.start for s in mfeats2 ])
            mf.end = max([ s.end for s in mfeats2 ])
            #
            mfeats[1:] = mfeats2

    def output (self) :
	def byTopLevel(a, b):
	    fa = a[0]
	    fb = b[0]
	    c = cmp(fa.seqid, fb.seqid)
	    if c : return c
	    c = cmp(fa.start, fb.start)
	    return c
	allFeats = self.mgi2feats.values()
	allFeats.sort(byTopLevel)
	for feats in allFeats:
	    mf = feats[0]
	    matches = filter( lambda x: not "_rejected" in x.attributes, feats[1:] )
	    if "_rejected" in mf.attributes or len(matches) == 0:
	        continue
	    model = [mf]
	    for m in matches:
	        model += gff3.flattenModel(m)
	    for f in model:
	        sys.stdout.write(str(f))

    def main(self):
	self.loadSeqidFile()
	self.loadPslFile()
        self.processAlignments()
        self.output()

MgiComputedMerger().main()
