#
# mgiComputedMerge.py
#
# % python mgiComputedMerge.py PSLFILE MGISEQIDFILE > OUTFILE
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
	self.mgi2feats = {}

    def usage(self):
	print "USAGE: python %s PSLFILE MGISEQIDFILE > OUTFILE" % sys.argv[0]
	sys.exit(-1)

    def loadSeqidFile(self):
	fd = open(self.mgiFile, 'r')
	for line in fd:
	    seqid, mgiid, symbol, name, mcv_type, so_type = line.strip().split("\t")
	    self.seqid2gene[seqid] = mgiid
	    feat = self.mgi2feats.get(mgiid, None)
	    if feat is None:
		feat = gff3.Feature()
		feat.ID = mgiid
		feat.source = "MGI"
		feat.type = "gene"	# FIXME
		feat.currie = mgiid
		feat.Name = symbol
		feat.description = name
		feat.mcv_type = mcv_type
		feat.so_term_name = so_type
		self.mgi2feats[mgiid] = [ feat ]

    def loadPslFile(self):
	for f in psl.toGff(self.pslFile):
	    # 
	    f.seqid = f.seqid.replace("chr","")
	    s = f.qname.split(".")[0]
	    mgiid = self.seqid2gene[s]
	    f.mgi_id = mgiid
	    if f.type == "match": f.Parent = [ mgiid ]
	    #
	    self.mgi2feats[mgiid].append(f)
	    #
	    r = self.mgi2feats[mgiid][0]
	    r.seqid = f.seqid
	    r.strand = f.strand
	    r.start = f.start if r.start=="." else min(r.start, f.start)
	    r.end = f.end if r.end=="." else max(r.end, f.end)


    def main(self):
	def byTopLevel(a, b):
	    fa = a[0]
	    fb = b[0]
	    c = cmp(fa.seqid, fb.seqid)
	    if c : return c
	    c = cmp(fa.start, fb.start)
	    return c
	self.loadSeqidFile()
	self.loadPslFile()
	maker = gff3.IdMaker()
	allFeats = self.mgi2feats.values()
	allFeats.sort(byTopLevel)
	for feats in allFeats:
	    if len(feats) == 1:
	        continue
	    gff3.reassignIDs(feats, maker)
	    for f in feats:
	        sys.stdout.write(str(f))

MgiComputedMerger().main()
