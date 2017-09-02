#
# mgiComputedMerge.py
# Merges the filtered blat results (PSLFILE) with the gene data from MGI (MGISEQIDFILE)
# and generates the MGI blatted model file.
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

    def log(self, m):
        sys.stderr.write(m)

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
	# Each PSL line is parsed and turned into a gff3 feature hierarchy (match->match_part*)
	# Here we iterate over the model roots (matches); the match_parse dangle below (f.children)
	for f in gff3.models(psl.toGff(self.pslFile)):
	    # 
	    seqid = f.qName.split(".")[0]	# get the seqid w/o version number
	    mgiid = self.seqid2gene[seqid]	# lookup the corresponding mgiid
	    mfeats = self.mgi2feats[mgiid]	# list containing the gene followed by its match features
	    mf = mfeats[0]			# the top-level mgi feature
	    #
	    # All the features n the model need tweaking:
	    for ff in gff3.flattenModel(f):
		# "chromosome: replace "chr5" for example with just "5"
		ff.seqid = ff.seqid.replace("chr","")
		# tag with the MGI#
		ff.mgi_id = mgiid

	    # Attach the match feature to the gene
	    f.Parent = [ mgiid ]

	    # If first match for the gene, set its chromosome, strand, start, end
	    if mf.seqid == ".":
		mf.seqid = f.seqid
		mf.strand = f.strand
		mf.start = f.start if mf.start=="." else min(mf.start, f.start)
		mf.end = f.end if mf.end=="." else max(mf.end, f.end)
	    else:
	        dupes = filter(lambda x: x.qName == f.qName, mfeats[1:])
		# Only accept sequences that have only one alignment
		if len(dupes) > 0:
		    f._rejected = True
		    if len(dupes) == 1:
			# First dupe discovered for this seqid. Log the error and invalidate the first guy
			self.log("REJECTING SEQUENCE - Seqid (%s) for gene (%s) has multiple alignments\n"%(seqid,mgiid))
			dupes[0]._rejected = True
		# If mult seqs align to multiple chroms, reject the lot.
		elif mf.seqid != f.seqid:
		    self.log("REJECTING GENE - Gene(%s) - Seqid (%s) matches on different chromosome than (%s)\n" \
		       %(mgiid, seqid, mfeats[1].qName))
		    mf._rejected = True
		# track the start/end range of subfeatures in the gene feature
		else:
		    mf.start = f.start if mf.start=="." else min(mf.start, f.start)
		    mf.end = f.end if mf.end=="." else max(mf.end, f.end)
	    # append match (even if rejected)
	    mfeats.append(f)


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
	    mf = feats[0]
	    matches = filter( lambda x: not "_rejected" in x.attributes, feats[1:] )
	    if "_rejected" in mf.attributes or len(matches) == 0:
	        continue
	    model = [mf]
	    for m in matches:
	        model += gff3.flattenModel(m)
	    gff3.reassignIDs(model, maker)
	    for f in model:
	        sys.stdout.write(str(f))

MgiComputedMerger().main()
