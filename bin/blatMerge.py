#
# blatMerge.py
#
# Merges the filtered blat results (PSLFILE) with the gene data from MGI (MGISEQIDFILE)
# and generates the MGI blatted model file.
# It also generates a file suitable to use as input to the coordinate load process.
#
# Applies additional filters on the results.
#
# % python blatMerge.py PSLFILE MGISEQIDFILE > OUTFILE C4AM_OUTFILE
#
# The file is sorted in the standard way.
#

import sys
import os
import gff3
import psl

MIN_PCT_LENGTH = float(os.environ['BLAT_MIN_LENGTH'])
BUILD=os.environ['MOUSE_ASSEMBLY']
MAP_COLLECTION_NAME=os.environ['BLAT_C4AM_COL6']
MAP_COLLECTION_ABBREV=os.environ['BLAT_C4AM_COL7']

class MgiComputedMerger:
    def __init__(self):
        if len(sys.argv) != 5:
            self.usage()
        self.pslFile = sys.argv[1]
        self.mgiFile = sys.argv[2]
        self.outGffFile = sys.argv[3]
        self.outC4amFile = sys.argv[4]
        self.seqid2gene = {}
        self.seqid2type = {}
        self.seqid2div = {}
        self.mgi2feats = {}

    def log(self, m):
        sys.stderr.write(m)

    def usage(self):
        print("USAGE: python %s PSLFILE MGISEQIDFILE OUTFILE C4AM_OUTFILE" % sys.argv[0])
        sys.exit(-1)

    def loadSeqidFile(self):
        fd = open(self.mgiFile, 'r')
        for line in fd:
            seqid, seq_type, division, mgiid, symbol, name, mgi_type, so_type, mgi_chr = line.strip().split("\t")
            self.seqid2gene[seqid] = mgiid
            self.seqid2type[seqid] = seq_type
            self.seqid2div[seqid] = division
            feat = self.mgi2feats.get(mgiid, None)
            if feat is None:
                feat = gff3.Feature()
                feat.seqid = mgi_chr
                feat.ID = mgiid.replace('MGI:', 'MGI_C57BL6J_')
                feat.source = "MGI"
                feat.type = "pseudogene" if mgi_type == "pseudogene" else "gene"
                feat.curie = mgiid
                feat.gene_id = mgiid
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
            seqid = m.qName.split(".")[0]       # get the seqid w/o version number
            mgiid = self.seqid2gene[seqid]      # lookup the corresponding mgiid
            mfeats = self.mgi2feats[mgiid]      # list containing the gene followed by its match features
            if m.pctLength < MIN_PCT_LENGTH:
                self.logRejects("REJECTING SEQUENCE (%s) for GENE (%s) - pctLength (%1.2f) less than minimum (%1.2f)" % (seqid,mgiid,m.pctLength,MIN_PCT_LENGTH))
                self.logRejects(str(m))
                continue
            mfeats.append(m)
            self.counts[seqid] = self.counts.setdefault(seqid, 0) + 1

    def qNameNoVersion (self, f) :
        return f.qName.split(".")[0] 

    def logRejects (self, msg, rejects = None, sep='\n') :
        self.log(msg)
        if rejects and len(rejects) > 0:
            self.log(sep.join([str(s) for s in rejects]))
        self.log('\n')
        
    def processAlignments (self) :    
        for mgiid, mfeats in list(self.mgi2feats.items()):
            # the top-level mgi feature
            mf = mfeats[0]
            # throw out mult matches for the same seqid (ie only keep single matches)
            singles = [m for m in mfeats[1:] if self.counts[self.qNameNoVersion(m)] == 1]
            multiples = [m for m in mfeats[1:] if self.counts[self.qNameNoVersion(m)] > 1]
            if len(multiples) :
                self.logRejects("REJECTING SEQUENCES for GENE (%s) - multiple matches per sequence: %s" % (mgiid, set([ m.qName for m in multiples ])))
            # Remove single best hits to chromosome different from mgi genetic chromosome.
            ss = []
            for m in singles:
                # chromosome: replace "chr5" for example with just "5"
                m.seqid = m.seqid.replace("chr","")
                if m.seqid != mf.seqid and mf.seqid != "UN":
                    self.logRejects("REJECTING SEQUENCE (%s) for GENE (%s) - matches to different chromosome (%s) than MGI genetic chromosome (%s)" % (m.qName,mgiid,m.seqid,mf.seqid))
                else:
                    ss.append(m)
            singles = ss

            # if no single matches, reject the gene
            if len(singles) == 0:
                if len(multiples):
                    self.logRejects("REJECTING GENE (%s) - Only nonunique matches." % mgiid)
                else:
                    self.logRejects("REJECTING GENE (%s) - No matches." % mgiid)
                mf._rejected = True
                continue

            # if singles do not all agree on chromosome, reject the gene
            chroms = set([ s.seqid for s in singles ])
            if len(chroms) > 1:
                self.logRejects("REJECTING GENE (%s) - Sequences match to multiple chromosomes: %s " % (mgiid, set([ m.qName for m in singles ])))
                mf._rejected = True
                continue

            # If MGI chr is UN, assign the chromosome based on the blat hit
            if mf.seqid == "UN":
                mf.seqid = list(chroms)[0]
                mf._assigned = True
                self.log("ASSIGNING CHROMOSOME (%s) to gene (%s) - was UN.\n" % (mf.seqid, mgiid))

            # tweak the models
            for s in singles:
                for ss in gff3.flattenModel(s):
                    ss.seqid = ss.seqid.replace("chr","")
                    # tag with the MGI#
                    ss.mgi_id = mgiid
                # Attach the match feature to the gene
                s.Parent = [ mf.ID ]

            # Divide into DNA and RNA and EST matches
            rnaSingles = []
            estSingles = []
            dnaSingles = []
            for s in singles:
              seqid = self.qNameNoVersion(s)
              tp = self.seqid2type[seqid]
              dv = self.seqid2div[seqid]
              if tp == "RNA" and dv == "EST":
                estSingles.append(s)
              elif tp == "RNA" and dv != "EST":
                rnaSingles.append(s)
              else:
                dnaSingles.append(s)
              
            # If gene has RNA matches:
            #   - ignore DNA and EST matches
            #     - the MGI gff gene model constructed only from RNA match results
            #   - all RNA matches must be to same strand (else reject the lot)
            #   - gene's strand inferred from RNA matches
            # Else if gene has EST matches:
            #   - ignore DNA matches
            #     - the MGI gff gene model constructed only from EST match results.
            #   - EST match strand is considered unreliable
            #   - do not infer gene's strand - all strands should be set to "."
            # Else (gene has only DNA matches):
            #   - match strand is considered meaningless
            #   - do not infer strand for gene - all strands should be set to "."

            if len(rnaSingles) > 0 :
                strands = set([ s.strand for s in rnaSingles ])
                if len(strands) > 1 :
                    self.logRejects("REJECTING GENE (%s) - RNA sequences match to both strands: %s" % (mgiid, set([r.qName for r in rnaSingles])))
                    mf._rejected = True
                    continue
                mf.strand = list(strands)[0]
                mfeats[1:] = rnaSingles
            elif len(estSingles) > 0 :
                self.log("GENE (%s) - only has EST%s matches. Setting strand to '.'\n" % (mgiid, " and DNA" if len(dnaSingles) else ''))
                mf.strand = "."
                for ff in estSingles:
                  for fff in gff3.flattenModel(ff):
                      fff.strand = "."
                mfeats[1:] = estSingles
            else:
                self.log("GENE (%s) - only has DNA matches. Setting strand to '.'\n" % mgiid)
                mf.strand = "."
                for ff in dnaSingles:
                  for fff in gff3.flattenModel(ff):
                      fff.strand = "."
                mfeats[1:] = dnaSingles

            # Update the gene-level feature
            mf.seqid = list(chroms)[0]
            mf.start = min([ s.start for s in mfeats[1:] ])
            mf.end = max([ s.end for s in mfeats[1:] ])

    def output (self) :
        def topLevelKey (fs) :
            s = fs[0].seqid
            ss = fs[0].start
            if len(s) == 1 and s.isdigit():
                s = "0" + s
            return (s, ss)
        gffofd = open(self.outGffFile, 'w')
        c4amOfd = open(self.outC4amFile, 'w')
        c4amOfd.write('build=%s\n' % BUILD)
        
        allFeats = self.mgi2feats.values()
        allFeats = list(filter(lambda x: "_rejected" not in x[0][8], allFeats))
        allFeats.sort(key=topLevelKey)
        for feats in allFeats:
            mf = feats[0]
            matches = feats[1:]
            #
            seqids = list(map(lambda m: m.Name, matches))
            attrs = mf.attributes
            assigned = attrs.pop("_assigned", False)
            if "Parent" in attrs and len(attrs["Parent"]) > 0:
                # f is not a top-level feature. Skip it.
                continue
            ###
            c4amRec = [
                attrs["curie"],
                mf.seqid,
                str(mf.start),
                str(mf.end),
                "" if mf.strand == "." else mf.strand,
                MAP_COLLECTION_NAME,
                MAP_COLLECTION_ABBREV,
                "",
                "UN" if assigned else mf.seqid,
                mf.Name,
                ','.join(seqids)
            ]
            c4amOfd.write("\t".join(c4amRec)+"\n")
            ###
            model = [mf]
            for m in matches:
                model += gff3.flattenModel(m)
            for f in model:
                gffofd.write(str(f))
            gffofd.write("###\n")

    def main(self):
        self.loadSeqidFile()
        self.loadPslFile()
        self.processAlignments()
        self.output()

MgiComputedMerger().main()
