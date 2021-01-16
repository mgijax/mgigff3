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
from lib import mgiadhoc as db

MIN_PCT_LENGTH = float(os.environ['BLAT_MIN_LENGTH'])
BUILD=os.environ['MOUSE_ASSEMBLY']
MAP_COLLECTION_NAME=os.environ['BLAT_C4AM_COL6']
MAP_COLLECTION_ABBREV=os.environ['BLAT_C4AM_COL7']

DOT = '.'
SPACE = ' '
COMMA = ','
TAB = '\t'
NL = '\n'

class MgiComputedMerger:
    def __init__(self):
        if len(sys.argv) != 5:
            self.usage()
        self.pslFile = sys.argv[1]
        self.mgiFile = sys.argv[2]
        self.outGffFile = sys.argv[3]
        self.outC4amFile = sys.argv[4]
        self.outDiffFile = self.outC4amFile.rsplit(DOT,1)[0] + ".diff.txt"
        self.seqid2gene = {}
        self.gene2seqids = {}
        self.seqid2type = {}
        self.seqid2div = {}
        self.mgi2feats = {}
        self.mgi2current = {}
        self.c4amOutBuffer = {}
        self.diffOutBuffer = {}

    def log(self, m):
        sys.stderr.write(m)

    def usage(self):
        print("USAGE: python %s PSLFILE MGISEQIDFILE OUTFILE C4AM_OUTFILE" % sys.argv[0])
        sys.exit(-1)

    def loadSeqidFile(self):
        fd = open(self.mgiFile, 'r')
        for line in fd:
            seqid, seq_type, division, mgiid, symbol, name, mgi_type, so_type, mgi_chr = line.strip().split(TAB)
            self.seqid2gene[seqid] = mgiid
            self.gene2seqids.setdefault(mgiid, []).append(seqid)
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
            seqid = m.qName.split(DOT)[0]       # get the seqid w/o version number
            mgiid = self.seqid2gene[seqid]      # lookup the corresponding mgiid
            mfeats = self.mgi2feats[mgiid]      # list containing the gene followed by its match features
            if m.pctLength < MIN_PCT_LENGTH:
                self.logRejects("REJECTING SEQUENCE (%s) for GENE (%s) - pctLength (%1.2f) less than minimum (%1.2f)" % (seqid,mgiid,m.pctLength,MIN_PCT_LENGTH))
                self.logRejects(str(m))
                continue
            mfeats.append(m)
            self.counts[seqid] = self.counts.setdefault(seqid, 0) + 1

    def loadCurrentModels (self) :
        q = '''
            SELECT a.accid as "mgiid", m.symbol, mc.chromosome as "chr", f.startcoordinate as "start", f.endcoordinate as "end", f.strand
            FROM
               map_coord_collection cc,
               map_coordinate c,
               map_coord_feature f,
               mrk_chromosome mc,
               mrk_marker m,
               acc_accession a
            WHERE cc.name = '%s'
            AND cc._collection_key = c._collection_key
            AND c._object_key = mc._chromosome_key
            AND c._map_key = f._map_key
            AND f._object_key = m._marker_key
            AND m._marker_key = a._object_key
            AND a._mgitype_key = 2
            AND a._logicaldb_key = 1
            AND a.preferred = 1
        ''' % MAP_COLLECTION_NAME
        for r in db.sql(q) :
            self.mgi2current[r['mgiid']] = (r['mgiid'], r['symbol'], [], r['chr'], int(r['start']), int(r['end']), r['strand'])

    def qNameNoVersion (self, f) :
        return f.qName.split(DOT)[0] 

    def logRejects (self, msg, rejects = None, sep=NL) :
        self.log(msg)
        if rejects and len(rejects) > 0:
            self.log(sep.join([str(s) for s in rejects]))
        self.log(NL)
        
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
                mf.strand = DOT
                for ff in estSingles:
                  for fff in gff3.flattenModel(ff):
                      fff.strand = DOT
                mfeats[1:] = estSingles
            else:
                self.log("GENE (%s) - only has DNA matches. Setting strand to '.'\n" % mgiid)
                mf.strand = DOT
                for ff in dnaSingles:
                  for fff in gff3.flattenModel(ff):
                      fff.strand = DOT
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
        self.gffofd = open(self.outGffFile, 'w')
        self.c4amOfd = open(self.outC4amFile, 'w')
        self.c4amOfd.write('build=%s\n' % BUILD)
        #
        self.openDiffLog()
        # 
        allFeats = self.mgi2feats.values()
        allFeats = list(filter(lambda x: "_rejected" not in x[0][8], allFeats))
        allFeats.sort(key=topLevelKey)
        #
        self.writtenIds = set()
        for feats in allFeats:
            mf = feats[0]
            matches = feats[1:]
            #
            seqids = list(map(lambda m: m.Name, matches))
            attrs = mf.attributes
            if "Parent" in attrs and len(attrs["Parent"]) > 0:
                # f is not a top-level feature. Skip it.
                continue
            self.writeC4am(mf, seqids)
            self.writeGff3(mf, matches)
        # end for loop
        #
        self.closeDiffLog()
        self.closeC4amFile()

    def writeGff3 (self, mf, matches) :
        model = [mf]
        for m in matches:
            model += gff3.flattenModel(m)
        for f in model:
            self.gffofd.write(str(f))
        self.gffofd.write("###\n")

    def writeC4am (self, mf, seqids) :
        self.writtenIds.add(mf.curie)
        ###
        newCoords = (mf.curie, mf.Name, seqids, mf.seqid, mf.start, mf.end, None if mf.strand == DOT else mf.strand)
        currCoords = self.mgi2current.get(mf.curie, None)
        action = "Unchanged"
        if currCoords is None or newCoords is None or newCoords[3:] != currCoords[3:]:
            action = "Added" if not currCoords else "Deleted" if not newCoords else "Changed"
            self.writeDiff(action, currCoords, newCoords)
        ###
        assigned = mf.attributes.pop("_assigned", False)
        c4amRec = [
            mf.curie,
            mf.seqid,
            str(mf.start),
            str(mf.end),
            "" if mf.strand == DOT else mf.strand,
            MAP_COLLECTION_NAME,
            MAP_COLLECTION_ABBREV,
            "",
            action,
            "UN" if assigned else mf.seqid,
            mf.Name,
            ','.join(seqids)
        ]
        self.c4amOutBuffer.setdefault(action,[]).append(c4amRec)

    def closeC4amFile(self):
        ks = list(self.c4amOutBuffer.keys())
        ks.sort()
        for action in ks:
            for r in self.c4amOutBuffer[action]:
                self.c4amOfd.write(TAB.join(r)+NL)

    def openDiffLog (self):
        self.diffofd = open(self.outDiffFile, 'w')
        self.noDiffs = True
        self.diffofd.write(TAB.join([
            "action",
            "MGI_id",
            "symbol",

            "prod_seqids",
            "prod_chr",
            "prod_start",
            "prod_end",
            "prod_strand",

            "blat_seqids",
            "blat_chr",
            "blat_start",
            "blat_end",
            "blat_strand",
            ]) + NL)

    def writeDiff (self, action, currCoords, newCoords) :
        self.noDiffs = False
        #
        mgiid  = (newCoords or currCoords)[0]
        symbol = (newCoords or currCoords)[1]
        currentCoords = currCoords[3:] if currCoords else 4*[DOT]
        newSeqids = SPACE.join(newCoords[2]) if newCoords else DOT
        allSeqids = SPACE.join(self.gene2seqids[mgiid]) if mgiid in self.gene2seqids else DOT
        newCoords = newCoords[3:] if newCoords else 4*[DOT]
        #
        self.diffOutBuffer.setdefault(action, []).append(TAB.join([
                action,
                mgiid,
                symbol,
                allSeqids,
                TAB.join(map(str, currentCoords)),
                newSeqids,
                TAB.join(map(str, newCoords)),
            ]) + NL)

    def closeDiffLog (self):
        # check for/report deletes.
        for mgiid, m in self.mgi2current.items():
            if mgiid not in self.writtenIds:
                self.writeDiff('Deleted', m, None)
        #
        ks = list(self.diffOutBuffer.keys())
        if len(ks) == 0:
            self.diffofd.write("No differences detected.\n")
        else:
            ks.sort()
            for action in ks:
                for line in self.diffOutBuffer[action]:
                    self.diffofd.write(line)
        self.diffofd.close()

    def main(self):
        self.loadCurrentModels()
        self.loadSeqidFile()
        self.loadPslFile()
        self.processAlignments()
        self.output()

MgiComputedMerger().main()
