#
# convertNCBI.py
#
# Filters/tweaks the NCBI file.
# * converts column 1 values from sequence ids to chromosome (or contig) ids.
# * filters out unwanted feature types
# * converts features of type gene where pseudo=true, into features of type pseudogene.
#   Also converts their transcripts and exon to their pseudogenic cousins.
# * removes "unwanted" attributes
#
# Usage:
#  python convertNCBI.py < /path/to/ncbidatafile.gff3 > output.gff3
# 
# 
# 

import sys
import os
import types
import gff3
import optparse

TAB	= '\t'
HASH	= '#'
EXCLUDE_TYPES = set([
  "D_loop",
  "cDNA_match",
  "enhancer",
  "match",
  "origin_of_replication",
  "region",
  "sequence_alteration",
  "sequence_feature",
])

class ConvertNCBI:
    def __init__(self):
	self.currentRegionId = None
	self.currentRegion = None
	self.transcriptCount = 0
	
    def getGeneID(self, f):
	ids = f.attributes.get('Dbxref',[])
	if type(ids) is types.StringType:
	    ids = [ids]
	for x in ids:
	    if x.startswith('GeneID:'):
		return x[7:]
	#raise RuntimeError("No GeneID found in " + str(f))
	return None


    # Processes the next feature from the NCBI gff3 file.
    # Returns the feature, or None if the feature should be skipped.
    # (One of the main things this code does is decided which features 
    # we're interested in.)
    #
    # The NCBI file uses sequence ids in col 1 to refer to chromosomes (and contigs).
    # (E.g. "NC_000067.6").
    # In addition, the features on a given chromosome are preceded in the file by a feature 
    # of type chromosome which contains the chromosome name (e.g. "1")
    #
    def process(self, f):
	# NCBI file has multi-level sort, with first level being by region 
	# (e.g. a chromosome or a contig)
	# A feature with 'region' in col 3 introduces a new region, and its features
	# will follow. See if it's one we want.
	if f[2] == 'region':
	  if f[0] == self.currentRegionId:
	    # A region feature within the current region. Skip.
	    return None
	  # Region with a different ID. See if it's one we care about.
	  self.currentRegionId = f[0]
	  chr = f.attributes.get('chromosome',None)
	  map = f.attributes.get('map', None)
	  genome = f.attributes.get('genome', None)
	  if f.attributes.get('strain',None) != 'C57BL/6J':
	    # skip anything not on B6
	    self.currentRegion = None
	  elif chr == 'Unknown':
	    # unplaced contig
	    self.currentRegion = f[0]
	  elif genome == 'genomic':
	    # unlocalized contig
	    self.currentRegion = chr + '|' + f[0]
	  elif genome == 'chromosome':
	    # regular ol' chromosome
	    self.currentRegion = chr
	  elif genome == 'mitochondrion':
	    # mitochondrion
	    self.currentRegion = 'MT'
	  else:
	    # something else
	    self.currentRegion = None
	  return None
	#
	# A feature in the current region
	if f[0] != self.currentRegionId:
	  raise RuntimeError('Internal error. Region id mismatch detected.' + str(f))
	if self.currentRegion is None:
	  # don't care about this region
	  return None
	#  
	if f[2] in EXCLUDE_TYPES:
	  return None
	#
	f[0] = self.currentRegion
	f[1] = "NCBI"
	if "Parent" not in f.attributes:
	    # set the curie for top level nodes
	    xrs = f.attributes.get('Dbxref', [])
	    if type(xrs) is types.StringType:
		xrs = [xrs]
	    xrs = filter(lambda x: x.startswith("GeneID:"), xrs)
	    if len(xrs) == 1:
		f.attributes["curie"] = xrs[0].replace("GeneID:","NCBI_Gene:")
	    # set the so_term_name for top level nodes
	    biotype = f.attributes.get("gene_biotype", None)
	    if biotype:
		if biotype == "protein_coding":
		    biotype += "_gene"
		elif biotype == "pseudogene":
		    f[2] = "pseudogene"
		f.attributes["so_term_name"] = biotype
	#
	f.attributes.pop('gene_biotype',None)
	f.attributes.pop('Dbxref',None)
	f.attributes.pop('gene_synonym',None)
	f.attributes.pop('product',None)
	f.attributes.pop('model_evidence',None)
	f.attributes.pop('gbkey',None)
	f.attributes.pop('gene',None)
	f.attributes.pop('description',None)
	f.attributes.pop('pseudo',None)
	if f.type == "exon":
	    f.attributes.pop("transcript_id",None)
	    f.attributes.pop("ncrna_class",None)
	return f

    def pre(self, inp):
	skipped = {}
        for f in gff3.iterate(inp):
	    if self.process(f):
	       yield f
	    else:
	       n = skipped.get(f.type, 0)
	       skipped[f.type] = n + 1
	if len(skipped):
	    self.log("Counts of skipped features, by type:\n")
	    tps = skipped.keys()
	    tps.sort()
	    for t in tps:
		self.log("%s\t%d\n"%(t, skipped[t]))

    def log(self, s):
        sys.stderr.write(s)

    # NCBI provides 2-level pseudogenes. We need 3-levels. 
    # Checks the provided model to see if it's a pseudogene with
    # only exon children. If so, it inserts a pseudogenic_transcript
    # in between.
    #
    def checkPseudogene(self, m):
        if not m.type == "pseudogene" or len(m.children) == 0:
	    return
	for c in m.children:
	    if c.type != "exon":
	        return
	# If here, its a pseudogene with only exons for children.
	# Insert the transcript
	self.transcriptCount += 1
        t = gff3.Feature(m)
	t.attributes = {}
	t.ID = "inserted_transcript_%d" % self.transcriptCount
	t.type = "pseudogenic_transcript"
	t.Parent = [ m.ID ]
	for c in m.children:
	    c.Parent = [ t.ID ]
	    c.type = "pseudogenic_exon"
	gff3.crossReference([m, t] + list(m.children))
        
    def main(self):
	for m in gff3.models(self.pre(sys.stdin)):
	   self.checkPseudogene(m)
	   for f in gff3.flattenModel(m):
	       print str(f),

#
if __name__ == "__main__":
    ConvertNCBI().main()
