#
# convertNCBI.py
#
# Filters NCBI file for features on one of the chromosomes in the chromosome file.
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
	
    def getGeneID(self, f):
	ids = f.attributes.get('Dbxref',[])
	if type(ids) is types.StringType:
	    ids = [ids]
	for x in ids:
	    if x.startswith('GeneID:'):
		return x[7:]
	#raise RuntimeError("No GeneID found in " + str(f))
	return None


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
	xrs = f.attributes.get('Dbxref', [])
	if type(xrs) is types.StringType:
	    xrs = [xrs]
	xrs = filter(lambda x: x.startswith("GeneID:"), xrs)
	if len(xrs) == 1:
	    f.attributes["curie"] = "NCBI_Gene:" + xrs[0].split(":")[1]
	biotype = f.attributes.get("gene_biotype", None)
	if biotype:
	    f.attributes["so_term_name"] = biotype
	    del f.attributes["gene_biotype"]
	f.attributes.pop('product',None)
	f.attributes.pop('model_evidence',None)
	f.attributes.pop('gbkey',None)
	f.attributes.pop('gene',None)
	f.attributes.pop('Dbxref',None)
	f.attributes.pop('description',None)
	return f

    def pre(self, inp):
        for f in gff3.iterate(inp):
	    if self.process(f):
	       yield f

    def main(self):
	for m in gff3.models(self.pre(sys.stdin)):
	   for f in gff3.flattenModel(m):
	       print str(f),

#
if __name__ == "__main__":
    ConvertNCBI().main()
