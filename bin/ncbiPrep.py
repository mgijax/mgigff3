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

TAB     = '\t'
HASH    = '#'

# Column 3 types to filter out.
# May need to update this set when there is a new NCBI data release.
#
EXCLUDE_TYPES = set([
  "centromere",
  "D_loop",
  "cDNA_match",
  "enhancer",
  "match",
  "origin_of_replication",
  "region",
  "sequence_alteration",
  "sequence_feature",
  "CAAT_signal",
  "DNAseI_hypersensitive_site",
  "TATA_box",
  "TSS",
  "biological_region",
  "conserved_region",
  "enhancer_blocking_element",
  "imprinting_control_region",
  "locus_control_region",
  "matrix_attachment_site",
  "mobile_genetic_element",
  "nucleotide_motif",
  "promoter",
  "protein_binding_site",
  "replication_regulatory_region",
  "response_element",
  "sequence_alteration",
  "sequence_feature",
  "silencer",
  "transcriptional_cis_regulatory_region",
])

class ConvertNCBI:
    def __init__(self):
        self.currentRegionId = None
        self.currentRegion = None
        self.transcriptCount = 0
        self.exonCount = 0
        
    def getGeneID(self, f):
        ids = f.attributes.get('Dbxref',[])
        if type(ids) is str:
            ids = [ids]
        for x in ids:
            if x.startswith('GeneID:'):
                return x[7:]
        #raise RuntimeError("No GeneID found in " + str(f))
        return None


    # Processes the next feature from the NCBI gff3 file.
    # Returns the feature, or None if the feature should be skipped.
    # (One of the main things this code does is decide which features 
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
        elif f[2] == "lnc_RNA":
          f[2] = "lncRNA"
        #
        xrs = f.attributes.get('Dbxref', [])
        if type(xrs) is str:
            xrs = [xrs]

        #
        f[0] = self.currentRegion
        f[1] = "NCBI"
        if "Parent" not in f.attributes:
            # set the curie for top level nodes
            xrs = list(filter(lambda x: x.startswith("GeneID:"), xrs))
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
            tps = list(skipped.keys())
            tps.sort()
            for t in tps:
                self.log("%s\t%d\n"%(t, skipped[t]))

    def log(self, s):
        sys.stderr.write(s)

    # NCBI provides 1- and 2-level pseudogenes. We need 3-levels. 
    # Checks the provided model to see if it's a pseudogene with no or
    # only exon children. For 2 level models, (pseudogene->exon) insert
    # pseudogenic transcript. For 1 level model, adds pseudogenic
    # transcript and exon. Always results in 3 level model.
    #
    def checkPseudogene(self, m):
        if not m.type == "pseudogene":
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
        if len(m.children) > 0:
            # convert 2 level model by inserting a transcript
            for c in m.children:
                c.Parent = [ t.ID ]
                c.type = "pseudogenic_exon"
            gff3.crossReference([m, t] + list(m.children))
        else:
            # convert a 1-level model by appending a transcript and an exon
            self.exonCount += 1
            e = gff3.Feature(t)
            e.ID = "inserted_exon_%d" % self.exonCount
            e.type = "pseudogenic_exon"
            e.Parent = [t.ID]
            gff3.crossReference([m, t, e])

        
    #
    def checkTranscriptNames(self, m):
        curie = m.curie
        for f in gff3.flattenModel(m):
            if len(f.parents) and len(f.children) and "Name" not in f.attributes:
                try:
                    #f.Name = list(f.parents)[0].curie + "_" + f.type
                    f.Name = curie + "_" + f.type
                except:
                    print("ERROR:", str(f))
                    print(str(gff3.flattenModel(m)))
                    sys.exit(-1)

    # Converts NCBI structure into that required for the Alliance
    # NCBI structure for miRNA genes:
    #   gene
    #      primary_transcript
    #          exon
    #          miRNA
    #              exon
    #          miRNA
    #              exon
    # into Alliance specificied format. Note change in
    # parentage of the miRNAs, and the addition of Derives_from
    # relationships
    #   gene
    #      primary_transcript
    #          exon
    #      miRNA (Derives_from the primary transcr)
    #          exon
    #      miRNA (Derives_from the primary transcr)
    #          exon
    def checkMiRnas(self, m):
        if m.attributes.get('so_term_name','') != 'miRNA':
            return
        for f in gff3.flattenModel(m):
            if f.type == 'miRNA':
                # move miRNA to be direct child of gene
                p = list(f.parents)[0]
                pid = p.ID
                p.children.remove(f)
                #
                m.children.add(f)
                f.parents.clear()
                f.parents.add(m)
                f.Parent = [ m.ID ]
                f.Derives_from = [ pid ]
            elif f.type == 'primary_transcript':
                f.type = 'pre_miRNA'
    #
    def main(self):
        for m in gff3.models(self.pre(sys.stdin)):
           self.checkMiRnas(m)
           self.checkPseudogene(m)
           self.checkTranscriptNames(m)
           for f in gff3.flattenModel(m):
               print(str(f), end='')

#
if __name__ == "__main__":
    ConvertNCBI().main()
