#
# prepMgi.py
#
# Generates a GFF3 file of genes and pseudogenes in MGI.
# Only gene/pseudogene features; no subfeatures (obviously).
# These genes/pseudogenes will for the roots of model hierarchies in the 
# output. 
# Each feature includes the MGI id, symbol, name, actual SO type (e.g., tRNA_gene) 
# and genomic coordinates.
# Each feature also includes Dbxrefs to their various gene models.
# These Dbxrefs will direct the merging of models downstream.
#
# IMPLEMENTATION: Uses MouseMine web services. Retrieves all SequenceFeatures 
# in the Mouse Genome Catalog (along with any model IDs)
# Filters for genes and pseudogenes.
# Writes GFF3 file.
#

import gff3
import sys
from intermine.webservice import Service
#
TAB	= "\t"
#
pnameMap = {	# map from provider name to prefixes.
    "NCBI Gene Model" : "NCBI_Gene",
    "Ensembl Gene Model" : "ENSEMBL",
    "miRBase" : "miRBase",
    }

def getFeatures():
    service = Service("http://www.mousemine.org/mousemine/service")
    #
    query = service.new_query("SequenceFeature")
    query.add_view(
	"primaryIdentifier", 
	"symbol", 
	"name", 
	"sequenceOntologyTerm.name",
	"mgiType", 
	"chromosomeLocation.locatedOn.primaryIdentifier",
	"chromosomeLocation.start",
	"chromosomeLocation.end",
	"chromosomeLocation.strand",
	"crossReferences.identifier", 
	"crossReferences.source.name"
    )
    query.add_sort_order("SequenceFeature.chromosomeLocation.locatedOn.primaryIdentifier", "ASC")
    query.add_sort_order("SequenceFeature.chromosomeLocation.start", "ASC")
    #
    query.add_constraint("organism.taxonId", "=", "10090", code = "B")
    query.add_constraint("dataSets.name", "=", "Mouse Gene Catalog from MGI", code = "C")
    query.add_constraint("crossReferences.source.name", "ONE OF", 
         ["Ensembl Gene Model", "NCBI Gene Model", "miRBase"], code = "D")
    query.outerjoin("crossReferences")
    #
    return query

#
def main():
    # FIXME: due to a bug in the Intermine client lib, the sort orders in the query are being ignored.
    # This forces us to accumulate results and sort internally before outputting anything.
    # Once the intermine bug is fixed, we should be able to output features in the order returned 
    # by the query.
    feats = []
    for f in getFeatures():
	soterm = f.sequenceOntologyTerm.name
	if not soterm:
	    continue
	if not ("gene" in soterm or "pseudo" in soterm or "lnc_RNA" in soterm or "lncRNA" in soterm):
	    continue
	col3 = "pseudogene" if "pseudo" in soterm else "gene"
	s = f.chromosomeLocation.strand
	strand = "+" if s == "+1" else "-" if s == "-1" else "."
	dbxrefs = [ pnameMap[xr.source.name] + ":" + xr.identifier for xr in f.crossReferences ]
	#
	soterm = f.sequenceOntologyTerm.name
	if soterm in ["lncRNA","lnc_RNA"]:
	    soterm = "lncRNA_gene"
	elif soterm == "antisense_lncRNA":
	    # SO has no term for "antisense lncRNA gene"
	    soterm = "lncRNA_gene"
	# filter out genes with no model IDs
	if len(dbxrefs) > 0:
	    g = gff3.Feature([
		f.chromosomeLocation.locatedOn.primaryIdentifier,
		"MGI",
		col3,
		str(f.chromosomeLocation.start),
		str(f.chromosomeLocation.end),
		".",
		strand,
		".",
		{
		    "ID": f.primaryIdentifier,
		    "curie":  f.primaryIdentifier,
		    "gene_id": f.primaryIdentifier,
		    "so_term_name" : soterm,
		    "Name": f.symbol,
		    "description": f.name,
		    "Dbxref" : dbxrefs

		}
	    ])
	    feats.append(g)
    # end for loop

    feats.sort( lambda a,b: cmp((a.seqid,a.start), (b.seqid,b.start)) )
    for f in feats:
        sys.stdout.write(str(f))

#
main()


