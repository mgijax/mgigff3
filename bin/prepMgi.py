
from intermine.webservice import Service
import gff3
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
    query.add_constraint("organism.taxonId", "=", "10090", code = "B")
    query.add_constraint("dataSets.name", "=", "Mouse Gene Catalog from MGI", code = "C")
    query.add_constraint("crossReferences.source.name", "ONE OF", 
         ["Ensembl Gene Model", "NCBI Gene Model", "miRBase"], code = "D")
    query.outerjoin("crossReferences")
    #
    return query

#
def main():
    for f in getFeatures():
	soterm = f.sequenceOntologyTerm.name
	if "gene" in soterm or "pseudo" in soterm or "lnc_RNA" in soterm or "lncRNA" in soterm:
	    col3 = "pseudogene" if "pseudo" in soterm else "gene"
	    s = f.chromosomeLocation.strand
	    strand = "+" if s == "+1" else "-" if s == "-1" else "."
	    dbxrefs = [ pnameMap[xr.source.name] + ":" + xr.identifier for xr in f.crossReferences ]
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
		    "so_term_name" : f.sequenceOntologyTerm.name,
		    "Name": f.symbol,
		    "description": f.name,
		    "Dbxref" : dbxrefs

		}
	    ])
	    print str(g),
#
main()
