#
# prepMgi.py
#
# Generates a GFF3 file of genes and pseudogenes in MGI.
# Writes to stdout.
# Only gene/pseudogene features; no subfeatures (obviously).
# These genes/pseudogenes will form the roots of model hierarchies in the 
# output. 
# Each feature includes the MGI id, symbol, name, actual SO type (e.g., tRNA_gene) 
# and genomic coordinates.
# Each feature also includes Dbxrefs to their various gene models.
# These Dbxrefs will direct the merging of models downstream.
#

import sys
from lib import mgiadhoc as db
from lib import gff3
from lib import mcv2so

geneModelLdbKeys = '59,60,83' # Ensembl gene model, NCBI gene model, miRBase
ldbkey2prefix = {
  59: 'NCBI_Gene',
  60: 'ENSEMBL',
  83: 'miRBase'
}

mgiGenes = '''
  select 
    aa.accid, 
    mm.symbol,
    mm.name,
    mt.name as markertype,
    mmc.term as mcvtype,
    mlc.chromosome,
    mlc.startcoordinate,
    mlc.endcoordinate,
    mlc.strand
  from
    mrk_marker mm,
    acc_accession aa,
    mrk_location_cache mlc,
    mrk_mcv_cache mmc,
    mrk_types mt
  where mm._marker_key = aa._object_key
    and aa._logicaldb_key = 1
    and aa._mgitype_key = 2
    and aa.private = 0
    and aa.preferred = 1
    and mm._marker_key = mmc._marker_key
    and mmc.qualifier = 'D'
    and mm._marker_key = mlc._marker_key
    and mm._marker_type_key = mt._marker_type_key
    and mt.name in ('Gene','Pseudogene')
  order by mlc.chromosome, mlc.startcoordinate
'''

# returns MGI id / provider model id pairs
mgiModelIds = '''
  select 
    aa.accid as mgiid,
    aa2.accid as modelid,
    aa2._logicaldb_key
  from
    acc_accession aa,
    acc_accession aa2
  where aa._logicaldb_key = 1
    and aa._mgitype_key = 2
    and aa.private = 0
    and aa.preferred = 1
    and aa._object_key = aa2._object_key
    and aa2._logicaldb_key in (%s)
    and aa2._mgitype_key = 2
    and aa2.preferred = 1
''' % geneModelLdbKeys

def main () :
    # Read the MGI/model id pairs. Build index from MGI id -> list of model ids
    ix = {}
    for rec in db.sql(mgiModelIds):
      xref = ldbkey2prefix[rec['_logicaldb_key']] + ':' + rec['modelid']
      ix.setdefault(rec['mgiid'],[]).append(xref)
    # Read all MGI genes and pseudogenes, and convert them to GFF3.
    for rec in db.sql(mgiGenes):
      mgiid = rec['accid']
      if not mgiid in ix:
        # only want things with model ids
        continue
      symbol = rec['symbol']
      name = rec['name']
      chr = rec['chromosome']
      start = '.' if rec['startcoordinate'] is None else int(rec['startcoordinate'])
      end   = '.' if rec['endcoordinate']   is None else int(rec['endcoordinate'])
      strand = rec['strand']
      if strand is None:
        sys.stderr.write("null strand: " + str(rec) + "\n")
        strand = '.'
      markertype = rec['markertype']
      mcvtype = rec['mcvtype']
      so_term_name = mcv2so.mcv2so[mcvtype]
      attrs = dict([
        # mint a B6 strain-specific ID from the MGI ID
        ['ID', 'MGI_C57BL6J_' + mgiid[4:]],
        # by convention, GFF3 Name attr is what a browser displays
        ['Name', symbol],
        # also by convention, GFF3 description attr used for "long" names
        ['description', name.replace(';',',')],
        # convention
        ['gene_id', mgiid],
        # for the Alliance
        ['curie', mgiid],
        # list the model ids
        ['Dbxref', ix[mgiid]],
        # what MGI calls it
        ['mgi_type', mcvtype],
        # a bona fide SO term corresponding to mgi_type
        # Throw an error if there is no mapping for this MCV term.
        ['so_term_name', so_term_name]
      ])
      gffrec = [
        chr,
        'MGI',
        markertype.lower(),
        str(start),
        str(end),
        '.',
        strand,
        '.',
        gff3.formatColumn9(attrs)
      ]
      print('\t'.join(gffrec))

####
#
main()
