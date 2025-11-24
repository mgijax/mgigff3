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
import re
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
    mm._marker_key,
    aa.accid, 
    mm.symbol,
    mm.name,
    mt.name as markertype,
    mmc.term as mcvtype,
    mlc.chromosome,
    mlc.genomicchromosome,
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

# returns regulates/regulated_by relationships (WTS2-1756/SPRT-128)
mgiRegulates = '''
    select distinct m2._marker_key, m2.symbol, m1.symbol as regulator, c.pubmedid, c.mgiid
    from MGI_Relationship r, MRK_Marker m1, MRK_Marker m2, BIB_Citation_Cache c
    where r._category_key = 1013
    and r._object_key_1 = m1._marker_key
    and r._object_key_2 = m2._marker_key
    and r._refs_key = c._refs_key
'''

SA_re = re.compile(r'(\d+)')
def smartAlphaKey (a) :
    key = []
    for ap in SA_re.split(a):
        if ap.isdigit():
            ap = int(ap)
        key.append(ap)
    return key

def main () :
    # Read the MGI/model id pairs. Build index from MGI id -> list of model ids
    ix = {}
    for rec in db.sql(mgiModelIds):
        xref = ldbkey2prefix[rec['_logicaldb_key']] + ':' + rec['modelid']
        ix.setdefault(rec['mgiid'],[]).append(xref)
    #
    rx = {}
    for rec in db.sql(mgiRegulates):
        refid = ('PMID:'+rec['pubmedid']) if rec['pubmedid'] else rec['mgiid']
        t = rec['regulator'] + '[Ref_ID:' + refid + ']'
        rx.setdefault(rec['_marker_key'],[]).append(t)
    # Read all MGI genes and pseudogenes, and convert them to GFF3.
    for rec in db.sql(mgiGenes):
      try:
          mgiid = rec['accid']
          mkey = rec['_marker_key']
          if not mgiid in ix:
            # only want things with model ids
            continue
          symbol = rec['symbol']
          name = rec['name']
          chr = rec['genomicchromosome'] if rec['genomicchromosome'] else rec['chromosome']
          start = '.' if rec['startcoordinate'] is None else int(rec['startcoordinate'])
          end   = '.' if rec['endcoordinate']   is None else int(rec['endcoordinate'])
          strand = rec['strand']
          if strand is None:
            sys.stderr.write("null strand: " + str(rec) + "\n")
            strand = '.'
          markertype = rec['markertype']
          mcvtype = rec['mcvtype']
          so_term_name = mcv2so.mcv2so[mcvtype]
          attrs = [
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
            ['so_term_name', so_term_name],
          ]
          if mkey in rx:
            attrs.append(['Regulated_by', sorted(rx[mkey], key=smartAlphaKey)])
          #
          gffrec = [
            chr,
            'MGI',
            markertype.lower(),
            str(start),
            str(end),
            '.',
            strand,
            '.',
            gff3.formatColumn9(dict(attrs))
          ]
          print('\t'.join(gffrec))
      except:
        print ("ERROR:")
        print (str(rec))
        sys.exit(1)


####
#
main()
