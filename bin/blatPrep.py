#
# prepMgiComputed.py
#
# Outputs sequence IDs for genes that do not have associated gene models.
#
# Sequences for a gene must satisfy:
# - from Genbank/EMBL/DDBJ or RefSeq only
# - only DNA and RNA (no protein)
# - if DNA, must be < 10 kb in length
# - must be associated with only one gene.
# 
# Outputs a 7-column file:
#       seqid   ID of the sequence
#       seqtype The type of the sequence. Either "DNA" or "RNA"
#       division The genbank division for the sequence
#       mgiid   MGI id of the associated gene
#       symbol  The gene's symbol
#       name    The gene's name
#       mcvtype The gene's MCV type
#       sotype  The gene's SO type

from lib import mgiadhoc as db

geneModelLdbKeys = '59,60,83,85'

genesWithModels = '''
  SELECT distinct m._marker_key
  FROM MRK_Marker m, ACC_Accession a
  WHERE m._marker_status_key = 1
  AND m._organism_key = 1
  AND m._marker_type_key = 1
  AND m._marker_key = a._object_key
  AND a._mgitype_key = 2
  AND a._logicaldb_key in (%s)
  ''' % geneModelLdbKeys

genesWithoutModels = '''
  SELECT m._marker_key
  FROM MRK_Marker m
  WHERE m._marker_status_key = 1
  AND m._organism_key = 1
  AND m._marker_type_key = 1
  AND m._marker_key not in (%s)
  ''' % genesWithModels

seqsWithOneGene = '''
  SELECT _sequence_key
  FROM SEQ_Marker_Cache
  GROUP BY _Sequence_key
  HAVING count(*) = 1
  '''

gwomSequences = '''
  SELECT distinct c.accid as sequenceid, st.term as seq_type, s.division, a.accid as markerid, m.symbol, m.name, v.term as mcv_type
  FROM SEQ_Marker_Cache c, SEQ_Sequence s, ACC_Accession a, MRK_Marker m, MRK_MCV_Cache v, VOC_Term st
  WHERE c._sequence_key = s._sequence_key
  AND c._LogicalDB_key in (9,27)
  AND c._marker_key in (%s)
  AND a._object_key = c._marker_key
  AND a._mgitype_key = 2
  AND a._logicaldb_key = 1
  AND a.preferred = 1
  AND c._marker_key = m._marker_key
  AND m._Marker_key = v._Marker_key
  AND v.qualifier = 'D'
  AND s._sequencetype_key = st._term_key
  /* Seq assoc with only one gene */
  AND s._sequence_key in (%s)
  /* RNA or DNA <= 10kb in len */
  AND (s._sequencetype_key = 316346 OR (s._sequencetype_key = 316347 AND s.length <= 10000))
  ''' % (genesWithoutModels, seqsWithOneGene)


typemap = {
  'antisense lncRNA gene'       : 'lncRNA_gene',
  'gene segment'                : 'gene_segment',
  'heritable phenotypic marker' : 'heritable_phenotypic_marker',
  'intronic lncRNA gene'        : 'lncRNA_gene',
  'lincRNA gene'                : 'lincRNA_gene',
  'lncRNA gene'                 : 'lncRNA_gene',
  'miRNA gene'                  : 'miRNA_gene',
  'protein coding gene'         : 'protein_coding_gene',
  'tRNA gene'                   : 'tRNA_gene',
  'unclassified gene'           : 'gene',
  'unclassified non-coding RNA gene' : 'ncRNA_gene',
  'sense intronic lncRNA gene' : 'lncRNA_gene',
}

for r in db.sql(gwomSequences):
    mcvtype = r["mcv_type"]
    sotype  = typemap.get(mcvtype, mcvtype)
    line = [
      r["sequenceid"],
      r["seq_type"],
      r["division"],
      r["markerid"],
      r["symbol"],
      r["name"],
      mcvtype,
      sotype,
    ]
    print("\t".join(line))
