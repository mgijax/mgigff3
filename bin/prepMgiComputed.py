#
# prepMgiComputed.py
#
# Outputs a 2-column file of
#	mgiid	seqid
# pairs, for genes that do not have associated gene models
# For each such gene, outputs one row for each sequence associated
# with the gene that also satisfies:
# - from Genbank/EMBL/DDBJ or RefSeq only
# - only DNA and RNA (no protein)
# - if DNA, must be < 10 kb in length
# - must be associated with only one gene.
# 

import mgiadhoc as db

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
  SELECT distinct c.accid as sequenceid, a.accid as markerid
  FROM SEQ_Marker_Cache c, SEQ_Sequence s, ACC_Accession a
  WHERE c._sequence_key = s._sequence_key
  AND c._LogicalDB_key in (9,27)
  AND c._marker_key in (%s)
  AND a._object_key = c._marker_key
  AND a._mgitype_key = 2
  AND a._logicaldb_key = 1
  AND a.preferred = 1
  /* Seq assoc with only one gene */
  AND s._sequence_key in (%s)
  /* RNA or DNA <= 10kb in len */
  AND (s._sequencetype_key = 316346 OR (s._sequencetype_key = 316347 AND s.length <= 10000))
  ''' % (genesWithoutModels, seqsWithOneGene)


for r in db.sql(gwomSequences):
    print r["markerid"] + "\t" + r["sequenceid"]
