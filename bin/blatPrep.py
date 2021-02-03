#
# blatPrep.py
#
# Outputs sequence IDs for genes that do not have associated gene models.
#
# Sequences for a gene must satisfy:
# - from Genbank/EMBL/DDBJ or RefSeq only
# - only DNA and RNA (no protein)
# - if DNA, must be < 10 kb in length
# - must be associated with only one gene.
# - must not be a strain-specific gene
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

import sys
from lib import mgiadhoc as db
from lib import mcv2so

# Genes and pseudogenes.
# Exclude miRNA and tRNA genes.
# Exclude genes with a strain-specificity note
genes = '''
  SELECT m._marker_key
  FROM MRK_Marker m, MRK_MCV_Cache mcv
  WHERE m._marker_status_key = 1
  AND m._organism_key = 1
  AND m._marker_type_key in (1,7)
  AND m._marker_key = mcv._marker_key
  AND mcv.qualifier = 'D'
  AND mcv.term not in ('miRNA gene', 'tRNA gene')
  AND m._marker_key not in (
      select _object_key
      from mgi_note
      where _notetype_key = 1035)
  '''

# Sequences in MGI associated with exactly one gene
seqsWithOneGene = '''
  SELECT _sequence_key
  FROM SEQ_Marker_Cache
  GROUP BY _Sequence_key
  HAVING count(*) = 1
  '''

# All qualified sequences in MGI
# - associated with exactly one gene (gene or pseudogene, no miRNA, no tRNA)
# - type = DNA or RNA
# - logical db = Genbank or RefSeq
# - if type = DNA, then length <= 10 kb
# - not low quality
# - not in problem alignment set
sequences = '''
  SELECT distinct
      c.accid as sequenceid,
      st.term as seq_type,
      s.division,
      a.accid as markerid,
      m.symbol,
      m.name,
      v.term as mcv_type,
      mc.chromosome
  FROM
      SEQ_Marker_Cache c,
      SEQ_Sequence s,
      ACC_Accession a,
      MRK_Marker m,
      MRK_MCV_Cache v,
      VOC_Term st,
      MRK_Location_Cache mc
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
  AND c._marker_key = mc._marker_key
  AND s._sequencetype_key = st._term_key
  /* sequence quality != low */
  AND s._sequencequality_key != 316340
  /* Seq assoc with only one gene */
  AND s._sequence_key in (%s)
  /* RNA or DNA <= 10kb in len */
  AND (s._sequencetype_key = 316346 OR (s._sequencetype_key = 316347 AND s.length <= 10000))
  /* not in problem alignment set */
  AND c._sequence_key not in (
      SELECT _object_key
      FROM MGI_SetMember
      WHERE _Set_key = 1058)
  ''' % (genes, seqsWithOneGene)

genesWithModels = set([ mid.strip() for mid in sys.stdin ])

for r in db.sql(sequences):
    if r["markerid"] in genesWithModels:
        continue
    mcvtype = r["mcv_type"]
    sotype  = mcv2so.mcv2so[mcvtype]
    line = [
      r["sequenceid"],
      r["seq_type"],
      r["division"],
      r["markerid"],
      r["symbol"],
      r["name"],
      mcvtype,
      sotype,
      r["chromosome"]
    ]
    print("\t".join(line))
