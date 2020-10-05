#
# getSOterms.py
#
#

import sys
import mgiadhoc as db

SO_vocab_key = 138
SO_logicaldb_key = 145

SO_terms = '''
  SELECT t.term, a.accid
  FROM VOC_Term t, ACC_ACCESSION a
  WHERE t._vocab_key = %d
  AND t._term_key = a._object_key
  AND a._logicaldb_key = %d
  AND a.preferred = 1
''' % (SO_vocab_key, SO_logicaldb_key)

for r in db.sql(SO_terms):
  print(r['accid'] + '\t' + r['term'])
