#
# getSOterms.py
#
#

import sys
import mgiadhoc as db

SO_vocab_key = 138
SO_logicaldb_key = 145

SO_terms = '''
  SELECT t.name
  FROM VOC_Term t
  WHERE t._vocan_key = %d
  LIMIT 25
''' % SO_vocab_key

for r in db.sql(SO_terms):
  print(r)
