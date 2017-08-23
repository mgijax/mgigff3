
import sys
import gff3

for m in gff3.models(sys.stdin, flatten=True):
  print m
