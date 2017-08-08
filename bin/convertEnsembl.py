
import sys
import lib.gff3 as gff3

for (model, x) in gff3.models(sys.stdin):
    for f in model:
        sys.stdout.write(str(f))
