import gff3
import sys
idMaker = gff3.IdMaker()
for m in gff3.models(sys.stdin):
    feats = gff3.flattenModel(m)
    gff3.reassignIDs(feats, idMaker)
    for f in feats:
        sys.stdout.write(str(f))
    sys.stdout.write("###\n")

