#
# gffMerge.py
#
# Merges n sorted GFF3 files into one.
#
# Usage:
#   python3 gffMerge.py file1.gff3 file2.gff3 [...] > merged.gff3
#

import sys
from lib import gff3

def xRoot (src) :
    for model in src:
        model[0][8]["children"] = model[1:]
        yield model[0]

gffFiles = sys.argv[1:]
gffIters = [ xRoot(gff3.iterate(f, returnGroups=True)) for f in gffFiles ]

for m in gff3.merge(*gffIters):
    children = m[8].pop("children",[])
    sys.stdout.write(str(m))
    for f in children:
        sys.stdout.write(str(f))
sys.stdout.write("###\n")
