#
# gffSort.py
#
# Usage:
#   python3 gffSort.py < input.gff > output.gff
#

import sys
from lib import gff3

gffStream = gff3.iterate(sys.stdin, returnHeader=True, returnGroups=True)
header = next(gffStream)
allModels = list(gffStream)
# sort by models by chromosome, then start position of the top-level feature
allModels.sort(key = lambda m: (m[0][0], m[0][3]))
for model in allModels:
    for f in model:
        print(gff3.format(f), end='')
    print('###')

