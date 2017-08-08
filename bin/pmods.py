import lib.gff3 as gff3
import sys

for (grp, ix) in gff3.models(sys.argv[1]):
    for f in grp:
	sys.stdout.write(str(f))
