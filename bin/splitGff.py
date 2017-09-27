#
# splitGff.py
#
# Splits a gff file containing multiple chromosomes into multiple files each containing one chromosome.
#

import sys
import gff3
import argparse

parser = argparse.ArgumentParser(description='Split a GFF file into multiple files, based on chromosome.')
parser.add_argument('-d', '--dir', default='.', help='template for naming the files')
parser.add_argument('-t', '--tmplt', default='chr%s.gff', help='template for naming the files')

args = parser.parse_args()
gff3.splitFile(sys.stdin, args.dir, args.tmplt)

