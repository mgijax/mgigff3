#
# psl.py
#
# Classes for handling .psl format. 
#
# CONTENTS:
#    class Alignment: one alignment/line. Contains 21 attributes
#	representing the 21 fields of a psl line. (See attribute names
#	below.) Can also be accessed as a list (indices 0-20). Attributes
#	are integers, strings, or lists as appropriate. The three multivalued
#	fields (blockSizes, qStarts, tStarts) are represented as lists of integers.
#
#    iterate(): for iterating through a psl file. Each iteration returns next Alignment object.
# 
# PSL FORMAT:
# Psl is the format output by the Blat sequence alignment tool.
#
# Psl is a tab-delimited, tabular format.
# Each line represents one alignment, which may consist of multiple blocks.
# An alignment gives the areas (blocks) of correspondence between a query sequence and a target.
# Each line has 21 columns:
# 1. matches - Number of bases that match that aren't repeats
# 2. misMatches - Number of bases that don't match
# 3. repMatches - Number of bases that match but are part of repeats
# 4. nCount - Number of 'N' bases
# 5. qNumInsert - Number of inserts in query
# 6. qBaseInsert - Number of bases inserted in query
# 7. tNumInsert - Number of inserts in target
# 8. tBaseInsert - Number of bases inserted in target
# *9. strand - '+' or '-' for query strand. For translated alignments, second '+'or '-' is for genomic strand
# *10. qName - Query sequence name
# 11. qSize - Query sequence size
# 12. qStart - Alignment start position in query (0-based)
# 13. qEnd - Alignment end position in query (actually, 1 beyond last aligned base)
# *14. tName - Target sequence name
# 15. tSize - Target sequence size
# 16. tStart - Alignment start position in target
# 17. tEnd - Alignment end position in target
# *18. blockCount - Number of blocks in the alignment (a block contains no gaps)
# *19. blockSizes - Comma-separated list of sizes of each block
# 20. qStarts - Comma-separated list of starting positions of each block in query
# *21. tStarts - Comma-separated list of starting positions of each block in target
#
# *Asterisks indicate the columns we are usually concerned with
#
# A psl file may optionally have a 5-line header.
#

import sys
import types
import gff3

TAB	= '\t'
NL	= '\n'
COMMA	= ','

def psl2gff(seqlen,start,end,strand):
    return (start+1,end)
    '''
    if strand == "+":
      gstart = start + 1
      gend = end
    elif strand == "-":
      gstart = seqlen - end + 1
      gend = seqlen - start
    else:
        raise RuntimeError("Unrecognized strand value.")
    return (gstart,gend)
    '''

class Alignment(types.ListType):

    # These are the standard field names.
    fields = [
	"matches",	# 0
	"misMatches",	# 1
	"repMatches",	# 2
	"nCount",	# 3
	"qNumInsert",	# 4
	"qBaseInsert",	# 5
	"tNumInsert",	# 6
	"tBaseInsert",	# 7
	"strand",	# 8
	"qName",	# 9
	"qSize",	# 10
	"qStart",	# 11
	"qEnd",		# 12
	"tName",	# 13
	"tSize",	# 14
	"tStart",	# 15
	"tEnd",		# 16
	"blockCount",	# 17
	"blockSizes",	# 18
	"qStarts",	# 19
	"tStarts",	# 20
	]
    
    # a dict that maps field names to indices
    field2index = dict(map(lambda a : (a[1],a[0]), enumerate(fields)))

    #
    def __init__(self, arg=None):
	if arg is None:
	    arg = [0,0,0,0,0,0,0,0,'','',0,0,0,'',0,0,0,0,[],[],[]]
	elif type(arg) is types.StringType:
	    arg = self.__parse(arg)

	if len(arg) != len(self.fields):
	    raise ValueError("Invalid initializer for Alignment: " \
		+ (" %d fields\n" % len(arg)) + str(arg))
	types.ListType.__init__(self,arg)

    def __parse(self,line):
	fields = line.split(TAB)
	if fields[-1].endswith(NL):
	    fields[-1] = fields[-1][:-1]
	for (i,f) in enumerate(fields):
	    if i in [8,9,13]:
		pass
	    elif i >= 18:
		if f.endswith(COMMA):
		    f = f[:-1]
		fields[i] = map(int,f.strip().split(COMMA))
	    else:
		fields[i] = int(f)
	return fields

    def __getattr__(self, name):
	i = Alignment.field2index.get(name,None)
	if i is None:
	    raise AttributeError(name)
	else:
	    return self[i]

    def __setattr__(self, name, value):
	if (name=="start" or name=="end") and value != ".":
	    value = int(value)
	i = Alignment.field2index.get(name,None)
	if i is None:
	    raise AttributeError(name)
	else:
	    self[i]=value
    
    def __str__(self):
	x = self[:]
	for i in (18,19,20):
	    x[i] = COMMA.join(map(str, x[i])) + COMMA
	return TAB.join(map(str,x)) + NL

    def matchLength(self):
        return self.qEnd - self.qStart

    def percentIdentity(self):
        return (100.0 * self.matches) / self.matchLength()

    def percentLength(self):
        return (100.0 * self.matchLength()) / self.qSize

def iterate(input):
    #
    # Set up the input
    #
    closeit = False
    if type(input) is types.StringType:
	if input=="-":
	    input = sys.stdin
	else:
	    input = open(input, 'r')
	    closeit = True
    for a in input:
	if type(a) is types.StringType:
	    a = Alignment(a)
	yield a
    if closeit:
	input.close()

def toGff(input):
    # for notes on converting between psl and gff coordinates see: http://genome.ucsc.edu/FAQ/FAQformat.html#format2
    alignments = iterate(input)
    idMaker = gff3.IdMaker()
    for a in alignments:
	aid = idMaker.next("match")
	(gs,ge) = psl2gff(a.tSize, a.tStart, a.tEnd, "+") # the alignment coords are always w.r.t + strand
        root = gff3.Feature([
		a.tName,
		"BlatAlignment",
		"match",
		gs,
		ge,
		".",
		a.strand,
		'.',
	        {
		    'ID'    : aid,
		    'Name'  : a.qName,
		    'qName' : a.qName,
		    'matchLen' : a.matchLength(),
		    'pctIdentity' : a.percentIdentity(),
		    'pctLen' : a.percentLength(),
		    'qStart' : a.qStart,
		    'qEnd' : a.qEnd,
		    'qSize' : a.qSize,
		    'matches' : a.matches,
		}
               ])
        yield root 
	for i,tstart in enumerate(a.tStarts):
	    # alignment seg coordinates are w.r.t. their actual strand (i.e., neg strand coords run in reverse)
	    tend = tstart + a.blockSizes[i]
	    (gs,ge) = psl2gff(a.tSize, tstart, tend, a.strand)
	    part = gff3.Feature([
		a.tName,
		"BlatAlignment",
		"match_part",
		gs,
		ge,
		".",
		a.strand,
		'.',
		{
		    'ID'    : idMaker.next("match_part"),
		    'Parent': [ aid ],
		    'qName' : a.qName
		}
	       ])
	    yield part


    
if __name__ == "__main__":
    for a in toGff(sys.stdin):
	sys.stdout.write(str(a))
