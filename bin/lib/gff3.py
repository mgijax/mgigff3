#
# gff3.py
#
# Utility functions for working with GFF3 files.
# This is NOT a complete implementation of GFF3! (E.g., it
# does not do validation.)
#
# Example usage: Print the ID (col 9 attribute) and length
# of all the genes in a gff file.
#	import gff3
#	for feature in gff3.iterate("mydata.gff"):
#	    if feature.type == "gene":
#	        length = feature.end - feature.start + 1
#	        print feature.attributes['ID'], length
#	    
# Feature objects. This module defines a class called 'Feature'. 
# A Feature object corresponds to a line in a GFF3 file and provides
# convenient access to the field values, as well as access to attribute/value
# pairs in column 9.
#
# A feature can be created several ways:
#	- With no arguments, every field is initialized to ".", except
#	col 9, which is initialized to {}.
#	- With a string argument (e.g. a line from a GFF file), the 
#	string is parsed and the fields initialized accordingly.
#	- With a Feature object, the new Feature is a copy.
#	
# Fields of a feature can be accessed:
#	- by index: f[i], where is 0 <= i <= 8
#	  Beware that indexes are 0-based
#	  while common usage in GFF documentation/discussion is 1-based.
#	  E.g., f[2] is the type column (GFF column 3).
#	- by name: f.xxx, where xxx is one of the defined column names:
#	    seqid, source, type, start, end, score, strand, phase, attributes
#	  Thus, for example, f.start is equivalent to f[3]
#
# Features also provide support for attributes in column 9.
# Attributes are maintained in a dict that maps names
# to values. A value is either a string or a list of strings.
# You can directly manipulate the attributes as a normal Python dict,
# e.g.,
#	if not f.attributes.has_key("Name"):
#	    f.attributes["Name"] = "foo"
#
# As long as an attribute name (1) does not conflict with a 
# predefined field name (above) and (2) is a valid python
# identifier, you can also access attributes as if they were
# direct attributes of the feature, e.g., f.Name is
# equivalent to f.attributes["Name"]
#
# CDS features. In most cases, the one feature / one GFF3 line correspondence
# works. CDS features are one notable exception. (There may be others. I don't know.)
# A single (real) CDS comprises a _set_ of GFF3 lines because a real CDS comprises
# a sequence of spliced bits of the genome. In GFF3, this becomes a sequence of individual
# lines that describe the bits all with the same ID and type. The general rule in GFF3
# is that it is OK for one "feature" to consist of multiple discontiguous pieces.
#
# The standard functions hasattr, getattr, and setattr work
# consistently with the above semantics.
#
#----------------------------------------------------
import sys
import types
import urllib
import re
from OrderedSet import OrderedSet

#----------------------------------------------------
HEADER = '##gff-version 3\n'

C9SEP = ';'
QUOTECHARS_RE = re.compile(r'[\t\n\r;=%&,]')
COMMENT_CHAR = '#'
GROUPSEP = "###\n"

#----------------------------------------------------
HASH	= '#'
TAB	= '\t'
NL	= '\n'
SP	= ' '
SEMI	= ';'
EQ	= '='
COMMA	= ','

WSP_RE = re.compile(r'^\s*$')

#----------------------------------------------------
#
class ParseError(RuntimeError):
    pass

#----------------------------------------------------
#
# A Feature encapsulates one row of a GFF3 file. The 9 fields
# of a Feature can be accessed either by index (0-9) or by name,
# as defined by the GFF3 definition document at:
#	http://www.sequenceontology.org/gff3.shtml
#
# The Feature constructor can be called in several ways:
#	- with no arguments. All fields initialized to the
#	empty value ("."), except the 9th, which is initialized
#	to an empty dict.
#	- a string, i.e., a line from a GFF3 file. (NOTE: 
#	comment lines are NOT accepted.) The string is parsed.
#	- a list of values. The list must have length 9, and
#	the 9th value must either be a string (which will be
#	parsed) or a dict of name->value mappings.
#	- another Feature. The new feature is a copy.
# 
# A feature is a List of 9 values corresponding to the 9 GFF3 columns.
# A feature can be accessed by index.
# A feature can also be accessed by the 9 defined field names:
#  seqid, source, type, start, end, score, strand, phase, attributes.
# For example, for feature f, f[0] is equivalent to f.seqid, etc.
# The 9th column (f[8] or f.attributes) is a dict containing 
# name->value mappings. All names are strings. A value is either a
# string, or a list of strings. For example, 
#	f.attributes['Name']='foo'
# sets the feature's Name attribute to foo. 
# Attributes can also be accessed directly as attributes of the
# feature itself, e.g.,
#	f.Name = "foo"
#
# In general, the expression
#	f.xxx
# is equivalent to
#	f.attributes["xxx"]
# provided "xxx" in not one of the 9 GFF3 column names.
#
class Feature(types.ListType):

    # These are the standard field names.
    fields = [
	"seqid",
	"source",
	"type",
	"start",
	"end",
	"score",
	"strand",
	"phase",
	"attributes",
	]
    
    # a dict that maps field names to indices
    field2index = dict(map(lambda a : (a[1],a[0]), enumerate(fields)))

    #
    def __init__(self, arg=None):
	if arg is None:
	    arg = ['.'] * 9
	    arg[-1] = []
	elif type(arg) is types.StringType:
	    arg = parse(arg)
	elif len(arg) != 9:
	    raise ValueError("Invalid initializer for GFFFeature: " \
		+ (" %d fields\n" % len(arg)) + str(arg))
	types.ListType.__init__(self,arg)
	#
	t8 = type(self[8])
	if t8 is types.StringType:
	    self[8] = parseColumn9(self[8])
	elif t8 is types.ListType:
	    self[8] = dict(self[8])
	elif t8 is types.DictType:
	    # Copy the dict.
	    # If there are list valued attrs, make sure we
	    # don't share the list obj.
	    d = {}
	    for k,v in self[8].iteritems():
		if type(v) is types.ListType:
		    d[k] = v[:]
		else:
		    d[k] = v
	    self[8] = d
	if self.start != ".":
	    self.start = int(self.start)
	if self.end != ".":
	    self.end = int(self.end)

    def __hash__(self):
	return hash(self.attributes.get('ID',None))

    def __getattr__(self, name):
	i = Feature.field2index.get(name,None)
	if i is None:
	    v = self[8].get(name,None)
	    if v is None:
		raise AttributeError(name)
	    else:
		return v
	else:
	    return self[i]

    def __setattr__(self, name, value):
	if (name=="start" or name=="end") and value != ".":
	    value = int(value)
	i = Feature.field2index.get(name,None)
	if i is None:
	    self[8][name] = value
	else:
	    self[i]=value
    
    def __str__(self):
	return format(self)

#----------------------------------------------------
# A very simple file iterator that yields a sequence
# of GFF3 Features. 
#
# Args:
#  input (file name or open file) If file name is "-", reads
#	from standard input.
#  returnGroups (boolean) If True, groups Features into lists
#	before yielding. This only makes sense if the GFF3 file
#	uses the "###" construct. (See GFF3 spec.) If False,
#	(the default), yields each Feature individually.
#	If returnGroups is true and the input has no group separator lines,
#	all the features are returned as a single list. (Don't do this if
#	there are too many features.
#  returnHeader (boolean) If True, returns the (possibly empty) list of comment
#	lines at the top of the file as the first element in the iteration.
#	Default is False (first element is first feature or group).
#
def iterate(source, returnGroups=False, returnHeader=False):
    #
    # Set up the input
    #
    closeit = False
    if type(source) is types.StringType:
	if source=="-":
	    source = sys.stdin
	else:
	    source = open(source, 'r')
	    closeit = True
    group = []
    #
    # Iterate through file.
    #
    header = ([] if returnHeader else None)
    for lineNum, line in enumerate(source):
	if line.startswith(COMMENT_CHAR):
	    if header is not None:
	        header.append(line)
	    elif returnGroups and line == GROUPSEP and len(group) > 0:
		yield group
		group = []
	    else:
		continue
	else:
	    if header is not None:
	        yield header
		header = None
	    try:
		f = Feature(line)
	    except:
	        raise RuntimeError("GFF3 parse error at line %d:\n%s" % (lineNum+1, line))
	    if returnGroups:
		group.append(f)
	    else:
		yield f

    if returnGroups and len(group) > 0:
	yield group
	group = []

    #
    # Close input.
    #
    if closeit:
	source.close()

#----------------------------------------------------
def xmodels(source):
    for grp in iterate(source, returnGroups=True):
	crossReference(grp)
        yield grp

#----------------------------------------------------
# Iterator that yields a sequence of models. Each model is the sequence of features 
# that make it up. The ID/Parent relationships within the group forms a tree (or a dag)
# with the first element being the root.
# ASSUMES: the incoming features are assumed to be sorted;
# 
# Args:
#   features	the individual gff3 features, in order
# Yields:
#   sequence of models, each of which is a list of Features.
#
def models(features):
    if type(features) is types.StringType \
    or type(features) is types.FileType:
        features = iterate(features)
    #
    id2feature = {}
    models = []
    #
    def addChild(p, c):
	p.children.add(c)
	c.parents.add(p)
    #
    def overlaps( f1, f2 ):
        return f1.seqid == f2.seqid and f1.end >= f2.start and f1.start <= f2.end
    #
    def sortCmp(f1, f2):
	# root comes first
	if not len(f1.parents) and len(f2.parents): return -1
	if len(f1.parents) and not len(f2.parents): return 1
	# leaves come last
	if len(f1.children) == 0 and len(f2.children) > 0: return 1
	if len(f1.children) > 0 and len(f2.children) == 0: return -1
	# sort by start position
	c = cmp(f1.start, f2.start)
	if c: return c
	# 
	c = cmp(f1.end, f2.end)
	return -c
    #
    def flushFeature(f):
	id2feature.pop(f.attributes.get('ID',None),None)
	return f
    #
    def flushModel(f):
        m = []
	def _(f):
	    m.append(f)
	    for c in f.children:
	        _(c)
	_(f)
	m.sort(sortCmp)
	return [flushFeature(x) for x in m]

    # Flushes models. If no argument, flushes all models. If a Feature is passed, flush
    # models based on comparison with that feature.
    def flush(f=None):
	flushed = []
	newmodels = []
	for r in models:
	    if not f or not overlaps(r, f):
	        flushed.append(flushModel(r))
	    else:
	        newmodels.append(r)
        models[:] = newmodels
	return flushed

    # Main loop. Iterate over input features. 
    # Each feature either: starts a new model,
    # or attaches to an existing model.
    # Each feature may trigger some models to
    # be flushed.
    for i,f in enumerate(features):
	# print ">>>", str(f)
	# attach direct xref attributes for each 
	f.parents = OrderedSet()
	f.children = OrderedSet()
        if hasattr(f, 'ID'):
	    id2feature[f.ID] = f
        if hasattr(f, 'Parent'):
	    pts = [f.Parent] if type(f.Parent) is types.StringType else f.Parent
	    for pid in pts:
		try:
		    p = id2feature[pid]
		    addChild(p, f)
		except KeyError:
		    raise RuntimeError("Orphan detected. Feature #%di\n%s" % (i, str(f)))
	else:
	    models.append(f)
	    # only try to flush if current feature is a root
	    for m in flush(f):
		yield m
    # flush all remaining models
    for m in flush():
        yield m

#----------------------------------------------------
# Builds and returns an index from feature.ID to feature.
# Args:
#  features: (enumerable) the features to index
#  id2feature (dict, optional) If provided, adds entries. 
#		Otherwise, creates a new index.
# Returns:
#  A dictionary { ID -> Feature }
#
def index(features, id2feature=None):
    id2feature = id2feature if id2feature else {}
    for f in features:
	id = f.attributes.get("ID",None)
	if id:
	    id2feature[id] = f
    return id2feature

#----------------------------------------------------
# Turns all Parent attributes in a group of features in actual object references.
# Each feature gets two new attributes: parents and children: each feature points to
# 0 or more parent features, and is the parent of 0 or more children.
# Parents and children attributes are implemented as OrderedSets.
# 
# All referenced Parent objects must be among the given set of features (error otherwise).
#
def crossReference(features):
    id2feature = index(features)
    #
    for f in features:
	f.parents = OrderedSet()
	f.children = OrderedSet()
	pIds = f.attributes.get("Parent",[])
	pIds = ([pIds] if type(pIds) is types.StringType else pIds)
	for pid in pIds:
	    parent = id2feature[pid]
	    f.parents.add(parent)
	    parent.children.add(f)
    return id2feature

#----------------------------------------------------
# Parses one line from a GFF3 file.
# Returns None if the line is a comment line. Otherwise,
# returns a list of values. If parseCol9 is True,
# (the default), the 9th column is parsed into a dict
# of name-value mappings. If False, the 9th column
# is the unparsed string.
#
def parse(line, parseCol9=True):
    if line[0:1] == HASH:
        return None
    tokens = line.split(TAB)
    if len(tokens) != 9:
        raise ParseError("Wrong number of columns (%d)\n%s" % (len(tokens),line))
    if tokens[8][-1:] == NL:
        tokens[8] = tokens[8][:-1]
    if parseCol9:
	tokens[8] = parseColumn9(tokens[8])
    return tokens

#----------------------------------------------------
# Parses a string of name-value attributes, as defined by GFF3. 
# Returns the corresponding dictionary. 
# 
def parseColumn9(value):
    if value == ".":
	return {}
    c9 = {}
    for t in value.split(SEMI):
	if WSP_RE.match(t):
	    continue
	tt = t.split(EQ)
	if len(tt) != 2:
	    raise ParseError("Bad column 9 format near '%s'."%t)
	n = unquote(tt[0].strip())
	v = map(unquote, tt[1].strip().split(COMMA))
	if len(v) == 1:
	    v = v[0]
	c9[n] = v
    return c9

#----------------------------------------------------
#
# Substitutes the %XX hex code for the special characters: 
# tab, newline, formfeed, ampersand, equals, semicolon,
# percent, comma.
#
def quote(v):
    return QUOTECHARS_RE.sub(lambda m:"%%%0x"%ord(m.group(0)), str(v))

#----------------------------------------------------
#
# Unquotes all hex quoted characters.
#
def unquote(v):
    return urllib.unquote(str(v)) 

#----------------------------------------------------
#
PRE = ['ID','Name','Parent']
EXCLUDE = ['parents','children']

#----------------------------------------------------
#
# Formats a dictionary of name/value pairs appropriately
# for column 9.
#
def formatColumn9(vals):
    if type(vals) is types.StringType:
	return quote(vals)
    parts = []
    for n in PRE:
	x = vals.get(n, None)
	if x:
	    parts.append(formatAttribute(n,x))
    for n,v in vals.iteritems():
	if n not in PRE and n not in EXCLUDE:
	    parts.append(formatAttribute(n,v))
    ret = C9SEP.join(parts)
    return ret

#----------------------------------------------------
#
# Formats one name/value pair appropriate for inclusion in
# column 9.
#
def formatAttribute(n, v):
    if type(v) is types.ListType:
	return "%s=%s" % (quote(n), COMMA.join(map(quote,v)))
    else:
	return "%s=%s" % (quote(n), quote(v))

#----------------------------------------------------
#
# Formats a list into a GFF3 line (newline included)
#
def format(tokens):
    lt = len(tokens)
    if lt > 9:
	tokens2 = tokens[0:9]
    elif lt < 9:
	tokens2 = tokens + ['.']*(lt-9)
	tokens2[-1] = {}
    else:
	tokens2 = tokens[:]
    tokens2[8] = formatColumn9(tokens[8])
    return TAB.join(map(str,tokens2)) + NL

#----------------------------------------------------
#
#
if __name__=="__main__":
    def printeval(expr, ns):
	v = eval(expr,ns)
	print expr, "\t=", v
	return v

    def selftest():
	f = printeval('Feature()', globals())
	f = printeval("Feature('12	MGI	gene	12345678	12347890	.	+	.	ID=MGI:222222;Name=Abc')", globals())
	ns = {'f':f}
	printeval('f[0]', ns)
	printeval('f[0:4]', ns)
	printeval('f.seqid', ns)
	printeval('f.end - f.start', ns)
	printeval('f.attributes', ns)
	printeval('f.attributes["ID"]', ns)
	printeval('f.ID', ns)
	printeval('Feature.fields', globals())
	return f

    if len(sys.argv) == 2:
	if sys.argv[1] == "-":
	    fd = sys.stdin
	else:
	    fd = open(sys.argv[1])
	for f in iterate(sys.argv[1], returnHeader=True):
	    print f,
	fd.close()
    else:
	f=selftest()
