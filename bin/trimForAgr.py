#
# trimForAgr.py
#
import sys
import gff3

sys.stdout.write(gff3.HEADER)
for m in gff3.models(sys.stdin):
    if m.so_term_name != "protein_coding_gene":
        continue
    for t in list(m.children):
	for e in list(t.children):
	    if not e.type in ["exon","CDS"]:
	        t.children.remove(e)
        if len(t.children) == 0 or not t.type in ["transcript", "mRNA"]:
	    m.children.remove(t)
	    continue
    if len(m.children) == 0:
        continue
    for f in gff3.flattenModel(m):
        sys.stdout.write(str(f))
    #
    sys.stdout.write("###\n")
