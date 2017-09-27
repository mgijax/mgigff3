#
# trimForAgr.py
#
import sys
import gff3

sys.stdout.write(gff3.HEADER)
for m in gff3.models(sys.stdin):
    if m.type == "gene":
	if m.so_term_name == "protein_coding_gene":
	    # protein coding genes. Trim any parts that don't conform to the spec:
	    # - 3 levels
	    # - only mRNA or transcript at middle level
	    # - only exon or CDS at leaf
	    ttypes = ["transcript", "mRNA"]
	    ltypes = ["exon","CDS"]
	elif m.so_term_name == "miRNA_gene":
	    # miRNA genes: restricted to gene->pre_miRNA->miRNA
	    ttypes = ["pre_miRNA"]
	    ltypes = ["miRNA"]
	else:
	    # else any mid- or leaf-level types are ok
	    ttypes = None
	    ltypes = None
	#
	for t in list(m.children):
	    for e in list(t.children):
		if ltypes and not e.type in ltypes:
		    t.children.remove(e)
	    if len(t.children) == 0 or (ttypes and not t.type in ttypes):
		m.children.remove(t)
	#
	for f in gff3.flattenModel2(m):
	    sys.stdout.write(str(f))
    else:
	# pseudogenes
        sys.stdout.write(str(m))
