#
# trimForAgr.py
#
import sys
import os
import gff3
import time

###
# Loads a tab-delimited file of SO its and terms.
# Returns a mapping from term to id
#
def loadSOTerms () :
    soterm2id = {}
    fd = open(os.environ["SO_TERM_FILE"], 'r')
    for line in fd:
        soid,soterm = line.strip().split('\t')
        soterm2id[soterm] = soid
    fd.close()
    return soterm2id

###
def processModel (m, soterm2id) :
    # map so_term_name to its SO id and store it
    m.attributes["Ontology_term"] = soterm2id[m.attributes["so_term_name"]]
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
            else:
                try:
                    curie = t.source + ":" + t.attributes["transcript_id"]
                    t.attributes["curie"] = curie
                    t.attributes["Dbxref"] = curie
                    t.attributes["ontology_term"] = soterm2id[t.type]
                except KeyError:
                    m.children.remove(t)
        #
        for f in gff3.flattenModel2(m):
            sys.stdout.write(str(f))
    else:
        # pseudogenes - no substructure at alliance
        sys.stdout.write(str(m))

###
def writeHeader (timestamp, build) :
    sys.stdout.write(gff3.HEADER)
    sys.stdout.write("#!data-source MGI\n")
    sys.stdout.write("#!date-produced %s\n" % timestamp)
    sys.stdout.write("#!assembly %s\n" % build)

###
def main ():
    #
    soterm2id = loadSOTerms()
    #
    writeHeader(time.asctime(), os.environ["ENSEMBLbuild"])
    #
    for m in gff3.models(sys.stdin):
        processModel(m, soterm2id)

###
main()
