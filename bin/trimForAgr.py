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
def log (s) :
    sys.stderr.write(s)
    sys.stderr.write('\n')

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
            ttypes = ["transcript", "mRNA", "lnc_RNA", "pseudogenic_transcript", "unconfirmed_transcript"]
            ltypes = ["exon","CDS"]
        elif m.so_term_name == "miRNA_gene":
            # miRNA genes: restricted to gene->pre_miRNA->miRNA
            ttypes = ["pre_miRNA", "miRNA"]
            ltypes = ["exon"]
        else:
            # else any mid- or leaf-level types are ok
            ttypes = None
            ltypes = None
        #
        for t in list(m.children):
            for e in list(t.children):
                if ltypes and not e.type in ltypes:
                    t.children.remove(e)
                # (schema-1.0.1.3) Make sure CDS feature has curie-form protein_id
                if e.type == "CDS" and "protein_id" in e.attributes:
                    source = "RefSeq" if e.source == "NCBI" else e.source
                    curie = source + ":" + e.attributes["protein_id"]
                    e.attributes["protein_id"] = curie
                #
                e.attributes.pop("transcript_id", None)
                e.attributes.pop("gene_id", None)
            if len(t.children) == 0 or (ttypes and not t.type in ttypes) or not t.attributes.get("transcript_id", None):
                if t.type != 'miRNA':
                    m.children.remove(t)
                    log("Removing %s" % str(t))
                else:
                    t.attributes.pop("gene_id", None)
            else:
                try:
                    # Transcripts already have non-curie transcript_id at this point.
                    #
                    # (schema-1.0.1.0) Make sure transcripts have curies
                    source = "RefSeq" if t.source == "NCBI" else t.source
                    curie = source + ":" + t.attributes["transcript_id"]
                    t.attributes["curie"] = curie
                    # (schema-1.0.1.3) Make sure transcripts have curie-form transcript_id in addition to curie
                    t.attributes["transcript_id"] = curie
                    #
                    t.attributes["Dbxref"] = curie
                    t.attributes["Ontology_term"] = soterm2id[t.type]
                    if t.type == "mRNA" and m.attributes["so_term_name"] != "protein_coding_gene" :
                        # for AGR, if a gene has a coding transcript, it must be protein_coding_gene
                        log("WARNING: Changing gene to protein coding because mRNA detected: %s %s" % (m.ID, m.attributes["Name"]))
                        m.attributes["so_term_name"] = "protein_coding_gene"
                        m.attributes["Ontology_term"] = "SO:0001217"
                    t.attributes.pop("gene_id", None)
                except KeyError:
                    log("KeyError (%s) for transcript: %s" % (sys.exc_info()[1],str(t)))
                    m.children.remove(t)
        #
        for f in gff3.flattenModel2(m):
            sys.stdout.write(str(f))
    else:
        # pseudogenes - no substructure at alliance
        sys.stdout.write(str(m))

###
def writeHeader (attrs) :
    sys.stdout.write(gff3.HEADER)
    for (n,v) in attrs:
        sys.stdout.write("#!%s %s\n" % (n, v))

###
def main ():
    log("Starting trimForAgr...")
    #
    soterm2id = loadSOTerms()
    #
    writeHeader([
        ('data-source', 'MGI'),
        ('date-produced', time.asctime()),
        ('assembly', os.environ["ENSEMBLbuild"]),
        ('annotationSource RefSeq', os.environ["NCBIver"]),
        ('annotationSource ENSEMBL', os.environ["ENSEMBLver"]),
        ])
    #
    for m in gff3.models(sys.stdin):
        processModel(m, soterm2id)

###
main()
