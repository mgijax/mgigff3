
#
# mirbase2gff3.py
#
# The GFF3 file from mirBASE has two feature types: miRNA_primary_transcript and miRNA.
# miRNAs are mature products derived from primary transcripts. They have a "Derives_from" 
# attribute that points to their transcript.
# 
# The mapping to the Alliance/MGI.gff3 structure for miRNAs goes like this:
# - each miRNA_primary_transcript generates a gene (necessary but temporary placeholder - will
#   go away downstream at merge step). The transcript is made a child of the gene.
# - each miRNA_primary_transcript also generates an exon as its sole child
# - each miRNA is made a child of the gene of its primary transcript
#   (the Derives_from attribute is maintained)
# - each miRNA generates an exon as its sole child
#
#
import sys
import gff3

seenIds = set()
currGene = None
currTrans = None

for f in gff3.iterate( sys.stdin ):
    i = f.ID
    if i in seenIds:
        sys.stderr.write('DUPLICATE ID (skipped): ')
        sys.stderr.write(str(f))
        sys.stderr.write('\n')
        continue
    else:
        seenIds.add(i)

    if f[0].startswith('chr'):
        f[0] = f[0][3:]
    f[1] = "miRBase"

    f.attributes.pop("Alias", None)

    if f[2] == "miRNA_primary_transcript":
        #
        # Gene feature
        #
        g = currGene = gff3.Feature(f)
        g.ID += "_G"
        g[2] = 'gene'   # 
        g.attributes["curie"] = "miRBase:" + f.ID
        g.attributes["so_term_name"] = "miRNA_gene"
        sys.stdout.write(str(g))
        #
        # Transcript feature
        #
        currTrans = f
        f[2] = 'pre_miRNA'
        f.Parent = g.ID
        f.transcript_id = f.ID
        sys.stdout.write(str(f))
        #
        # Exon feature
        #
        e = gff3.Feature(f)
        e[2] = 'exon'
        e.attributes = {}
        e.ID = f.ID + "_exon"
        e.Parent = f.ID
        sys.stdout.write(str(e))
    elif f[2] == "miRNA": # per AGR
        #
        # miRNA feature
        #
        f.Parent = currGene.ID
        sys.stdout.write(str(f))
        #
        # Exon feature
        #
        e = gff3.Feature(f)
        e[2] = 'exon'
        e.attributes = {}
        e.Parent = f.ID
        e.ID = f.ID + "_exon"
        sys.stdout.write(str(e))
