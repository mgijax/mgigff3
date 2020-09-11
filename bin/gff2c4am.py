#
# Converts GFF3 features into a file suitable for input to the C4AM (coordinates for any marker) load.
#
# Reads from stdin and writes to stdout.
#
# Only outputs lines for top level features (anything without a Parent)
#
import sys
import gff3

for f in gff3.iterate(sys.stdin) :
    attrs = f.attributes
    if "Parent" in attrs and len(attrs["Parent"]) > 0:
        # f is not a top-level feature. Skip it.
        continue
    c4amRec = [
        attrs["curie"],
        f.seqid,
        str(f.start),
        str(f.end),
        "" if f.strand == "." else f.strand,
        "MGI_Blat",
        "MGI",
        ""
    ]
    sys.stdout.write("\t".join(c4amRec)+"\n")
