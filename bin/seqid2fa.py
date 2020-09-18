#
# seqid2fa.py
#
# Fetches sequences from entrez by id. Input is a single-column file of ids.
# Output is a fasta-formatted file of sequences from Genbank (technically, 
# from the nucleotide database accessible via NCBI eutils).
# Usage:
#       % python seqid2fa.py INFILE OUTFILE
#
#

import sys
import os
import types
import time
import urllib.request, urllib.parse, urllib.error
import xml.dom.minidom

# FIXME. Move all these into config.sh and access via os.environ
FETCHURL  = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi'
BATCHSIZE = 500
SLEEPTIME = 3
TOOL      = "MGI"
EMAIL     = "Joel.Richardson@jax.org"
NTRIES    = 3

def log (msg) :
    sys.stderr.write(msg)
    sys.stderr.flush()

cfile = None
cfd = None
def processLine (line, yielded, cacheDir) :
    global cfile
    global cfd
    if line.startswith(">"):
        ident = line.split()[0][1:]
        identNoV = ident.split(".")[0]
        yielded.add(ident)
        yielded.add(identNoV)
        if cacheDir:
            cfile = os.path.join(cacheDir, identNoV + '.fa')
            if cfd:
                cfd.close()
            cfd = open(cfile, 'w')
    if cfd:
        cfd.write(line)
    return line

def openSequenceFetch(
        ids,
        tool,
        email,
        db='nucleotide',
        retmode='text',
        rettype='fasta',
        batchsize=BATCHSIZE,
        sleeptime=SLEEPTIME,
        cacheDir=None):
    # 
    if cacheDir:
        os.makedirs(cacheDir, exist_ok=True)
        #
        # First yield all sequences currently in the cache
        ids2 = []
        for ident in ids:
            cachefile = os.path.join(cacheDir, ident + '.fa')
            if os.path.exists(cachefile):
                fd = open(cachefile, 'r')
                for line in fd:
                    yield line
                fd.close()
            else:
                ids2.append(ident)
        # Now go on and request the remainder from NCBI
        ids = ids2

    #
    for i in range(0, len(ids), batchsize):

        idBatch = ids[i:i+batchsize]
        yielded = set()

        # Get the next batch of ids.
        params = urllib.parse.urlencode(
           {'db': db,
            'retmode' : retmode,
            'rettype' : rettype,
            'id' : ",".join(idBatch),
            'tool' : tool,
            'email' : email
            }).encode()

        # For each batch, try up to NTRIES time to get the sequences. Provides
        # some protection from intermittant errors from the eUtils server.
        for ntry in range(NTRIES):
            try:
                time.sleep(sleeptime)
                log("Fetching sequence batch %d, %d sequences, attempt %d\n" % (i/batchsize, len(idBatch), ntry+1))
                fd = urllib.request.urlopen(FETCHURL, params)
            except:
                log("Error: %s\n" % str(sys.exc_info()[1]))
                continue

            line = fd.readline().decode('utf-8')
            if not line:
                log("Could not read first line.\n")
                break
            elif not line.startswith(">"):
                log("First line is not a Fasta header: %s\n" % line)
                fd.close()
                continue
            # success!
            yield processLine(line, yielded)
            for line in fd:
                yield processLine(line.decode('utf-8'), yielded)
            fd.close()
            for ident in idBatch:
                if not ident in yielded:
                    log("Requested seqid not returned: %s\n" % ident)
            break
        else:
            # if we exhaust the loop, there was an error
            log("Failed to get data from eUtils after %d tries.\n" % NTRIES)
            sys.exit(1)

def fetchSequences(
        infile,
        outfile,
        tool,
        email,
        cacheDir,
        db='nucleotide',
        retmode='text',
        rettype='fasta',
        batchsize=BATCHSIZE,
        sleeptime=SLEEPTIME):

    if infile == "-" or infile is sys.stdin:
        fd = sys.stdin
        infile = "<stdin>"
    else:
        fd = open(infile, 'r')
    lines = fd.readlines()
    fd.close()
    #
    ids = [l.split()[0] for l in lines]
    if outfile == "-" or outfile is sys.stdout:
        fd = sys.stdout
    else:
        fd = open(outfile, 'w')
    #
    for line in openSequenceFetch(ids,tool,email,db,retmode,rettype,batchsize,sleeptime,cacheDir):
        if line[0] == ">":
            seqidv = line[1:].split()[0] # everything after ">" up to the first whitespace char
            seqid = seqidv.split('.')[0]
        fd.write(line)
    #
    if fd is not sys.stdout:
        fd.close()

if __name__ == "__main__":
    if len(sys.argv) == 1:
        fetchSequences('-', '-', TOOL, EMAIL, None)
    elif len(sys.argv) == 2:
        fetchSequences(sys.argv[1], '-', TOOL, EMAIL, None)
    elif len(sys.argv) == 3:
        fetchSequences(sys.argv[1], sys.argv[2], TOOL, EMAIL, None)
    elif len(sys.argv) == 4:
        fetchSequences(sys.argv[1], sys.argv[2], TOOL, EMAIL, sys.argv[3])
