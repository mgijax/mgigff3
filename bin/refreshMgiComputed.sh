#!/usr/bin/bash

source config.sh

# For the MGI genes that don't have a gene model from a provider, we generate ersatz models by:
# 1. Finding "good" sequences
# 2. Blatting those sequences against the mouse genome
# 3. Filtering the blat results and turning them into little models.
#

#
# 1. generate a file of mgi_id/seq_id pairs.
logit "Phase 2: getting MGI data..."
${PYTHON24} ${BIN}/prepMgiComputed.py | ${SORT} > ${WORKINGDIR}/mgiComputed.seqids.txt
checkExit

exit

# 2a. feed the sequence ids (2nd column) to the fetch script
# this gets the sequences from ncbi and writes to a fasta file
${CUT} -f 1  ${WORKINGDIR}/mgiComputed.seqids.txt | ${PYTHON} seqid2fa.py - ${WORKINGDIR}/mgiComputed.seqs.fa

# 2b. Blat the sequences against the mouse genome assembly
${GFCLIENT} -nohead -minIdentity=95 ${BLAT_HOST} ${BLAT_PORT} / ${WORKINGDIR}/mgiComputed.seqs.fa ${WORKINGDIR}/mgiComputed.blat.psl

# 3a. Use pslreps to filter the results to single best hits
${PSLREPS} -singleHit -nohead ${WORKINGDIR}/mgiComputed.blat.psl ${WORKINGDIR}/mgiComputed.pslreps.psl ${WORKINGDIR}/mgiComputed.pslreps.psr

# Extracts the list of sequence ids for records in a fasta file
#${GREP} "^>" ${WORKINGDIR}/mgiComputed.seqs.fa | ${SED} 's/>\([^. ]*\).*/\1/' > ${WORKINGDIR}/ids.txt

# 3b. Convert the alignments to GFF3 hierarchies and attach the corresponding MGI id
${PYTHON} mgiComputedMerge.py ${WORKINGDIR}/mgiComputed.pslreps.psl ${WORKINGDIR}/mgiComputed.seqids.txt | ${SPLITCMD} -t "mgicomputed.chr%s.gff"


