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

# 2a. feed the sequence ids (2nd column) to the fetch script
# this gets the sequences from ncbi and writes to a fasta file
logit "Phase 2: retrieving sequences from NCBI..."
${CUT} -f 1  ${WORKINGDIR}/mgiComputed.seqids.txt | ${PYTHON} seqid2fa.py - ${WORKINGDIR}/mgiComputed.seqs.fa
checkExit

# 2b. Blat the sequences against the mouse genome assembly
logit "Phase 2: Blat'ing the sequences..."
${GFCLIENT} -nohead -minIdentity=95 ${BLAT_HOST} ${BLAT_PORT} / ${WORKINGDIR}/mgiComputed.seqs.fa ${WORKINGDIR}/mgiComputed.blat.psl
checkExit

# 3a. Use pslreps to filter the results to single best hits
logit "Phase 2: Filtering blat hits with pslreps..."
${PSLREPS} -singleHit -nohead ${WORKINGDIR}/mgiComputed.blat.psl ${WORKINGDIR}/mgiComputed.pslreps.psl ${WORKINGDIR}/mgiComputed.pslreps.psr
checkExit

# 3b. Convert the alignments to GFF3 hierarchies and attach the corresponding MGI id
logit "Phase 2: Generating models..."
${PYTHON} mgiComputedMerge.py ${WORKINGDIR}/mgiComputed.pslreps.psl ${WORKINGDIR}/mgiComputed.seqids.txt | ${SPLITCMD} -t "mgicomputed.chr%s.gff"
checkExit

logit "Phase 2: finished."


