#!/usr/bin/bash
#
# refreshMgiComputed.sh
#
# For the MGI genes that don't have a gene model from a provider, generates ersatz models.
# Writes to models to stdout.
# 
#
# 1. Finding "good enough" sequences for each gene
#	1a. Query MGI for the genes wihtout models and their qualifying seuqnece IDs
#	1b. Download the sequences from NCBI into a fasta file.
# 2. Blatting those sequences against the mouse genome
# 3. Filtering the blat results. 
#	3a. Only sequences with a single match across the genome 
#	3b. Only matches meeting %length and %identity thresholds
# 4. Turning the blat match results into little models
#	4a. Merges the MGI data from 1a with the filtered blat results from 3b
#	4b. Creates models with the structure: gene -> match -> match_part
# 5. Writing the results to GFF3 file
#
#

source config.sh

# 1. generate a file of mgi_id/seq_id pairs.
logit "Phase 2: getting MGI data..."
${PYTHON24} ${BIN}/prepMgiComputed.py > ${WORKINGDIR}/mgiComputed.seqids.txt
checkExit

# 2a. feed the sequence ids (2nd column) to the fetch script
# this gets the sequences from ncbi and writes to a fasta file
logit "Phase 2: retrieving sequences from NCBI..."
${CUT} -f 1  ${WORKINGDIR}/mgiComputed.seqids.txt | ${PYTHON} seqid2fa.py - ${WORKINGDIR}/mgiComputed.seqs.fa
checkExit

# 2b. Blat the sequences against the mouse genome assembly
logit "Phase 2: Blat'ing the sequences..."
${GFCLIENT} -nohead -minIdentity=95 ${BLAT_HOST} ${BLAT_PORT} / ${WORKINGDIR}/mgiComputed.seqs.fa ${WORKINGDIR}/mgiComputed.blat.psl >> ${LOGFILE} 2>&1
checkExit

# 3a. Use pslreps to filter the results to single best hits
logit "Phase 2: Filtering blat hits with pslreps..."
${PSLREPS} -singleHit -nohead ${WORKINGDIR}/mgiComputed.blat.psl ${WORKINGDIR}/mgiComputed.pslreps.psl ${WORKINGDIR}/mgiComputed.pslreps.psr >> ${LOGFILE} 2>&1
checkExit

# 3b. Convert the alignments to GFF3 hierarchies and attach the corresponding MGI id
logit "Phase 2: Generating models..."
${PYTHON} mgiComputedMerge.py ${WORKINGDIR}/mgiComputed.pslreps.psl ${WORKINGDIR}/mgiComputed.seqids.txt 
checkExit

logit "Phase 2: finished."


