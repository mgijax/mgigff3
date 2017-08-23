#!/usr/bin/bash

PYTHON=python2.4
CUT=cut
SED=sed
GFCLIENT=gfClient
PSLREPS=pslReps

#${PYTHON} ${BIN}/prepMgiComputed.py | ${SORT} > ${WORKINGDIR}/mgiComputed.tsv
#${CUT} -f 2  ${WORKINGDIR}/mgiComputed.tsv | head -100 | ${PYTHON} seqid2fa.py - - > ${CACHEDIR}/sequences.fa
#${GREP} "^>" ${CACHEDIR}/sequences.fa | ${SED} 's/>\([^. ]*\).*/\1/' > ${CACHEDIR}/ids.txt

BLAT_HOST="bhmgiapp01.jax.org"
BLAT_PORT="9038"

#${GFCLIENT} -nohead -minIdentity=95 ${BLAT_HOST} ${BLAT_PORT} / ${CACHEDIR}/sequences.fa ${WORKINGDIR}/blatresults.psl
#${PSLREPS} -singleHit -nohead ${WORKINGDIR}/blatresults.psl ${WORKINGDIR}/pslreps.psl ${WORKINGDIR}/pslreps.psr
${PYTHON} phase2merge.py ${WORKINGDIR}/mgiComputed.tsv ${WORKINGDIR}/pslreps.psl > ${WORKINGDIR}/mgiComputed.gff


