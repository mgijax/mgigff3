#!/usr/bin/bash

PYTHON=python
CURL=curl
CP=cp
GUNZIP=gunzip
GREP=grep

# https://stackoverflow.com/questions/59895/getting-the-source-directory-of-a-bash-script-from-within
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PYTHONPATH=${DIR}/lib:${PYTHONPATH:-.}
export PYTHONPATH

BIN=${DIR}
WORKINGDIR=${DIR}/../work

DATADIR=${DIR}/../data
#DATADIR=/data/research/mouse_build_38_external/original_annotations

SORTCMD='sort -k 1,1 -k 4n,4n -k 5nr,5nr'
SPLIT="${PYTHON} ${BIN}/splitGff.py -d ${WORKINGDIR}"

LOGFILE=${WORKINGDIR}/LOG
truncate -s 0 ${LOGFILE}

function logit {
    echo `date` $1 >> ${LOGFILE}
}

logit "Starting refresh..." 

########
# MGI
#logit 'prepMgi...'
#${PYTHON} ${BIN}/prepMgi.py | ${SORTCMD} | ${SPLIT} -t "mgi.chr%s.gff" 

########
# MGI computed
# logit "prepMgiComputed..."
#${PYTHON} ${BIN}/prepMgiComputed.py > ${WORKINGDIR}/mgi.gff3

########
# NCBI
#logit 'prepNcbi...'
#${PYTHON} ${BIN}/prepNcbi.py < ${DATADIR}/ref_GRCm38.p4_top_level.gff3 | ${SPLIT} -t "ncbi.chr%s.gff"

########
# miRBase
#logit 'prepMirbase...'
#${PYTHON} ${BIN}/prepMirbase.py < ${DATADIR}/miRBase21_mmu.gff3 | ${SPLIT} -t "mirbase.chr%s.gff"

########
# ENSEMBL
ENSEMBLver=89
ENSEMBLfile=Mus_musculus.GRCm38.${ENSEMBLver}.gff3
ENSEMBLurl=ftp://ftp.ensembl.org/pub/release-${ENSEMBLver}/gff3/mus_musculus/${ENSEMBLfile}.gz
#
#${CURL} ${ENSEMBLurl} | ${GUNZIP} > ${DATADIR}/${ENSEMBLfile}
#logit 'prepEnsembl...'
#${PYTHON} ${BIN}/prepEnsembl.py < ${DATADIR}/${ENSEMBLfile} | ${SPLIT} -t "ensembl.chr%s.gff"

########

chrs=( 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 X Y MT )
for i in "${chrs[@]}"
do : 
    logit
    logit "merging chr${i}..."
    ${PYTHON} ${BIN}/merge.py ${WORKINGDIR}/*.chr${i}.gff > ${WORKINGDIR}/chr${i}.gff 2>> ${LOGFILE}
done

########
logit 'Refresh finished.'
