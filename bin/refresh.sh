#!/usr/bin/bash

PYTHON=python
CURL=curl
CP=cp
GUNZIP=gunzip
GREP=grep
SORT=sort
MKDIR=mkdir
DATE=date
TOUCH=touch

# The chromosomes and their order.
CHRS=( 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 X Y MT )

# https://stackoverflow.com/questions/59895/getting-the-source-directory-of-a-bash-script-from-within
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
BIN=${DIR}
DATADIR=${DIR}/../data
#DATADIR=/data/research/mouse_build_38_external/original_annotations
WORKINGDIR=${DIR}/../work
PYTHONPATH=${DIR}/lib:${PYTHONPATH:-.}

#
DATESTAMP=`${DATE} +"%Y-%m-%d"`

#
SORTCMD="sort -k 1,1 -k 4n,4n -k 5nr,5nr"
SPLITCMD="${PYTHON} ${BIN}/splitGff.py -d ${WORKINGDIR}"

#
${MKDIR} -p ${WORKINGDIR}

LOGFILE=${WORKINGDIR}/LOG.${DATESTAMP}
${TOUCH} ${LOGFILE}

function logit {
    echo `date` $1 >> ${LOGFILE}
}

export DIR BIN PYTHONPATH WORKINGDIR DATADIR LOGFILE SORT GREP SORTCMD SPLITCMD

nargs=$#

domgi=F
domgicomputed=F
doncbi=F
doensembl=F
domirbase=F
domerge=F
doexome=F

until [ -z "$1" ]  # Until all parameters used up . . .
do
    case "$1" in
    mgi)
	domgi=T
        ;;
    mgic*)
	domgicomputed=T
        ;;
    ncbi)
	doncbi=T
        ;;
    ensembl)
	doensembl=T
        ;;
    mirbase)
	domirbase=T
        ;;
    merge)
	domerge=T
        ;;
    exome)
	doexome=T
        ;;
    *)
        echo "Unrecognized option:" $1
	exit -1
    esac
    shift
done


########

logit "================================================================================="
logit "Starting refresh..." 

########
# MGI
if [ $nargs -eq 0 -o $domgi == T ]; then
    logit 'prepMgi...'
    ${PYTHON} ${BIN}/prepMgi.py 2>> ${LOGFILE} | ${SORTCMD} | ${SPLITCMD} -t "mgi.chr%s.gff" 
fi

########
# MGI computed
if [ $nargs -eq 0 -o $domgicomputed == T ]; then
    logit "prepMgiComputed..."
    refreshMgiComputed.sh 2>> ${LOGFILE}
fi

########
# NCBI
if [ $nargs -eq 0 -o $doncbi == T ]; then
    logit 'prepNcbi...'
    ${PYTHON} ${BIN}/prepNcbi.py 2>> ${LOGFILE} < ${DATADIR}/ref_GRCm38.p4_top_level.gff3 | ${SPLITCMD} -t "ncbi.chr%s.gff"
fi

########
# miRBase
if [ $nargs -eq 0 -o $domirbase == T ]; then
    logit 'prepMirbase...'
    ${PYTHON} ${BIN}/prepMirbase.py 2>> ${LOGFILE} < ${DATADIR}/miRBase21_mmu.gff3 | ${SPLITCMD} -t "mirbase.chr%s.gff"
fi

########
# ENSEMBL
if [ $nargs -eq 0 -o $doensembl == T ]; then
    ENSEMBLver=89
    ENSEMBLfile=Mus_musculus.GRCm38.${ENSEMBLver}.gff3
    ENSEMBLurl=ftp://ftp.ensembl.org/pub/release-${ENSEMBLver}/gff3/mus_musculus/${ENSEMBLfile}.gz
    #
    #${CURL} ${ENSEMBLurl} | ${GUNZIP} > ${DATADIR}/${ENSEMBLfile}
    logit 'prepEnsembl...'
    ${PYTHON} ${BIN}/prepEnsembl.py 2>> ${LOGFILE} < ${DATADIR}/${ENSEMBLfile} | ${SPLITCMD} -t "ensembl.chr%s.gff"
fi

########
# MERGE phase
if [ $nargs -eq 0 -o $domerge == T ]; then
    for i in "${CHRS[@]}"
    do : 
	logit
	logit "merging chr${i}..."
	${PYTHON} ${BIN}/merge.py ${WORKINGDIR}/*.chr${i}.gff > ${WORKINGDIR}/chr${i}.gff 2>> ${LOGFILE}
    done
    logit "catting files..."
    cat ${WORKINGDIR}/chr*.gff | ${PYTHON} ${BIN}/reassignIDs.py > ${WORKINGDIR}/MGI.gff3 2>> ${LOGFILE}
fi

########
# EXOME phase
if [ $nargs -eq 0 -o $doexome == T ]; then
    logit "exome phase not implemented yet"
    #logit "creating MGI exome file..."
    #${PYTHON} ${BIN}/exome.py ${WORKINGDIR}/MGI.gff3 > ${WORKINGDIR}/MGI.exome.gff3 2>> ${LOGFILE}
fi

########
logit 'Refresh finished.'
