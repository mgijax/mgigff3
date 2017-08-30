#!/usr/bin/bash

PYTHON=python
CURL=curl
CP=cp
GUNZIP=gunzip
GREP=grep
SORT=sort
MKDIR=mkdir

# https://stackoverflow.com/questions/59895/getting-the-source-directory-of-a-bash-script-from-within
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
PYTHONPATH=${DIR}/lib:${PYTHONPATH:-.}

#
BIN=${DIR}
WORKINGDIR=${DIR}/../work
${MKDIR} -p ${WORKINGDIR}
CACHEDIR=${DIR}/../cache
${MKDIR} -p ${CACHEDIR}

DATADIR=${DIR}/../data
#DATADIR=/data/research/mouse_build_38_external/original_annotations

SORTCMD='sort -k 1,1 -k 4n,4n -k 5nr,5nr'
SPLIT="${PYTHON} ${BIN}/splitGff.py -d ${WORKINGDIR}"

LOGFILE=${WORKINGDIR}/LOG
truncate -s 0 ${LOGFILE}

function logit {
    echo `date` $1 >> ${LOGFILE}
}

export DIR BIN PYTHONPATH WORKINGDIR CACHEDIR DATADIR LOGFILE SORT GREP SORTCMD SPLIT

nargs=$#

domgi=F
domgicomputed=F
doncbi=F
doensembl=F
domirbase=F
domerge=F

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
    *)
        echo "Unrecognized option:" $1
	exit -1
    esac
    shift
done


########

logit "Starting refresh..." 

########
# MGI
if [ $nargs -eq 0 -o $domgi == T ]; then
    logit 'prepMgi...'
    ${PYTHON} ${BIN}/prepMgi.py | ${SORTCMD} | ${SPLIT} -t "mgi.chr%s.gff" 
fi

########
# MGI computed
if [ $nargs -eq 0 -o $domgicomputed == T ]; then
    logit "prepMgiComputed..."
    prepMgiComputed.sh
fi

########
# NCBI
if [ $nargs -eq 0 -o $doncbi == T ]; then
    logit 'prepNcbi...'
    ${PYTHON} ${BIN}/prepNcbi.py < ${DATADIR}/ref_GRCm38.p4_top_level.gff3 | ${SPLIT} -t "ncbi.chr%s.gff"
fi

########
# miRBase
if [ $nargs -eq 0 -o $domirbase == T ]; then
    logit 'prepMirbase...'
    ${PYTHON} ${BIN}/prepMirbase.py < ${DATADIR}/miRBase21_mmu.gff3 | ${SPLIT} -t "mirbase.chr%s.gff"
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
    ${PYTHON} ${BIN}/prepEnsembl.py < ${DATADIR}/${ENSEMBLfile} | ${SPLIT} -t "ensembl.chr%s.gff"
fi

########
# MERGE phase
if [ $nargs -eq 0 -o $domerge == T ]; then
    chrs=( 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 X Y MT )
    for i in "${chrs[@]}"
    do : 
	logit
	logit "merging chr${i}..."
	${PYTHON} ${BIN}/merge.py ${WORKINGDIR}/*.chr${i}.gff > ${WORKINGDIR}/chr${i}.gff 2>> ${LOGFILE}
    done
fi

########
logit 'Refresh finished.'
