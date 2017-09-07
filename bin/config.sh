#!/usr/bin/bash

set -o pipefail

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

function checkExit {
    if [ $? -ne 0 ]
    then
	logit "Caught error code. Exiting..."
	exit 1
    fi
}

export DIR BIN PYTHONPATH WORKINGDIR DATADIR LOGFILE SORT GREP SORTCMD SPLITCMD

