#!/usr/bin/bash

set -o pipefail

# ---------------------
PYTHON=python
PYTHON24="python2.4"
CURL=curl
CP=cp
GUNZIP=gunzip
GREP=grep
SORT=sort
MKDIR=mkdir
DATE=date
TOUCH=touch
CUT=cut
SED=sed
GFCLIENT=gfClient
PSLREPS=pslReps
GREP=grep

# ---------------------
BLAT_HOST="bhmgiapp01.jax.org"
BLAT_PORT="9038"

# The chromosomes and their order.
CHRS=( 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 X Y MT )

# ---------------------
# https://stackoverflow.com/questions/59895/getting-the-source-directory-of-a-bash-script-from-within
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
BIN=${DIR}
DATADIR=${DIR}/../data
#DATADIR=/data/research/mouse_build_38_external/original_annotations
WORKINGDIR=${DIR}/../work
DISTRIBDIR=${DIR}/../dist

PYTHONPATH=${DIR}/lib:${PYTHONPATH:-.}

#
DATESTAMP=`${DATE} +"%Y-%m-%d"`

# ---------------------
SORTCMD="sort -k 1,1 -k 4n,4n -k 5nr,5nr"
SPLITCMD="${PYTHON} ${BIN}/splitGff.py -d ${WORKINGDIR}"

# ---------------------
${MKDIR} -p ${WORKINGDIR}

# ---------------------
LOGFILE=${WORKINGDIR}/LOG.${DATESTAMP}
${TOUCH} ${LOGFILE}

# ---------------------
# Echos its arguments to the log file. Prepends a datetime stamp.
#
function logit {
    echo `date` "$*" >> ${LOGFILE}
}

# ---------------------
# Logs a message and exits with error code 1.
#
function die {
    logit "$*"
    exit 1
}

# ---------------------
# Tests an assertion. 
# If success, logs OK and returns 0.
# If failure, logs FAILED, and return 1.
# Arguments:
#    label	$1 = label to print for test (e.g., title)
#    val1	$2 = the computed valued to test (e.g., the count of protein coding genes)
#    op		$3 = the operator to use. One of: -eq -ne -gt -ge -lt -le
#    val2	$4 = the value to test against, e.g., a sanity threshhold value for number of protein coding genes
#
function assert {
    test $2 $3 $4
    if [ $? -ne 0 ]; then
        logit "FAILED ASSERTION: $1: $2 $3 $4"
	return 1
    else
	logit "OK: $1: $2 $3 $4"
	return 0
    fi
}

# ---------------------
# If the exit code ($?) is not zero, exits with a message.
#
function checkExit {
    c=$?
    if [ $c -ne 0 ]; then
        die "ERROR: Caught error exit code." 
    fi
    return 0
}

# ---------------------
export DIR BIN PYTHONPATH WORKINGDIR DATADIR LOGFILE SORT GREP SORTCMD SPLITCMD
export DISTRIBDIR
export CP TOUCH CUT SED GFCLIENT PSLREPS GREP BLAT_HOST  BLAT_HOST

