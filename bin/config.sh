#!/usr/bin/bash

set -o pipefail

# ---------------------
export PATH="${PATH}:/opt/python3.7/bin"
# ---------------------
export PYTHON=python3
export CURL=curl
export CP=cp
export GUNZIP=gunzip
export GZIP=gzip
export GREP=grep
export SORT=sort
export MKDIR=mkdir
export DATE=date
export TOUCH=touch
export CUT=cut
export SED=sed
export GFCLIENT=gfClient
export PSLREPS=pslReps
export GREP=grep
export RM=rm
export LN=ln
export MV=mv
export FIND=find

# ---------------------
#
export DATESTAMP=`${DATE} +"%Y-%m-%d"`
export DATESTAMP2=`${DATE} +"%Y%m%d.%H%M%S"`
export YEAR=`${DATE} +"%Y"`
export YEARMONTH=`${DATE} +"%Y%m"`

# ---------------------
export BLAT_HOST="bhmgiapp01.jax.org"

export BLAT_PORT="9038"
export BLAT_ARGS="-nohead -minIdentity=98"
#
export PSLREPS_ARGS="-singleHit -nohead "
#
# The chromosomes and their order.
export CHRS=( 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 X Y MT )

# ---------------------
export BIN=${DIR}/bin
if [ "${DATADIR}" == "" ]; then
    export DATADIR=${DIR}/data
fi
if [ "${WORKINGDIR}" == "" ]; then
    export WORKINGDIR=${DIR}/work
fi
if [ "${DISTRIBDIR}" == "" ]; then
    export DISTRIBDIR=${DIR}/dist
fi
export ARCHIVEDIR=${DISTRIBDIR}/archive
export MONTHLYDIR=${ARCHIVEDIR}/monthly
export ANNUALDIR=${ARCHIVEDIR}/annual

# Age limit in days. Archived files older than this are culled to one per year.
export ARCHIVEAGELIMIT=365

export PYTHONPATH=${PYTHONPATH:-.}:${BIN}/lib:${BIN}/lib/intermine-1.09.09-py2.7.egg
# ---------------------
#export MGI_HOST="mgi-adhoc.jax.org"
export MGI_HOST="bhmgidevdb01.jax.org"
export MGI_DATABASE="prod"
export MGI_USER="mgd_public"
export MGI_PASSWORD="mgdpub"    # FIXME: should be name of file containing password

# ---------------------
export SORTCMD="sort -k 1,1 -k 4n,4n -k 5nr,5nr"
export SPLITCMD="${PYTHON} ${BIN}/splitGff.py -d ${WORKINGDIR}"
export COUNTCMD="${PYTHON} ${BIN}/profileGff.py"

# ---------------------
#
export NCBIver=108.20200622
export NCBIbuild=GCF_000001635.26_GRCm38.p6
export NCBIfile=${NCBIbuild}_genomic.gff
export NCBIurl=https://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Mus_musculus/annotation_releases/${NCBIver}/${NCBIbuild}/${NCBIfile}.gz
export NCBIprep="${PYTHON} ${BIN}/ncbiPrep.py"
#
export MIRfile=mmu.gff3
export MIRurl=ftp://mirbase.org/pub/mirbase/CURRENT/genomes/${MIRfile}
export MIRprep="${PYTHON} ${BIN}/mirbasePrep.py"
#
export ENSEMBLver=99
export ENSEMBLbuild=GRCm38
export ENSEMBLfile=Mus_musculus.${ENSEMBLbuild}.${ENSEMBLver}.gff3
export ENSEMBLurl=ftp://ftp.ensembl.org/pub/release-${ENSEMBLver}/gff3/mus_musculus/${ENSEMBLfile}.gz
export ENSEMBLprep="${PYTHON} ${BIN}/ensemblPrep.py"
#
export LOGFILE=${WORKINGDIR}/refresh.log
#
export SO_TERM_FILE="${WORKINGDIR}/so_terms.tsv"

# ---------------------
# Echos its arguments to the log file. Prepends a datetime stamp.
#
function logit {
    echo `${DATE}` "$*" >> ${LOGFILE}
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

