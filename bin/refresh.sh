#!/usr/bin/bash

PYTHON=python
CURL=curl
CP=cp
GUNZIP=gunzip
GREP=grep

# https://stackoverflow.com/questions/59895/getting-the-source-directory-of-a-bash-script-from-within
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

BIN=${DIR}
WORKINGDIR=${DIR}/../output

DATADIR=${DIR}/../data
#DATADIR=/data/research/mouse_build_38_external/original_annotations

########
# NCBI
#${PYTHON} ${BIN}/convertNCBI.py ${DATADIR}/ref_GRCm38.p4_top_level.gff3 > ${WORKINGDIR}/ncbi.gff3

########
# miRBase
${PYTHON} ${BIN}/convertMirbase.py < ${DATADIR}/miRBase21_mmu.gff3 > ${WORKINGDIR}/mirbase.gff3

########
# ENSEMBL
ENSEMBLver=89
ENSEMBLfile=Mus_musculus.GRCm38.${ENSEMBLver}.gff3
ENSEMBLurl=ftp://ftp.ensembl.org/pub/release-${ENSEMBLver}/gff3/mus_musculus/${ENSEMBLfile}.gz
#
#${CURL} ${ENSEMBLurl} | ${GUNZIP} > ${DATADIR}/${ENSEMBLfile}
#${GREP} -v biological_region ${DATADIR}/${ENSEMBLfile} > ${WORKINGDIR}/ensembl.gff3

########
