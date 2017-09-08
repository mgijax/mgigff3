#!/usr/bin/bash

source config.sh

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
    checkExit
fi

########
# MGI computed
if [ $nargs -eq 0 -o $domgicomputed == T ]; then
    logit "prepMgiComputed..."
    ${BIN}/refreshMgiComputed.sh 2>> ${LOGFILE}
    checkExit
fi

########
# NCBI
if [ $nargs -eq 0 -o $doncbi == T ]; then
    logit 'prepNcbi...'
    ${PYTHON} ${BIN}/prepNcbi.py 2>> ${LOGFILE} < ${DATADIR}/ref_GRCm38.p4_top_level.gff3 | ${SPLITCMD} -t "ncbi.chr%s.gff"
    checkExit
fi

########
# miRBase
if [ $nargs -eq 0 -o $domirbase == T ]; then
    logit 'prepMirbase...'
    ${PYTHON} ${BIN}/prepMirbase.py 2>> ${LOGFILE} < ${DATADIR}/miRBase21_mmu.gff3 | ${SPLITCMD} -t "mirbase.chr%s.gff"
    checkExit
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
    checkExit
fi

########
# MERGE phase
if [ $nargs -eq 0 -o $domerge == T ]; then
    for i in "${CHRS[@]}"
    do : 
	logit
	logit "merging chr${i}..."
	${PYTHON} ${BIN}/merge.py ${WORKINGDIR}/*.chr${i}.gff > ${WORKINGDIR}/chr${i}.gff 2>> ${LOGFILE}
	checkExit
    done
    logit "catting files..."
    echo '##gff-version 3' > ${WORKINGDIR}/MGI.gff3
    cat ${WORKINGDIR}/chr*.gff >> ${WORKINGDIR}/MGI.gff3 2>> ${LOGFILE}


    logit "Generating sample file..."
    ${BIN}/sample.sh < ${WORKINGDIR}/MGI.gff3 > ${WORKINGDIR}/MGI.sample.gff3
    checkExit

    # Generate feature type profile
    logit "Generating feature type profile..."
    ${PYTHON} ${BIN}/countPCrels.py < ${WORKINGDIR}/MGI.gff3 > ${WORKINGDIR}/MGI.counts.txt
    checkExit
fi

########
# EXOME phase
if [ $nargs -eq 0 -o $doexome == T ]; then
    #logit "creating MGI exome file..."
    logit "exome phase not yet implemented"
    #${PYTHON} ${BIN}/exome.py ${WORKINGDIR}/MGI.gff3 > ${WORKINGDIR}/MGI.exome.gff3 2>> ${LOGFILE}
    #checkExit
fi

########
logit 'Refresh finished.'
