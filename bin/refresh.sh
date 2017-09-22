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
dodistrib=F
doagr=F

until [ -z "$1" ]  # Until all parameters used up . . .
do
    case "$1" in
    mgi)
	domgi=T
        ;;
    mgicom*)
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
    dist*)
	dodistrib=T
        ;;
    agr)
	doagr=T
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
    # Queries MGI data. Writes genes and pseudogenes into gff3 file. Each has Dbxrefs containing
    # any provider model ids. Downstream merging will be directed by this data.
    logit 'prepMgi...'
    ${PYTHON} ${BIN}/prepMgi.py 2>> ${LOGFILE} | ${SORTCMD} | ${SPLITCMD} -t "mgi.chr%s.gff" 
    checkExit

    logit 'Counting mgi after prep...'
    ${COUNTCMD} ${WORKINGDIR}/mgi.chr*.gff > ${WORKINGDIR}/mgi.counts.txt
    checkExit

fi

########
# MGI computed
if [ $nargs -eq 0 -o $domgicomputed == T ]; then
    # Computes Blat'ed models for MGI genes that don't have models but do have good enough
    # sequences. The sequences are blatted against the mouse genome, and the best scoring
    # hits are made into models. 
    logit "prepMgiComputed..."
    ${BIN}/refreshMgiComputed.sh 2>> ${LOGFILE}  | ${SPLITCMD} -t "mgicomputed.chr%s.gff"
    checkExit

    logit 'Counting mgicomputed after prep...'
    ${COUNTCMD} ${WORKINGDIR}/mgicomputed.chr*.gff > ${WORKINGDIR}/mgicomputed.counts.txt
    checkExit

fi

########
# NCBI
if [ $nargs -eq 0 -o $doncbi == T ]; then
    #
    NCBIfile=ref_GRCm38.p4_top_level.gff3
    NCBIurl=ftp://ftp.ncbi.nlm.nih.gov/genomes/Mus_musculus/GFF/${NCBIfile}.gz
    #
    logit 'Downloading ncbi'
    ${CURL} ${NCBIurl} | ${GUNZIP} > ${DATADIR}/${NCBIfile}
    checkExit

    logit 'Counting ncbi downloaded...'
    ${COUNTCMD} ${DATADIR}/${NCBIfile} > ${DATADIR}/ncbi.counts.txt
    checkExit

    logit 'Prep ncbi...'
    ${PYTHON} ${BIN}/prepNcbi.py 2>> ${LOGFILE} < ${DATADIR}/${NCBIfile} | ${SPLITCMD} -t "ncbi.chr%s.gff"
    checkExit

    logit 'Counting ncbi after prep...'
    ${COUNTCMD} ${WORKINGDIR}/ncbi.chr*.gff > ${WORKINGDIR}/ncbi.counts.txt
    checkExit

fi

########
# miRBase
if [ $nargs -eq 0 -o $domirbase == T ]; then
    #
    MIRfile=mmu.gff3
    MIRurl=ftp://mirbase.org/pub/mirbase/CURRENT/genomes/${MIRfile}
    #
    logit 'Downloading mirbase'
    ${CURL} ${MIRurl} > ${DATADIR}/${MIRfile}
    checkExit

    logit 'Counting mirbase downloaded...'
    ${COUNTCMD} ${DATADIR}/${MIRfile} > ${DATADIR}/mirbase.counts.txt
    checkExit

    logit 'prepMirbase...'
    ${PYTHON} ${BIN}/prepMirbase.py 2>> ${LOGFILE} < ${DATADIR}/${MIRfile} | ${SPLITCMD} -t "mirbase.chr%s.gff"
    checkExit

    logit 'Counting mirbase after prep...'
    ${COUNTCMD} ${WORKINGDIR}/mirbase.chr*.gff > ${WORKINGDIR}/mirbase.counts.txt
    checkExit
fi

########
# ENSEMBL
if [ $nargs -eq 0 -o $doensembl == T ]; then
    #
    # ftp://ftp.ensembl.org/pub/release-90/gff3/mus_musculus/Mus_musculus.GRCm38.90.gff3.gz
    # Note that the ensembl version number increases each release.
    #
    ENSEMBLver=90
    ENSEMBLfile=Mus_musculus.GRCm38.${ENSEMBLver}.gff3
    ENSEMBLurl=ftp://ftp.ensembl.org/pub/release-${ENSEMBLver}/gff3/mus_musculus/${ENSEMBLfile}.gz
    #
    logit 'Downloading ensembl...'
    ${CURL} ${ENSEMBLurl} | ${GUNZIP} > ${DATADIR}/${ENSEMBLfile}
    checkExit

    logit 'Counting ensembl downloaded...'
    ${COUNTCMD} ${DATADIR}/${ENSEMBLfile} > ${DATADIR}/ensembl.counts.txt
    checkExit

    logit 'Prep ensembl...'
    ${PYTHON} ${BIN}/prepEnsembl.py 2>> ${LOGFILE} < ${DATADIR}/${ENSEMBLfile} | ${SPLITCMD} -t "ensembl.chr%s.gff"
    checkExit

    logit 'Counting ensembl after prep...'
    ${COUNTCMD} ${WORKINGDIR}/ensembl.chr*.gff > ${WORKINGDIR}/ensembl.counts.txt
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

    logit "Running acceptance tests..."
    ${BIN}/acceptance.sh
    checkExit

    logit "Generating sample file..."
    ${BIN}/sample.sh < ${WORKINGDIR}/MGI.gff3 > ${WORKINGDIR}/MGI.sample.gff3
    checkExit

    # Generate feature type profile
    logit "Generating feature type profile..."
    ${COUNTCMD} < ${WORKINGDIR}/MGI.gff3 > ${WORKINGDIR}/MGI.counts.txt
    checkExit
fi

########
# EXOME phase
if [ $nargs -eq 0 -o $doexome == T ]; then
    logit "creating MGI exome file..."
    ${PYTHON} ${BIN}/exome.py < ${WORKINGDIR}/MGI.gff3 > ${WORKINGDIR}/MGI.exome.gff3 2>> ${LOGFILE}
    checkExit
fi

########
# AGR phase
if [ $nargs -eq 0 -o $doagr == T ]; then
    logit "Generating GFF3 file for AGR..."
    ${PYTHON} ${BIN}/trimForAgr.py < ${WORKINGDIR}/MGI.gff3 > ${WORKINGDIR}/MGI.agr.gff3 2>> ${LOGFILE}
    checkExit
fi

########
# DISTRIB phase
if [ $nargs -eq 0 -o $dodistrib == T ]; then
    function distrib {
	filename=$(basename "$1")
	extension="${filename##*.}"
	filenameNoExt="${filename%.*}"
	nfn=${DISTRIBDIR}/archive/${filenameNoExt}.${DATESTAMP2}.${extension}
	logit "Copying $1 to ${nfn} ..."
	${CP} $1 ${nfn}
	checkExit

	symLinkName=${DISTRIBDIR}/${filename}
	${RM} -f ${symLinkName}
	${LN} -s ${nfn} ${symLinkName}
    }
    distrib ${WORKINGDIR}/MGI.gff3
    distrib ${WORKINGDIR}/MGI.agr.gff3 
    distrib ${WORKINGDIR}/MGI.exome.gff3
    checkExit
fi

########
logit 'Refresh finished. No errors detected.'
