#!/usr/bin/bash

source config.sh

FILE=${WORKINGDIR}/MGI.gff3

# Run some basic sanity checks on the final result.

failed=0
function try {
    assert "$1" "$2" "$3" "$4"
    failed=$(( failed + $? ))
}

###
T="Number of lines in file"
v=`cat ${FILE} | wc -l`
try "$T" $v -ge 3000000

###
T="Number of protein coding genes"
v=`grep protein_coding_gene ${FILE} | wc -l`
try "$T" $v -gt 20000

###
T="Number of pseudogenes"
v=`grep "	pseudogene	" ${FILE} | wc -l`
try "$T" $v -gt 14000

###
T="Number of miRNA products"
v=`grep "	miRNA	" ${FILE} | wc -l`
try "$T" $v -gt 2500

###
T="Number of Bmp4 genes"
v=`grep MGI:88180 ${FILE} | grep '	gene	' | wc -l`
try "$T" $v -eq 1

###
T="Number of features in Bmp4 model"
v=`grep MGI:88180 ${FILE} | wc -l`
try "$T" $v -gt 30

###
#
if [ $failed -eq 0 ]; then
    logit "PASSED all tests."
    exit 0
else
    logit "FAILED" $failed "tests."
    exit $failed
fi
