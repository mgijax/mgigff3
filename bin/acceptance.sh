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
T="Number of Bmp4 genes = 1"
v=`grep MGI:88180 ${FILE} | grep '	gene	' | wc -l`
try "$T" $v -eq 1

###
T="Number of features in Bmp4 model > 30"
v=`grep MGI:88180 ${FILE} | wc -l`
try "$T" $v -gt 30

###
T="Lines in file >= 3M"
v=`cat ${FILE} | wc -l`
try "$T" $v -ge 3000000

###
T="Protein coding genes > 20000"
v=`grep protein_coding_gene ${FILE} | wc -l`
try "$T" $v -ge 20000

###
#
if [ $failed -eq 0 ]; then
    logit "PASSED all tests."
else
    logit "FAILED" $failed "tests."
fi
exit $failed
