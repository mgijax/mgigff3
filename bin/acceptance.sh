#!/usr/bin/bash
#
# WARNING: some of the tests below include literal TAB characters.
# Sometimes the TABs get converted to spaces, which messes up the tests. Something to watch out for,

source config.sh

FILE=${WORKINGDIR}/MGI.gff3

# Run some basic sanity checks on the final result.

failed=0
function try {
    assert "$1" "$2" "$3" "$4"
    failed=$(( failed + $? ))
}

###
try \
  "Number of lines in file" \
  `cat ${FILE} | wc -l` \
  -ge \
  3000000

###
try \
  "Number of protein coding genes" \
  `grep "so_term_name=protein_coding_gene" ${FILE} | wc -l` \
  -gt \
  20000

###
try \
  "Number of pseudogenes" \
  `grep "	pseudogene	" ${FILE} | wc -l` \
  -gt \
  13000

###
try \
    "Number of miRNA products" \
    `grep "	miRNA	" ${FILE} | wc -l` \
    -gt \
    2500

###
try \
    "Number of Bmp4 genes" \
    `grep MGI:88180 ${FILE} | grep "	gene	" | wc -l` \
    -eq \
    1

###
try \
    "Make sure B230334L07Rik only occurs once" \
    `grep MGI:2443922 ${FILE} | grep "	gene	" | wc -l` \
    -eq \
    1

###
try \
    "Number of features in Bmp4 model" \
    `grep MGI:88180 ${FILE} | wc -l` \
    -gt \
    30

###
#
if [ $failed -eq 0 ]; then
    logit "PASSED all tests."
    exit 0
else
    logit "FAILED" $failed "tests."
    exit $failed
fi
