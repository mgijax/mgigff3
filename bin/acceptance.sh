#!/usr/bin/bash

source config.sh

FILE=${WORKINGDIR}/MGI.gff3

# Run some basic sanity checks on the final result.

function try {
    msg=$1
    cmd=$2
    op=$3
    val=$4
    c=(`$cmd`)
    test $c $op $val
    b=$(( 1 - $?))
    assert $b $msg
    logit "OK:" $msg
}

v=`grep MGI:88180 ${FILE} | grep '	gene	' | wc -l`
assert "Number of Bmp4 gene features = 1" $v -eq 1

v=`cat ${FILE} | wc -l`
assert "Lines in file >= 3M" $v -ge 3000000

