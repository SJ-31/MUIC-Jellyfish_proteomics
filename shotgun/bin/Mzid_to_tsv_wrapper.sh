#!/usr/bin/env sh
#
mzid=$1
tsv=$2

mono $MzidToTsvConverter \
    -mzid $mzid \
    -tsv $tsv \
    -unroll \
