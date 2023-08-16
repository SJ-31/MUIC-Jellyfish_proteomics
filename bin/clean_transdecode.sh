#!/bin/bash

# $1 will be the transcriptome name
cat $2 | \
    sed -e "s/>.*:\(.*\)/>"$1"_\1/" \
    -e 's/(+)/_plus-TRANSCRIPTOME/' \
    -e 's/(-)/_minus-TRANSCRIPTOME/' \
    -e 's/-/_/' \
    -e 's/\*//' \
