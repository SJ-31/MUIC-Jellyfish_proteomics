#!/bin/bash
#
# Sed command for extracting the IDs from UniProt and Toxprot fasta files
cat $1 | \
    sed 's/>.*|\(.*\)|.*/>\1/'
