#!/usr/bin/env bash

file="$1"

echo "Number of decoys: $(grep 'rev_' $file | wc -l)"
echo "Number of targets: $(grep '>' $file | grep -v 'rev_' | wc -l)"
echo "Total sequences: $(grep '>' $file | wc -l)"
echo ""

echo "Decoy characteristics"
seqkit grep  -p "rev_*" $file -n -r | seqkit stat
echo ""

echo "Target characteristics"
seqkit grep  -p "rev_*" $file -n -r -v | seqkit stat
