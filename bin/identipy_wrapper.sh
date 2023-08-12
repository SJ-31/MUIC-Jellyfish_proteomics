#!/usr/bin/env sh
while getopts "g:p:c:" opt; do
    case $opt in
        g)
            config=$OPTARG ;;
        p)
            prefix=$OPTARG ;;
        c)
            converter=$OPTARG ;;
    esac
done
#
for mzML in *mzML
do
    identipy $mzML \
        -cfg ${config} \
        -at \
        -out .
done


for tsv in *tsv
do
    sample=$(echo $tsv | sed 's/\.tsv//')
    /home/shannc/anaconda3/bin/Rscript "$converter" $tsv \
    ${sample}.pin \
    ${sample}
done
