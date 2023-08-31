#!/usr/bin/env sh
# Only run in the smsnet work directory, with the EXTRACT_SMSNET process
grep "<s>" -v smsnet_temp.txt | grep -v "<unk>" | grep -v "</s>" | \
    sed -e 's/ //g' -e 's/m//g' | \
    awk 'length >= 7 {printf ">SMSNet%s-DENOVO\n%s\n",FNR,$0}' \
    > smsnet_normal.fasta
