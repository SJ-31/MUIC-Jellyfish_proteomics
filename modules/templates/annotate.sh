#!/usr/bin/env bash

annotate.py -i !{combined_tsv} \
    -m !{name}_meta.tsv \
    -a !{name}_anno.tsv

eggdir="!{name}_eggnog-unmatched"
ipdir="!{name}_interpro-unmatched"

if [ -e needs_annotating.fasta ]; then
    mkdir ${eggdir}; mv needs_annotating* ${eggdir}; cd ${eggdir}
    export EGGNOG_DATA_DIR="!{params}.eggnog_data_dir"
    emapper.py -i needs_annotating.fasta \
        --output temp \
        --report_orthologs

    Rscript !{params.bin}/sort_eggnog.r \
        -f needs_annotating.fasta \
        --blast needs_annotating.tsv \
        --output_unmatched !{name}_eggnog_unmatched.tsv \
        --output_meta !{name}_meta-eggnog.tsv \
        --output_anno !{name}_anno-eggnog.tsv \
        -seeds temp.emapper.seed_orthologs \
        -annotations temp.emapper.annotations
    # No need to get the sequences from the eggnog db, because
    # this is dealing with full-length proteins, not peptides
    cd ..

    annotate.py --merge_eggnog \
        --anno_tsv !{name}_anno.tsv \
        --eggnog_anno_tsv ${eggdir}/ !{name}_anno-eggnog.tsv

    if [ -e needs_annotating2.fasta ]; then
        interpro_header='query\tsequence_md5\tlength\tmember_db\tdb_accession\tdescription\tstart\tstop\tevalue\tstatus\tdate\tinterpro_accession\tinterpro_description\tGO\tpathways'
        mkdir ${ipdir}; mv needs_annotating2.fasta ${ipdir}; cd ${ipdir}
        interproscan.sh -i needs_2.fasta \
            -f tsv \
            -goterms \
            -pa \
            -b temp
        echo -e "$header" > header.txt
        cat header.txt temp.tsv > annotated.tsv
        rm temp.tsv
        cd ..

        annotate.py --merge_interpro \
            --interpro_query ${ipdir}/needs_annotating2.fasta \
            -r $params.bin \
            -i ${ipdir}/annotated.tsv \
            --anno_tsv !{name}_anno.tsv
        mv sorted.tsv ${ipdir}
    fi
fi

