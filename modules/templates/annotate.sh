#!/usr/bin/env bash

annotate.py -i !{combined_tsv} \
    -m download_meta-!{name}.tsv \
    -a download_anno-!{name}.tsv

eggdir="!{name}_eggnog-unmatched"
ipdir="!{name}_interpro-unmatched"
previous_eggnog="!{outdir}/${eggdir}"
previous_interpro="!{outdir}/${ipdir}"

if [ -e needs_annotating.fasta ]; then
    if [ ! -d "${previous_eggnog}" ]; then
    mkdir ${eggdir}; mv needs_annotating* ${eggdir}; cd ${eggdir}
    export EGGNOG_DATA_DIR=!{params.eggnog_data_dir}
    conda run -n eggnog emapper.py -i needs_annotating.fasta \
        --output temp \
        --report_orthologs
    else
        cp -r "${previous_eggnog}" .
        cd ${eggdir}
    fi

    Rscript !{params.bin}/sort_eggnog.r \
        --output_fasta needs_annotating2.fasta \
        --blast needs_annotating.tsv \
        --output_unmatched !{name}_eggnog_unmatched.tsv \
        --output_meta !{name}_meta-eggnog.tsv \
        --output_anno !{name}_anno-eggnog.tsv \
        --seeds temp.emapper.seed_orthologs \
        --annotations temp.emapper.annotations
    mv needs_annotating2.fasta ..
    # No need to get the sequences from the eggnog db, because
    # this is dealing with full-length proteins, not peptides
    cd ..

    annotate.py --merge_eggnog \
        --anno_tsv download_anno-!{name}.tsv \
        --eggnog_anno_tsv ${eggdir}/!{name}_anno-eggnog.tsv

    if [ -e needs_annotating2.fasta ]; then
        if [ ! -d "${previous_interpro}" ]; then
            interpro_header='query\tsequence_md5\tlength\tmember_db\tdb_accession\tdescription\tstart\tstop\tevalue\tstatus\tdate\tinterpro_accession\tinterpro_description\tGO\tpathways'
            mkdir ${ipdir}; mv needs_annotating2.fasta ${ipdir}; cd ${ipdir}
            interproscan.sh -i needs_annotating2.fasta \
                -f tsv \
                -goterms \
                -pa \
                -b temp
            echo -e "$interpro_header" > header.txt
            cat header.txt temp.tsv > annotated.tsv
            rm temp.tsv
            cd ..
        else
            cp -r "${previous_interpro}" .
        fi

        annotate.py --merge_interpro \
            --interpro_query ${ipdir}/needs_annotating2.fasta \
            -r !{params.bin} \
            -i ${ipdir}/annotated.tsv \
            --anno_tsv download_anno-!{name}.tsv
        mv sorted.tsv ${ipdir}
    fi
fi

