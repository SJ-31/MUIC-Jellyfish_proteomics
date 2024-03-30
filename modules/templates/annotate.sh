#!/usr/bin/env bash

annotate.py -i !{combined_tsv} \
    --output !{params.pref}_downloads_anno-1.tsv

eggdir="annotate_eggnog_unmatched"
ipdir="annotate_interpro_unmatched"
previous_eggnog="!{outdir}/${eggdir}"
previous_interpro="!{outdir}/${ipdir}"

if [ -e needs_annotating.fasta ]; then
    mkdir ${eggdir}
    if [ ! -d "${previous_eggnog}" ]; then
    mv -Z needs_annotating* ${eggdir}; cd ${eggdir}
    export EGGNOG_DATA_DIR=!{params.eggnog_data_dir}
    conda run -n eggnog emapper.py -i needs_annotating.fasta \
        -m mmseqs \
        --output temp \
        --report_orthologs
    elif [ -e annotation_complete.fasta ]; then
        exit 0
    else
        mv -Z "${previous_eggnog}" .
        mv -Z needs_annotating* ${eggdir}; cd ${eggdir}
        touch eggnog_copied.txt
    fi

    Rscript !{params.bin}/R/sort_eggnog.r \
        --output_fasta needs_annotating2.fasta \
        --blast needs_annotating.tsv \
        --output_unmatched eggnog_unmatched.tsv \
        --output temp.tsv \
        --seeds temp.emapper.seed_orthologs \
        --annotations temp.emapper.annotations

    eggnog_seq.py -e !{params.eggnog_db_fasta} \
        -i temp.tsv \
        -o eggnog_matched.tsv

    mv needs_annotating2.fasta ..
    # No need to get the sequences from the eggnog db, because
    # this is dealing with full-length proteins, not peptides
    cd ..

    annotate.py --merge_eggnog \
        --more_anno !{params.pref}_downloads_anno-1.tsv \
        --output !{params.pref}_downloads_anno-2.tsv \
        --eggnog_tsv ${eggdir}/eggnog_matched.tsv

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
            touch ${ipdir}/interpro_copied.txt
        fi

        annotate.py --merge_interpro \
            --interpro_query ${ipdir}/needs_annotating2.fasta \
            -r !{params.bin}/R \
            -i ${ipdir}/annotated.tsv \
            --more_anno !{params.pref}_downloads_anno-2.tsv \
            --output !{params.pref}_downloads_anno-3.tsv
        mv sorted.tsv ${ipdir}
    fi
fi

