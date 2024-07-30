#!/usr/bin/env bash

if [[ ! -d !{outdir} ]]; then
    mkdir !{outdir}
fi

if [ -e !{outdir}/!{params.pref}_downloads_anno-1.tsv ]; then
    cp !{outdir}/!{params.pref}_downloads_anno-1.tsv .
    cp !{outdir}/needs_annotating* .
else
    annotate.py -i !{combined_tsv} \
        --output !{params.pref}_downloads_anno-1.tsv \
        --previous_saved !{params.saved_dir}/annotations.tsv
fi


if [ -e annotation_complete.fasta ]; then
    mv !{params.pref}_downloads_anno-1.tsv !{params.pref}_downloads_anno-complete.tsv
    exit 0
else
    cp !{params.pref}_downloads_anno-1.tsv !{outdir}
    cp needs_annotating* !{outdir}
fi

eggdir="annotate_eggnog_unmatched"
ipdir="annotate_interpro_unmatched"
previous_eggnog="!{outdir}/${eggdir}"
previous_interpro="!{outdir}/${ipdir}"

mkdir ${eggdir}
if [ ! -d "${previous_eggnog}" ]; then
    mv -Z needs_annotating* ${eggdir}; cd ${eggdir}
    export EGGNOG_DATA_DIR=!{params.eggnog_data_dir}
    conda run -n eggnog emapper.py -i needs_annotating.fasta \
        -m mmseqs \
        --output temp \
        --report_orthologs
else
    cp -r "${previous_eggnog}" .
    cd ${eggdir}
    touch eggnog_copied.txt
fi

Rscript !{params.bin}/R/sort_eggnog.r \
    --output_fasta needs_annotating2.fasta \
    --blast needs_annotating.tsv \
    --output_unmatched eggnog_unmatched.tsv \
    --output eggnog_matched.tsv \
    --seeds temp.emapper.seed_orthologs \
    --annotations temp.emapper.annotations

mv needs_annotating2.fasta ..
# No need to get the sequences from the eggnog db, because
# this is dealing with full-length proteins, not peptides
cd ..
cp -r ${eggdir} ${outdir}

annotate.py --merge_eggnog \
    --more_anno !{params.pref}_downloads_anno-1.tsv \
    --output !{params.pref}_downloads_anno-2.tsv \
    --eggnog_tsv ${eggdir}/eggnog_matched.tsv

if [ -e annotation_complete.fasta ]; then
    mv !{params.pref}_downloads_anno-2.tsv !{params.pref}_downloads_anno-complete.tsv
    exit 0
fi
# ----------------------------

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
    cp -r ${ipdir} ${outdir}
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
mv !{params.pref}_downloads_anno-3.tsv !{params.pref}_downloads_anno-complete.tsv
