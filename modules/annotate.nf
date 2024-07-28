process ANNOTATE {
    publishDir "$outdir", mode: "copy", pattern: "{*.tsv,*.fasta}"
    publishDir "$outdir", mode: "copy", pattern: "annotate_eggnog_unmatched"
    publishDir "$outdir", mode: "copy", pattern: "annotate_interpro_unmatched"

    input:
    path(combined_tsv)
    val(outdir)
    //

    output:
    path("${params.pref}_downloads_anno-complete.tsv"), emit: annotations
    path("${params.pref}_downloads_anno*")
    path("annotate_eggnog_unmatched"), optional: true
    path("annotation_complete.fasta"), emit: complete, optional: true
    path("annotate_interpro_unmatched"), optional: true
    path("still_unannotated.fasta"), emit: unannotated, optional: true
    path("*fasta"), optional: true
    //

    shell:
    check = file("${outdir}/${params.pref}_downloads_anno-complete.tsv")
    if (check.exists()) {
        '''
        cp -r !{outdir}/* .
        '''
    } else {
        template 'annotate.sh'
    }
}

process ANNOTATE_1 {
    publishDir "$outdir", mode: "copy"

    input:
    path(combined_tsv)
    val(outdir)
    //

    output:
    tuple path("${params.pref}_downloads_anno-*.tsv"), \
        path("*.fasta"), emit: to_annotate
    //

    shell:
    completed = file("${outdir}/${params.pref}_downloads_anno-*.tsv")
    if (completed.exists()) {
        '''
        cp -r !{outdir}/* .
        '''
    } else {
        '''
        annotate.py -i !{combined_tsv} \
            --output !{params.pref}_downloads_anno-1.tsv
        if [ -e annotation_complete.fasta ]; then
            mv !{params.pref}_downloads_anno-1.tsv !{params.pref}_downloads_anno-complete.tsv
            exit 0
        fi
        '''
    }
}

process ANNOTATE_2 {
    publishDir "$outdir", mode: "copy"

    input:
    tuple path(combined_tsv), path(fastas)
    val(outdir)
    //

    output:
    tuple path("${params.pref}_downloads_anno-*.tsv"), \
        path("*.fasta"), emit: to_annotate
    path("annotate_eggnog_unmatched"), optional: true
    //

    shell:
    check = file("${outdir}/annotate_eggnog_unmatched")
    if (check.exists()) {
        '''
        cp -r !{check} .
        cp !{combined_tsv} !{combined_tsv}
        cp !{outdir}/*.fasta .
        '''
    } else {
    '''
    eggdir="annotate_eggnog_unmatched"
    mkdir ${eggdir}
    mv -Z needs_annotating* ${eggdir}; cd ${eggdir}
    export EGGNOG_DATA_DIR=!{params.eggnog_data_dir}
    conda run -n eggnog emapper.py -i needs_annotating.fasta \
        -m mmseqs \
        --output temp \
        --report_orthologs

    Rscript !{params.bin}/R/sort_eggnog.r \
        --output_fasta needs_annotating2.fasta \
        --blast needs_annotating.tsv \
        --output_unmatched eggnog_unmatched.tsv \
        --output eggnog_matched.tsv \
        --seeds temp.emapper.seed_orthologs \
        --annotations temp.emapper.annotations

    mv needs_annotating2.fasta ..
    cd ..

    annotate.py --merge_eggnog \
        --more_anno !{params.pref}_downloads_anno-1.tsv \
        --output !{params.pref}_downloads_anno-2.tsv \
        --eggnog_tsv ${eggdir}/eggnog_matched.tsv

    if [ -e annotation_complete.fasta ]; then
        mv !{params.pref}_downloads_anno-2.tsv !{params.pref}_downloads_anno-complete.tsv
        exit 0
    fi
    '''
    }
}

process ANNOTATE_3 {
    publishDir "$outdir", mode: "copy"

    input:
    tuple path(combined_tsv), path(fastas)
    val(outdir)
    //

    output:
    tuple path("${params.pref}_downloads_anno-*.tsv"), \
        path("*.fasta"), emit: to_annotate
    path("annotate_interpro_unmatched"), optional: true
    //

    shell:
    check = file("${outdir}/annotate_interpro_unmatched")
    if (check.exists()) {
        '''
        cp -r !{check} .
        cp !{combined_tsv} !{combined_tsv}
        cp !{outdir}/*.fasta .
        '''
    } else {
    '''
    ipdir="annotate_interpro_unmatched"
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

    annotate.py --merge_interpro \
        --interpro_query ${ipdir}/needs_annotating2.fasta \
        -r !{params.bin}/R \
        -i ${ipdir}/annotated.tsv \
        --more_anno !{params.pref}_downloads_anno-2.tsv \
        --output !{params.pref}_downloads_anno-3.tsv
    mv sorted.tsv ${ipdir}
    mv !{params.pref}_downloads_anno-3.tsv !{params.pref}_downloads_anno-complete.tsv
    '''
    }
}
