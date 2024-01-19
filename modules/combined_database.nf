process COMBINED_DATABASE {
    publishDir "$outdir", mode: "copy"
    publishDir "$text_loc", mode: "copy", pattern: "*.DB.txt"

    input:
    path(other_fasta)
    path(denovo_peptides)
    val(outdir)
    val(text_loc)
    //

    output:
    path("decoysWnormal.fasta")
    path("all_normal.fasta")
    path("interpro2go")
    path("pfam2go")
    path("pfam_entries.tsv")
    path("all_decoys.fasta")
    path("blast_db")
    path("seq-header_mappings.tsv")
    path("${params.pref}.DB.txt"), emit: listing
    //

    shell:
    '''
    outdir=!{outdir}
    cat !{other_fasta} !{denovo_peptides} | seqkit rmdup > intermediate.fasta
    create_decoys.py intermediate.fasta seq-header_mappings.tsv
    rm intermediate.fasta

    find . \\( -name "*fasta" -o -name "*tsv" \\) -type f > temp_list.txt
    cat temp_list.txt | sed 's;./;!{outdir}/;' > !{params.pref}.DB.txt

    mv downloaded.fasta temp; cat temp | seqkit rmdup > downloaded.fasta
    makeblastdb -in downloaded.fasta \
        -parse_seqids \
        -blastdb_version 5 \
        -title Cnidarian_proteins \
        -dbtype prot \
        -out downloaded
    mkdir blast_db; mv downloaded* blast_db

    interpro_api.py --get_pfam pfam_entries.tsv
    curl http://current.geneontology.org/ontology/external2go/pfam2go > pfam2go
    curl http://current.geneontology.org/ontology/external2go/interpro2go > interpro2go
    '''
    //
}
