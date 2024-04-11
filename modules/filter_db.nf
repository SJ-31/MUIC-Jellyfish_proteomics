process FILTER_DB {
    publishDir "$params.logdir", mode: "copy", pattern: "*.log"

    input:
    path(first_pass_tsv)
    val(outdir)
    //

    output:
    val(path_to_db), emit: path
    path("blast_db")
    path("*.log")
    //

    script:
    path_to_db = "$params.outdir/blast_db/filtered"
    def check = file("$params.logdir/filter_db.log")
    if (check.exists()) {
        """
        cp -r "$outdir/blast_db" .
        cp filter_db.log
        """
    } else {
        """
        get_fasta.py -i $first_pass_tsv \
            -s seq \
            -d ProteinId \
            -o first_pass.fasta
        makeblastdb -in first_pass.fasta \
            -parse_seqids \
            -blastdb_version 5 \
            -title First_pass \
            -dbtype prot \
            -out filtered
        mkdir blast_db; mv filtered* blast_db

        cp -r blast_db $outdir
        cp .command.out filter_db.log
        """
    }

    //
}
