process SPADES {
    publishDir "$outdir", mode: "copy"
    conda "/home/shannc/anaconda3/envs/spades"

    input:
    tuple val(variant), path(reads)
    val(outdir)
    //
    output:
    path "*Spades*"
    //
    exec:
    formatted = reads[0].baseName.replaceAll(/.*-/, '')
                                 .replaceAll(/[12].*/, '')
    script:
    """
    python /home/shannc/tools/SPAdes-3.14.1-Linux/bin/rnaspades.py -o ${variant}_${formatted}-Spades -1 ${reads[0]} -2 ${reads[1]}
    """
    //
}

process TRANSDECODE {
    publishDir "$outdir", mode: 'copy'

    input:
    tuple val(folder), path(assembly)
    val(outdir)
    //
    output:
    path("*")
    //
    script:
    formatted = "${assembly.baseName}_TR"
    """
    TransDecoder.LongOrfs -t $assembly \
    -m 35 \
    -O $formatted
    """
    // The minimum length has been lowered to 35 from the default of 100 despite the fact that this greatly increases false positive rate. Lowering it to this level is necessary though because venom peptides can be this small
    // Although the false positive rate increases with a lower minimum peptide
    // length, venom peptides
}
