process BLOOM {
    publishDir "$params.assembled/rnabloom/"

    input:
    tuple val(variant), path(reads)
    //
    output:
    path "*Bloom*"
    //
    exec:
    formatted = reads[0].baseName.replaceAll(/.*-/, '')
                                 .replaceAll(/[12].*/, '')

    script:
    """
    rnabloom -left ${reads[0]} -right ${reads[1]} -outdir ${variant}_${formatted}-Bloom
    """
    //
}

process SPADES {
    publishDir "$params.assembled/rnaspades", mode: "copy"
    conda "/home/shannc/anaconda3/envs/spades"

    input:
    tuple val(variant), path(reads)
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

process PLASS {
    publishDir "$params.translated/plass"
    conda "/home/shannc/anaconda3/envs/plass"
    //
    input:
    tuple val(variant), path(reads)
    //
    output:
    path "${variant}_${formatted}-Plass*"
    //
    exec:
    formatted = reads[0].baseName.replaceAll(/.*-/, '')
                                 .replaceAll(/[12].*/, '')
    //
    script:
    """
    plass assemble ${reads[0]} ${reads[1]} ${variant}_${formatted}-Plass.fas tmp
    """
    // TODO: tmp is probably just a temporary location for plass to hold required files???
}

process TRANSDECODE {
    publishDir "$params.translated/$folder", mode: 'copy'

    input:
    tuple val(folder), path(assembly)
    //
    output:
    path("*")
    //
    exec:
    formatted = "${assembly.baseName}_TR"
    script: // The minimum length has been lowered to 35 from the default of 100 despite the fact that this greatly increases false positive rate. Lowering it to this level is necessary though because venom peptides can be this small
    """
    TransDecoder.LongOrfs -t $assembly \
    -m 35 \
    -O $formatted
    """
    // Although the false positive rate increases with a lower minimum peptide
    // length, venom peptides
}

process QALIGN { // Calculate database coverage by raw reads
    publishDir "$params.alignlog/$assembler"
    conda "/home/shannc/anaconda3/envs/rnaquast"

    input:
    tuple val(filter), val(assembler), path(assembly), path(reads)
    //
    output:
    path(formatted)
    //
    exec:
    formatted = assembly.baseName.replaceAll(/.*fasta/, '')
    //
    script:
    """
    rnaQUAST.py -c $assembly \
    -1 ${reads[0]} \
    -2 ${reads[1]} \
    -o ${formatted}
    """
    //
}

process BALIGN { // Calculate database coverage by raw reads
    publishDir "$outdir"

    input:
    tuple val(filter), val(assembler), path(assembly), path(reads)
    val(outdir)
    //
    output:
    path("align_${formatted}")
    //
    exec:
    formatted = assembly.baseName.replaceAll(/.*fasta/, '')
    //
    script:
    """
    bowtie2-build $assembly $assembly
    bowtie2 \
    -p 10  \
    -q \
    --no-unal \
    -k 20 \
    -x $assembly \
    -1 ${reads[0]} -2 ${reads[1]} \
    2> align_${formatted}.txt
    """
    //
}

