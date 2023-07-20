process MAXQUANT {
    publishDir "$outdir", mode: "copy"

    input:
    path(raw_file)
    val(mq_config)
    val(outdir)
    //
    output:
    path("${raw_file.baseName}_*.txt")

    shell:
    '''
    cp !{mq_config} template.xml
    maxquant_wrapper.py !{raw_file} template.xml
    dotnet /home/shannc/tools/MaxQuant_v2.3.1.0/bin/MaxQuantCmd.exe mqconfig.xml
    mv combined/txt/*txt .
    mv 'Oxidation (M)Sites.txt' Oxidation_M_sites.txt
    for txt in *txt
        do
            mv $txt !{raw_file.baseName}_${txt}
        done
    '''
    //
}
