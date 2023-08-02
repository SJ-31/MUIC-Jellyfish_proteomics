process MAXQUANT {
    stageInMode "copy"
    publishDir "$outdir", mode: "copy"

    input:
    path(raw_file)
    val(outdir)
    //
    output:
    path("${raw_file.baseName}_combined")
    path("*txt")

    shell:
    '''
    cp !{params.mqPars} template.xml
    maxquant_wrapper.py !{raw_file} template.xml
    rm template.xml
    dotnet /home/shannc/tools/MaxQuant_2.4.2.0/bin/MaxQuantCmd.exe mqconfig.xml
    for file in combined/txt/*
        do
            new=$(echo $file | cut -d "/" -f 2)
            cp $file !{raw_file.baseName}_${new}
    done
    mv combined/txt ./!{raw_file.baseName}_combined
    '''
    //
}
