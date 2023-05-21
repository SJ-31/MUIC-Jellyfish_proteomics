process MAXQUANT {
    input:
    path(files)
    //
    script:
    """
    dotnet /home/shannc/tools/MaxQuant_v2.3.1.0/bin/MaxQuantCmd.exe $files
    """
    //
}
