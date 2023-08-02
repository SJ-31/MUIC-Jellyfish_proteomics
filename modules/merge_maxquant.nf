process  {
    publishDir "$outdir", mode: "copy"

    input:
    path(maxquant_out)
    val(outdir)
    //

    output:
    path("maxquant_combined*tsv")
    //

    script:
    msms_header="﻿Raw file	Scan number	Retention time	Ion injection time	Total ion current	Collision energy	Summations	Base peak intensity	Elapsed time	Identified	Matched	Reverse	MS/MS IDs	Sequence	Length	Filtered peaks	m/z	Mass	Charge	MS2 m/z	Type	Fragmentation	Mass analyzer	Parent intensity fraction	Fraction of total spectrum	Base peak fraction	Precursor full scan number	Precursor intensity	Precursor apex fraction	Precursor apex offset	Precursor apex offset time	Scan event number	Modifications	Modified sequence	Proteins	Score	PEP	Experiment	MS3 scan numbers	Reporter PIF	Reporter fraction	Intens Comp Factor	CTCD Comp	RawOvFtT	AGC Fill	Scan index	MS scan index	MS scan number"

    """
    merge_tables.sh -r "$msmsheader" \
        -o maxquant_combined_msms.tsv \
        -p msmsScans.txt

    """
    //
}
