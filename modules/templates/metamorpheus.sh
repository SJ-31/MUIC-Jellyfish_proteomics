mode=!{mode}
metamorpheus -s !{mzmls} \
    -o . \
    -t !{config} \
    -d !{database}

mv Task1*/[Aa]ll* .
for i in [Aa]ll*
    do
    mv $i metamorpheus${mode}_${i}
done

mv metamorpheus${mode}_AllPSMs_FormattedForPercolator.tab edits.tab
cat edits.tab | sed \
    -e "s/DECOY_/rev_/g" \
    -e "s/\[Common Variable:Oxidation on M\]/[15.9949]/g" \
    -e "s/\[Common Fixed:Carbamidomethyl on C\]//g" \
    > metamorpheus${mode}_AllPSMs_FormattedForPercolator.tab
