type = !{type}
case ${type} in
    default)
        config="metamorpheus_params.toml" ;;
    glyco)
        config="metamorpheus_glyco_params.toml" ;;
    *)
      exit 1 ;;
esac

metamorpheus -s !{mzmls} \
    -o . \
    -t !{params.config}/${config} \
    -d !{database}
mv Task1SearchTask/[Aa]ll* .
for i in [Aa]ll*
    do
    mv $i metamorpheus${type}_${i}
    done
mv metamorpheus${type}_AllPSMs_FormattedForPercolator.tab edits.tab
cat edits.tab | sed \
    -e "s/DECOY_/rev_/g" \
    -e "s/\[Common Variable:Oxidation on M\]/[15.9949]/g" \
    -e "s/\[Common Fixed:Carbamidomethyl on C\]//g" \
    > metamorpheus${type}_AllPSMs_FormattedForPercolator.tab
