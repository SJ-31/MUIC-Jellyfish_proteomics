for xml in *.pep.xml
    do
        name=$(echo $xml | sed 's/\..*//')
        identipy2pin $xml -prefix rev_
done

if (( $(ls *.pin | wc -l) > 1)); then
    !{params.bin}/identipy_combine_pins.py identipy_all_pins.temp
else
    mv *pin identipy_all_pins.temp
fi
