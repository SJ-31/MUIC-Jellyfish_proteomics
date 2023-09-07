for xml in *.pep.xml
    do
        name=$(echo $xml | sed 's/\..*//')
        identipy2pin $xml -prefix rev_
done

!{params.bin}/identipy_combine_pins.py identipy_all_pins.temp
