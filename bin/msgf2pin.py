#!/usr/bin/env python
import numpy as np
import pandas as pd
from pathlib import Path
from pyteomics.mzid import MzIdentML
import pyteomics.mzid as ptz


def join_proteins(bad_line: list[str]):
    fits = bad_line[:32]
    the_rest = bad_line[32:]
    if not the_rest:
        return fits
    return fits[:-1] + [','.join([fits[-1]] + the_rest)]


decoy_path = "../tests/results/msgf2pin/decoy_sample.mzid"
valid_path = "../tests/results/msgf2pin/valid_sample.mzid"
real_pin = "../tests/results/msgf2pin/sample_pin.tab"

# subset = valid[( valid["MS-GF:RawScore"] == 78 ) & ( valid["MS-GF:DeNovoScore"] == 79 )]
# subset[subset["calculatedMassToCharge"] <= 600].iloc[0]

compare = pd.read_table(real_pin, sep="\t", on_bad_lines=join_proteins,
                        engine="python")
pin = pd.DataFrame()
decoys = ptz.DataFrame(decoy_path)
valid = ptz.DataFrame(valid_path)
validR = pd.DataFrame(MzIdentML(valid_path))
# validR["SpectrumIdentificationItem"][1]
neutron: int = 1.0033548378
decoys["Label"] = -1
valid["Label"] = 1  # Set to 1 for target psms
valid["ExpMass"] = valid[""]


# pin["SpecId"] =
# pin["Label"] = #
# pin["ScanNr"]
pin["ExpMass"] = round(valid["experimentalMassToCharge"], 3)
pin["CalcMass"] = round(valid["calculatedMassToCharge"], 3)
pin["RawScore"] = valid["MS-GF:RawScore"]
pin["DeNovoScore"] = valid["MS-GF:DeNovoScore"]
pin["ScoreRatio"] = pin["RawScore"]/pin["DeNovoScore"]
pin["Energy"] = pin["DeNovoScore"] - pin["RawScore"]
# pin["lnEvalue"]
pin["IsotopeError"] = valid["IsotopeError"]
pin["lnExplainedIonCurrentRatio"] = (np.log(valid["ExplainedIonCurrentRatio"])
                                     + 0.0001)
pin["lnNTermIonCurrentRatio"] = (np.log(valid["NTermIonCurrentRatio"])
                                 + 0.0001)
pin["lnCTermIonCurrentRatio"] = (np.log(valid["CTermIonCurrentRatio"])
                                 + 0.0001)
pin["lnMS2IonCurrent"] = np.log(valid["MS2IonCurrent"])
pin["Mass"] = valid["experimentalMassToCharge"]
pin["PepLen"] = valid["PeptideSequence"].apply(lambda x: len(x))
pin["dM"] = ((valid["experimentalMassToCharge"] -
             (valid["IsotopeError"] * neutron / valid["chargeState"]) -
             valid["calculatedMassToCharge"]) /
             valid["experimentalMassToCharge"])
pin["absDM"] = abs(pin["dM"])
pin["MeanErrorTop7"] = valid["MeanErrorTop7"].apply(lambda x: -x)
pin["StdevErrorTop7"] = valid["StdevErrorTop7"]
pin["sqMeanErrorTop7"] = pin["MeanErrorTop7"] ** 2
pin["Charge1"] = valid["chargeState"].apply(lambda x: 1 if x == 1 else 0)
pin["Charge2"] = valid["chargeState"].apply(lambda x: 2 if x == 1 else 0)
pin["Charge3"] = valid["chargeState"].apply(lambda x: 3 if x == 1 else 0)
pin["Charge4"] = valid["chargeState"].apply(lambda x: 4 if x == 1 else 0)
pin["Charge5"] = valid["chargeState"].apply(lambda x: 5 if x == 1 else 0)
# pin["enzN"]
# pin["enzC"]
pin["enzInt"] = pin["PeptideSequence"].apply(lambda x:
                                             (lowered := x.lower()).count("k")
                                             + lowered.count("r"))
pin["Peptide"] = valid["pre"].combine(valid["PeptideSequence"],
                                      lambda x, y: f"{x}.{y}")
pin["Peptide"] = valid["Peptide"].combine(valid["post"],
                                          lambda x, y: f"{x}.{y}")
# pin["Proteins"]
