#!/usr/bin/env python
import numpy as np
import subprocess
import re
import pandas as pd
import pyteomics.mzid as ptz
NEUTRON: int = 1.0033548378


def trypsin_spec(peptide, side):
    sites = {'K', 'R'}
    limiter = 'P'
    if side == "N":
        residue = peptide[0]
        next_res = peptide[2]
        return int((residue in sites and next_res != limiter) or
                   (residue == '-'))
    # side == "C"
    residue = peptide[-3]
    next_res = peptide[-1]
    return int((residue in sites and next_res != limiter) or
               (next_res == '-'))


def cleave(peptide: str, enzyme: str, side: str):
    if enzyme == "Trypsin":
        return trypsin_spec(peptide, side)
        # Does trypsin cut before proline?!


def get_mods(row):
    pep = row["PeptideSequence"]
    modifications = row["Modification"]
    if isinstance(modifications, list) and len(modifications) > 0:
        for mod in modifications:
            mass = mod["monoisotopicMassDelta"]
            loc = mod["location"] - 1
            pep = pep[:loc+1] + f"[{mass}]" + pep[loc+1:]
    return pep


def get_spec_id(msgf_row):
    filename = re.findall("(.*)\\.",  msgf_row["name"])[0]
    specid = msgf_row.get("spectrumID")
    if specid:
        specid = "SII" + specid.replace("index=", "")
    else:
        specid = "R" + msgf_row.get("rank")
    return f"{filename}_{specid}"


def format_mzid(mzid: pd.DataFrame, is_decoy):
    pin = pd.DataFrame()
    mzid = mzid[~mzid["StdevErrorAll"].isna()]  # This should remove
    # hits without features
    pin["SpecId"] = mzid.apply(get_spec_id, axis="columns")
    if is_decoy:
        pin["Label"] = -1
    else:
        pin["Label"] = 1
    pin["ScanNr"] = mzid.index.to_series() + 1
    pin["ExpMass"] = round(mzid["experimentalMassToCharge"], 3)
    pin["CalcMass"] = round(mzid["calculatedMassToCharge"], 3)
    pin["RawScore"] = mzid["MS-GF:RawScore"]
    pin["DeNovoScore"] = mzid["MS-GF:DeNovoScore"]
    pin["ScoreRatio"] = pin["RawScore"]/pin["DeNovoScore"]
    pin["ScoreRatio"] = pin["ScoreRatio"].apply(lambda x: 0 if np.isinf(x)
                                                else x)
    pin["Energy"] = pin["DeNovoScore"] - pin["RawScore"]
    pin["lnEValue"] = -np.log(mzid["MS-GF:EValue"])
    pin["IsotopeError"] = mzid["IsotopeError"]
    pin["lnExplainedIonCurrentRatio"] = np.log(mzid["ExplainedIonCurrentRatio"]
                                               + 0.0001)
    pin["lnNTermIonCurrentRatio"] = np.log(mzid["NTermIonCurrentRatio"]
                                           + 0.0001)
    pin["lnCTermIonCurrentRatio"] = np.log(mzid["CTermIonCurrentRatio"]
                                           + 0.0001)
    pin["lnMS2IonCurrent"] = np.log(mzid["MS2IonCurrent"])
    pin["PepLen"] = mzid["PeptideSequence"].apply(lambda x: len(x))
    pin["dM"] = ((mzid["experimentalMassToCharge"] -
                  (mzid["IsotopeError"] * NEUTRON / mzid["chargeState"]) -
                  mzid["calculatedMassToCharge"]) /
                 mzid["experimentalMassToCharge"])
    pin["absdM"] = abs(pin["dM"])
    pin["MeanErrorTop7"] = -mzid["MeanErrorTop7"]
    pin["StdevErrorTop7"] = mzid["StdevErrorTop7"]
    pin["sqMeanErrorTop7"] = pin["MeanErrorTop7"] ** 2
    pin["Charge1"] = mzid["chargeState"].apply(lambda x: 1 if x == 1 else 0)
    pin["Charge2"] = mzid["chargeState"].apply(lambda x: 2 if x == 1 else 0)
    pin["Charge3"] = mzid["chargeState"].apply(lambda x: 3 if x == 1 else 0)
    pin["Charge4"] = mzid["chargeState"].apply(lambda x: 4 if x == 1 else 0)
    pin["Charge5"] = mzid["chargeState"].apply(lambda x: 5 if x == 1 else 0)
    # If the n-terminal peptide is known
    pin["Peptide"] = mzid.apply(get_mods, axis="columns")
    pin["Peptide"] = mzid["pre"].combine(pin["Peptide"],
                                         lambda x, y: f"{x}.{y}")
    pin["Peptide"] = pin["Peptide"].combine(mzid["post"],
                                            lambda x, y: f"{x}.{y}")
    pin["Mass"] = round(mzid["experimentalMassToCharge"], 3)
    pin["enzC"] = pin["Peptide"].apply(cleave, enzyme="Trypsin", side="C")
    pin["enzN"] = pin["Peptide"].apply(cleave, enzyme="Trypsin", side="N")
    pin["enzInt"] = pin["Peptide"].apply(lambda x: (lwr := x.lower())
                                         .count("k") + lwr.count("r"))
    pin["Proteins"] = mzid["accession"].apply(lambda x: '\t'.join(x))
    pin = pin[['SpecId', 'Label', 'ScanNr', 'ExpMass', 'CalcMass', 'RawScore',
               'DeNovoScore', 'ScoreRatio', 'Energy', 'lnEValue',
               'IsotopeError', 'lnExplainedIonCurrentRatio',
               'lnNTermIonCurrentRatio', 'lnCTermIonCurrentRatio',
               'lnMS2IonCurrent', 'Mass', 'PepLen', 'dM', 'absdM',
               'MeanErrorTop7', 'sqMeanErrorTop7', 'StdevErrorTop7', 'Charge1',
               'Charge2', 'Charge3', 'Charge4', 'Charge5', 'enzN', 'enzC',
               'enzInt', 'Peptide', 'Proteins']]
    return pin


def main(args):
    decoys = ptz.DataFrame(args["decoy_path"])
    valid = ptz.DataFrame(args["valid_path"])
    decoy_pin = format_mzid(decoys, True)
    valid_pin = format_mzid(valid, False)
    header = pd.DataFrame({"SpecId": "DefaultDirection", "Label": "-",
                           "ScanNr": "-", "ExpMass": "-", "CalcMass": "-",
                           "RawScore": 0, "DeNovoScore": -1,
                           "ScoreRatio": 0.0, "Energy": -2,
                           "lnEValue": 6.0, "IsotopeError": -2.5,
                           "lnExplainedIonCurrentRatio": 0.0,
                           "lnNTermIonCurrentRatio": 0.0,
                           "lnCTermIonCurrentRatio": 0.0,
                           "lnMS2IonCurrent": 0.0,
                           "Mass": 0.0, "PepLen": 0.0,
                           "dM": 0.0, "absdM": -1.0,
                           "MeanErrorTop7": 0.0, "sqMeanErrorTop7": 0.0,
                           "StdevErrorTop7": 0.0, "Charge1": 0,
                           "Charge2": 0, "Charge3": 0,
                           "Charge4": 0, "Charge5": 0, "enzN": 0,
                           "enzC": 0, "enzInt": 0, "Peptide": None,
                           "Proteins": None, },
                          index=[0])
    final = pd.concat([header, decoy_pin, valid_pin])
    final.to_csv("temp.tsv", sep="\t", index=False)
    subprocess.run(f'cat temp.tsv | sed \'s;";;g\' > {args["output"]}',
                   shell=True)


def parse_args():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-d", "--decoy_mzid")
    parser.add_argument("-v", "--valid_mzid")
    parser.add_argument("-o", "--output")
    args = vars(parser.parse_args())  # convert to dict
    return args


if __name__ == '__main__':
    args = parse_args()
    main(args)
