#!/usr/bin/env python

import pandas as pd
import re
import numpy as np


def unify_mods(peptide):
    if peptide == " ":
        return np.NaN
    peptide = peptide.replace("_", "")
    peptide = peptide.replace("(Oxidation (M))", "[15.9949]")
    if get := re.search("\\(Acetyl \\(Protein N-term\\)\\)", peptide):
        replacement = "n[42.0106]"
        peptide = peptide[get.span()[1]] + replacement + peptide[get.span()[1] + 1 :]
    peptide = f"-.{peptide}.-"
    return peptide


def split_prot(prot):
    individual = prot.split(";")
    prot = ";".join([p.replace("REV__", "rev_") for p in individual])
    return prot


def format_mq(file):
    sample = pd.read_csv(file, sep="\t")
    pin = pd.DataFrame()
    pin["SpecId"] = sample.agg(lambda x: f'{x["Raw file"]}.{x["Scan number"]}', axis=1)
    pin["Label"] = sample["Reverse"].apply(lambda x: -1 if x == "+" else 1)
    pin["ScanNr"] = sample["Scan number"]
    selection = sample.loc[
        :,
        [
            "Retention time",
            "Total ion current",
            "Base peak intensity",
            "Precursor intensity",
            "Precursor apex fraction",
            "Precursor apex offset",
            "m/z",
            "MS2 m/z",
            "Charge",
            "Length",
            "Filtered peaks",
            "Intens Comp Factor",
            "Precursor apex offset time",
            "CTCD Comp",
            "RawOvFtT",
            "AGC Fill",
        ],
    ].reset_index(drop=True)
    pin = pd.concat([pin, selection], axis=1)
    pin["Peptide"] = sample["Modified sequence"].apply(unify_mods)
    pin["Proteins"] = sample["Proteins"].apply(split_prot)
    return pin.dropna(axis="index").reset_index(drop=True)


def parse_args():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input")
    parser.add_argument("-o", "--output")
    args = vars(parser.parse_args())  # convert to dict
    return args


if __name__ == "__main__":
    args = parse_args()
    pin = format_mq(args["input"])
    pin.to_csv(args["output"], sep="\t", index=None)
