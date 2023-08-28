#!/usr/bin/env ipython
#

import pandas as pd
import re
from pyteomics import parser
from pyteomics import mass

# from argparse import ArgumentParser
# parser = ArgumentParser()
# parser.add_argument("-m", "--seq_mapping")
# parser.add_argument("-i", "--intersected_searches")
# parser.add_argument("-r", "--mass_range")  # e.g. 350-1600
# args = parser.parse_args()
# prot_df = pd.read_csv(args.intersected_searches, sep="\t")
# mapping = pd.read_csv(args.seq_mapping, sep="\t")
# mass_range = (int(m) for m in parser.mass_range.split("-"))

# Testing
# prot_df = pd.read_csv("../results/ref/intersected_searches.tsv", sep="\t")
# mapping = pd.read_csv("../results/ref/all_normal_mapping.tsv", sep="\t")


def clean_pep(peptide):
    peptide = peptide.upper()
    return ''.join(re.findall("[^BXZ]+", peptide))

def n_observed(pep_list):
    peps = set(pep_list.split(","))
    if "NA" in peps:
        peps.remove("NA")
    return len(peps)

def n_observable(protein_id, mapping, mass_range):
    try:
        seq = mapping.loc[protein_id]["seq"]
    except KeyError:
        print("Id not found")
        return None
    pep_df = pd.DataFrame({"obsv_pep": list(parser.cleave(seq, 'trypsin'))})
    pep_df["obsv_pep"] = pep_df["obsv_pep"].apply(clean_pep)
    pep_df["mass"] = (pep_df["obsv_pep"]
                      .apply(lambda x: mass.calculate_mass(sequence=x)))
    pep_df["In_range"] = (pep_df["mass"]
                        .apply(lambda x: mass_range[0] <= x <= mass_range[1]))
    if pep_df["In_range"].any():
        return(pep_df["In_range"].value_counts()[True])
    return None

def emPAI(observed, observable):
    return 10**(observed/observable) - 1

def calculate_emPAI(df, mapping, m_range):
    not_full = df["header"].str.contains("-DENOVO|-TRANSCRIPTOME", regex=True)
    full = df[~not_full]
    other = df[not_full]
    mapping.index = mapping["id"]
    full["observable_peps"] = full["ProteinId"].apply(n_observable,
                                                      mapping=mapping,
                                                      mass_range=m_range)
    full["observed_peps"] = full["peptideIds"].apply(n_observed)
    full["emPAI"] = full.apply(lambda x: emPAI(x.observed_peps,
                                            x.observable_peps), axis=1)
    full = full.drop(["observable_peps", "observed_peps"], axis="columns")
    df = pd.concat([full, other])
    return df
