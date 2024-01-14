import pandas as pd
import re
from pyteomics import parser
from pyteomics import mass


def cleanPep(peptide):
    peptide = peptide.upper()
    peptide = "".join(re.findall("[^BXZ]+", peptide))  # Leaving these
    # unknown residues will raise an error with emPAI
    peptide = "".join(re.findall("[A-Z]+", peptide))
    return peptide


def n_observed(pep_list):
    peps = set(pep_list.split(","))
    if "NA" in peps:
        peps.remove("NA")
    return len(peps)


def n_observable(seq, mass_range):
    seq = cleanPep(seq)
    pep_df = pd.DataFrame({"obsv_pep": list(parser.cleave(seq, "trypsin"))})
    pep_df["mass"] = pep_df["obsv_pep"].apply(
        lambda x: mass.calculate_mass(sequence=x)
    )
    pep_df["In_range"] = pep_df["mass"].apply(
        lambda x: mass_range[0] <= x <= mass_range[1]
    )
    if pep_df["In_range"].any():
        return pep_df["In_range"].value_counts()[True]
    return None


def emPAI(observed, observable):
    return 10 ** (observed / observable) - 1


def calculate_emPAI(df, m_range: tuple):
    not_full = (
        df["header"]
        .str.contains("-DENOVO|-TRANSCRIPTOME", regex=True)
        .fillna(False)
        .combine(df["seq"] == "-", lambda x, y: x or y)
    )
    full = df[~not_full & df["seq"].notna()]
    full["peptideIds"] = full["peptideIds"].apply(
        lambda x: ",".join([cleanPep(p) for p in x.split(",")])
    )
    other = df[not_full | df["seq"].isna()]
    full["observable_peps"] = full["seq"].apply(
        n_observable, mass_range=m_range
    )
    full["observed_peps"] = full["peptideIds"].apply(n_observed)
    full["emPAI"] = full.apply(
        lambda x: emPAI(x.observed_peps, x.observable_peps), axis=1
    )
    full = full.drop(["observable_peps", "observed_peps"], axis="columns")
    df = pd.concat([full, other])
    return df
