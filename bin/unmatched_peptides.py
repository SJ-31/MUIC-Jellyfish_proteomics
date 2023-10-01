#!/usr/bin/env ipython

from pathlib import Path
import re
from argparse import ArgumentParser
import pandas as pd

parser = ArgumentParser(
    prog="ls",  # Name of the CLI program
    description="",
    epilog="",  # Final statement printed to terminal
)
parser.add_argument("output")
parser.add_argument("pep_threshold")
parser.add_argument("q_threshold")
args = parser.parse_args()


def perc_row(series_row):
    return pd.DataFrame(
        {
            "PSMId": series_row[0],
            "score": series_row[1],
            "q-value": series_row[2],
            "posterior_error_prob": series_row[3],
            "peptide": series_row[4],
            "proteinIds": ';'.join(series_row[5:])
        },
        index=[1])


def clean_peptide(peptide):
    peptide = peptide.upper().replace("X", "")
    return ''.join(re.findall("[A-Z]+", peptide))


def read_percolator(filepath, p_threshold, q_threshold):
    percolator = (pd.read_csv(filepath, index_col=False).iloc[:, 0]
                  .apply(str.split, sep="\t")
                  .apply(perc_row)
                  .to_list())
    final = pd.concat(percolator)
    final["posterior_error_prob"] = final["posterior_error_prob"].astype(float)
    final = final[final["q-value"] <= q_threshold]
    final = final[final["posterior_error_prob"] <= p_threshold]
    final["peptide"] = final['peptide'].apply(clean_peptide)
    return final


def read_tide(filepath, p_threshold, q_threshold):
    tide = pd.read_csv(filepath, sep="\t")
    tide["sequence"] = tide["sequence"].apply(clean_peptide)
    tide = (tide[tide["percolator q-value"] <= q_threshold])
    tide = (tide[tide["percolator PEP"] <= p_threshold]
            .rename(columns={"sequence": "peptide",
                             "protein id": "proteinIds"}))
    return tide


def make_fasta(pep):
    return f">U{pep.name}\n{pep[0]}\n"


files = [file.absolute() for file in Path(".").glob("*_percolator_*")]
names = [re.sub(r"_percolator_.*", "", f.name) for f in files]
engine_files = dict(zip(names, files))

all_unmatched = set()
pep_threshold = float(args.pep_threshold)
q_threshold = float(args.q_threshold)

for engine, path in engine_files.items():
    if engine == "tide":
        current = read_tide(engine_files["tide"], pep_threshold, q_threshold)
    else:
        current = read_percolator(path, pep_threshold, q_threshold)
    unmatched = set(current.where(current["proteinIds"] == '')["peptide"])
    all_unmatched = all_unmatched | unmatched


all_unmatched = (pd.Series(list(all_unmatched))
                 .reset_index()
                 .apply(make_fasta, 1))

with open(args.output, "w") as f:
    f.write(''.join(all_unmatched))
