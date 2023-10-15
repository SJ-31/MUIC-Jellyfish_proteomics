#!/usr/bin/env python

from pathlib import Path
import re
import pandas as pd


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
    percolator = (pd.read_csv(filepath, index_col=False).iloc[:, 0].apply(
        str.split, sep="\t").apply(perc_row).to_list())
    final = pd.concat(percolator)
    final["posterior_error_prob"] = final["posterior_error_prob"].astype(float)
    final["q-value"] = final["q-value"].astype(float)
    final = final[final["q-value"] <= q_threshold]
    final = final[final["posterior_error_prob"] <= p_threshold]
    final["peptide"] = final['peptide'].apply(clean_peptide)
    return final


def read_tide(filepath, p_threshold, q_threshold):
    tide = pd.read_csv(filepath, sep="\t")
    tide["sequence"] = tide["sequence"].apply(clean_peptide)
    tide = (tide[tide["percolator q-value"] <= q_threshold])
    tide = (tide[tide["percolator PEP"] <= p_threshold].rename(
        columns={
            "sequence": "peptide",
            "protein id": "proteinIds"
        }))
    return tide


def make_fasta(pep):
    return f">U{pep.name}\n{pep[0]}\n"


def get_engine_files(path) -> dict:
    files = [file.absolute() for file in Path(path).glob("*_percolator_*")]
    names = [re.sub(r"_percolator_.*", "", f.name) for f in files]
    engine_files = dict(zip(names, files))
    return engine_files


def parse_args():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input_path")
    parser.add_argument("-o", "--output")
    parser.add_argument("-t", "--unmatched_tsv")
    parser.add_argument("-p", "--pep_threshold")
    parser.add_argument("-q", "--q_threshold")
    args = vars(parser.parse_args())  # convert to dict
    return args


if __name__ == "__main__":
    # args = {
    #     "pep_threshold": 1,
    #     "q_threshold": 0.05,
    #     "output": "unmatched_peptides.fasta",
    #     "input_path": "../tests/percolator_psms"
    # }
    args = parse_args()
    all_unmatched = set()
    unmatched_df = []
    pep_threshold = float(args["pep_threshold"])
    q_threshold = float(args["q_threshold"])
    path = args["input_path"]
    engine_files = get_engine_files(path)
    for engine, path in engine_files.items():
        if engine == "tide":
            current = read_tide(engine_files["tide"], pep_threshold,
                                q_threshold)
        else:
            current = read_percolator(path, pep_threshold, q_threshold)
        unmatched = pd.Series(
            current.where(current["proteinIds"] == '')["peptide"].unique())
        unmatched_df.append(current[current["peptide"].isin(unmatched)])
        all_unmatched = all_unmatched | set(unmatched)
    all_unmatched = (pd.Series(list(all_unmatched)).reset_index().apply(
        make_fasta, 1))
    all_unmatched_df = pd.concat(unmatched_df).groupby("peptide").sample()
    all_unmatched_df.to_csv(args["unmatched_tsv"])
    with open(args["output"], "w") as f:
        f.write(''.join(all_unmatched))
