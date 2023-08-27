#!/usr/bin/env ipython
from pathlib import Path
import pandas as pd

def make_manifest(template, file_path):
    if not isinstance(template, pd.DataFrame):
        template_df = pd.read_csv(template, sep="\t")
    else:
        template_df = template
    mgf = [str(f.absolute()) for f in Path(file_path).glob("*.mgf")]
    mzML = [str(f.absolute()) for f in Path(file_path).glob("*.mzML")]
    if template_df.columns.to_list() != ['Raw', 'indexed_mzML', 'mgf']:
        raise ValueError('''Unsupported format\nThe template file should have
        headers "Raw", "indexed_mzML" and "mgf"''')
    template_df["mgf"] = mgf
    template_df["indexed_mzML"] = mzML
    return template_df


def parse_args():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--file_path")
    parser.add_argument("-t", "--template")
    parser.add_argument("-o", "--output")
    args = vars(parser.parse_args())  # convert to dict
    return args

if __name__ == '__main__':
    args = parse_args()
    df = make_manifest(args["template"], args["file_path"])
    df.to_csv(args["output"], sep="\t", index=False)
