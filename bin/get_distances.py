#!/usr/bin/env python
import sys

# Compute cosine and euclidean distance matrices from given embeddings

import numpy as np
import h5py
from pathlib import Path
import pandas as pd

h5py.get_config().track_order = True
SHUFFLE_SEED = 312002
SAMPLE_SEED = 14324


def firstKey(d: dict) -> str:
    return list(d.keys())[0]


def pytorchEmbedding(path, name_list, embd_list) -> None:
    import torch

    loaded = torch.load(path)
    name_list.append(loaded["label"])
    rep = loaded["mean_representations"]
    embd_list.append(rep[firstKey(rep)].numpy())


def getEmbeddings(path) -> tuple:
    names = []
    embeddings = []
    if ".hdf5" in path:
        with h5py.File(path, "r") as h:
            for entry in h.keys():
                names.append(entry)
                embeddings.append(h[entry][:])
    elif Path(path).is_dir():
        for p in Path(path).glob("*.pt"):
            pytorchEmbedding(p, names, embeddings)

    embeddings = np.array(embeddings)
    return names, embeddings


def writeEmbeddingsHDF5(path, names, embeddings) -> None:
    with h5py.File(path, "w") as efile:
        for n, e in zip(names, embeddings):
            efile.create_dataset(n, data=e)


def filterEmbeddings(path, tsv, criteria, id_col="ProteinId") -> tuple:
    df = pd.read_csv(tsv, sep="\t").query(criteria)
    embd_df = hdf5ToDf(path)
    embd_df = embd_df[embd_df.index.isin(df[id_col])]
    return df, embd_df


def hdf5ToDf(path) -> pd.DataFrame:
    names, embeddings = getEmbeddings(path)
    embd_df = pd.DataFrame(embeddings).rename(lambda x: f"V{x}", axis="columns")
    embd_df.index = names
    return embd_df


def getSaved(embd_path, dist_path) -> dict:
    """
    Read in pre-computed embeddings and distances, arranging them into a
    pandas dataframe to send to R
    :param embd_path: path to HDF5 file containing embeddings
    :param dist_path: path to HDF5 file containing distance matrices
    :return: a dictionary with three entries
    """
    output = {}
    embd_df = hdf5ToDf(embd_path)
    output["embeddings"] = embd_df
    with h5py.File(dist_path, "r") as dfile:
        output["euclidean"] = pd.DataFrame(dfile["metric/euclidean"])
        output["euclidean"].index = embd_df.index
        output["cosine"] = pd.DataFrame(dfile["metric/euclidean"])
        output["cosine"].index = embd_df.index
    return output


def parse_args():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input")
    parser.add_argument("-o", "--output")
    parser.add_argument("-w", "--write_embd_file")

    # path to the metadata for the embeddings in "input"
    # meant to be the sample
    parser.add_argument("--sample_tsv")
    # Query string used to filter the sample embeddings
    parser.add_argument("--filter_criteria")
    # hdf5 containing protein embeddings to compare sample
    # proteins against
    parser.add_argument("-c", "--comparison_embd")
    # path to metadata for the comparison embeddings
    parser.add_argument("--comparison_tsv")

    args = vars(parser.parse_args())
    return args


def writeDistances(embeddings, names, file) -> None:
    """
    Compute distance matrix from `embeddings`, and save all output
    to hdf5 file `file`
    """
    embeddings = embeddings.astype(np.float32)
    from sklearn import metrics

    euclidean = metrics.pairwise_distances(embeddings, metric="euclidean")
    cosine = metrics.pairwise_distances(embeddings, metric="cosine")
    with h5py.File(file, "w") as dfile:
        dfile.create_dataset("metric/euclidean", data=euclidean)
        dfile.create_dataset("metric/cosine", data=cosine)
        dfile.create_dataset("names", data=list(names))


def shuffleSample(df, n) -> pd.DataFrame:
    if df.shape[0] <= n:
        return df
    shuffled = df.sample(frac=1, random_state=SHUFFLE_SEED)
    return shuffled.sample(n=n, random_state=SAMPLE_SEED)


if __name__ == "__main__" and len(sys.argv) > 1 and "radian" not in sys.argv[0]:
    # Prevents reticulate from entering this chunk

    names: list
    emeddings: np.array

    args = parse_args()
    if not args["comparison_embd"]:
        names, embeddings = getEmbeddings(args["input"])
    else:
        # For comparing the sample proteins against proteins
        # of different taxonomic lineages, obtained from UniProt
        sample_meta, sample_embd = filterEmbeddings(
            args["input"], args["sample_tsv"], args["filter_criteria"]
        )
        comparison_embd = hdf5ToDf(args["comparison_embd"])
        comparison_meta = pd.read_csv(args["comparison_tsv"], sep="\t")

        # If a given lineage has more proteins than the sample, obtain a
        # random slice that lineage so that the protein numbers are comparable
        min_length = sample_meta["length"].min()
        comparison_meta.query("Length >= @min_length", inplace=True)
        n_embeddings = sample_meta.shape[0]

        sampled = (
            comparison_meta.groupby("Taxon")
            .apply(
                lambda x: shuffleSample(x, n=n_embeddings),
                include_groups=False,
            )
            .reset_index()
        )
        combined = pd.concat(
            [
                comparison_embd.query("index.isin(@sampled['Entry'])"),
                sample_embd,
            ]
        )
        names = combined.index
        embeddings = np.array(combined)
    writeDistances(embeddings, names, args["output"])
    if args["write_embd_file"]:
        writeEmbeddingsHDF5(args["write_embd_file"], names, embeddings)
