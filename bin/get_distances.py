#!/usr/bin/env python

# Compute cosine and euclidean distance matrices from given embeddings

import numpy as np
import h5py
import pandas as pd

h5py.get_config().track_order = True


def getEmbeddings(path) -> tuple[list, np.array]:
    protein_ids = []
    embeddings = []
    if ".hdf5" in path:
        with h5py.File(path, "r") as h:
            for entry in h.keys():
                protein_ids.append(entry)
                embeddings.append(h[entry][:])
        embeddings = np.array(embeddings)
    return protein_ids, embeddings


def getSaved(embd_path, dist_path) -> dict[str, np.array]:
    """
    Read in pre-computed embeddings and distances, arranging them into a
    pandas dataframe to send to R
    :param embd_path: path to HDF5 file containing embeddings
    :param dist_path: path to HDF5 file containing distance matrices
    :return: a dictionary with three entries
    """
    output = {}
    names, embeddings = getEmbeddings(embd_path)
    embd_df = pd.DataFrame(embeddings).rename(
        lambda x: f"V{x}", axis="columns"
    )
    output["embeddings"] = embd_df
    embd_df.index = names
    with h5py.File(dist_path, "r") as dfile:
        output["euclidean"] = pd.DataFrame(dfile["metric/euclidean"])
        output["euclidean"].index = names
        output["cosine"] = pd.DataFrame(dfile["metric/euclidean"])
        output["cosine"].index = names
    return output


def parse_args():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input")
    parser.add_argument("-o", "--output")
    args = vars(parser.parse_args())  # convert to dict
    return args


if __name__ == "__main__":
    from sklearn import metrics

    args = parse_args()
    names, embeddings = getEmbeddings(args["input"])
    euclidean = metrics.pairwise_distances(embeddings, metric="euclidean")
    cosine = metrics.pairwise_distances(embeddings, metric="cosine")
    with h5py.File(args["output"], "a") as dfile:
        dfile.create_dataset("metric/euclidean", data=euclidean)
        dfile.create_dataset("metric/cosine", data=cosine)
        dfile.create_dataset("names", data=names)
