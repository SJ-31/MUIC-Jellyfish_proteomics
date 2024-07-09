#!/usr/bin/env python

import sys
import numpy as np
import os
from subprocess import run, CalledProcessError
import tempfile
import polars as pl
import igraph as ig
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram
import leidenalg as la


def save_dendogram(
    linkage_matrix: np.array,
    filename: str,
    cutoff: float,
    size: tuple = (20, 10),
):
    fig, ax = plt.subplots()
    dendrogram(linkage_matrix, no_labels=True, ax=ax)
    ax.set_ylabel(f"height (cut at {cutoff})")
    ax.axhline(0.1, linestyle="dotted", color="black")
    fig.set_size_inches(size[0], size[1])
    fig.savefig(filename, bbox_inches="tight")
    return fig, ax


def linkage_matrix(fitted_model):
    counts = np.zeros(fitted_model.children_.shape[0])
    n_samples = len(fitted_model.labels_)
    for i, merge in enumerate(fitted_model.children_):
        current_count = 0
        for child_idx in merge:
            if child_idx < n_samples:
                current_count += 1  # leaf node
            else:
                current_count += counts[child_idx - n_samples]
        counts[i] = current_count

    linkage_matrix = np.column_stack(
        [fitted_model.children_, fitted_model.distances_, counts]
    ).astype(float)
    return linkage_matrix


def create_graph(nd_array, names) -> ig.Graph:
    """
    Create weighted graph from numpy array with names
    :param nd_array:
    :param names:
    :return:
    """
    graph = ig.Graph.Weighted_Adjacency(nd_array, mode="undirected")
    graph.vs["name"] = names
    return graph


# # Modularity vertex partition doesn't work


def partition_metrics(partition) -> dict:
    metrics = {
        "n_clusters": len(partition),
        "n_elements": partition.n,
        "quality": partition.quality(),
        "modularity": partition.modularity,
    }
    return metrics


def optimise_partition(graph: ig.Graph, p_type, it=10) -> la.VertexPartition:
    """
    Create a partition and attempt to optimize it.
    It is not possible to decide if a partition is optimal
    :param graph: Graph
    :param p_type: One of the partition types provided by the leidenalg
    package
    :param it: Number of iterations to optimize for
    :return:
    """
    optimiser = la.Optimiser()
    p = la.find_partition(graph, partition_type=p_type)
    optimiser.optimise_partition(partition=p, n_iterations=it)
    return p


def write_fasta(headers: list[str], seqs: list[str], filename: str) -> None:
    text = "\n".join([f">{h}\n{s}" for h, s in zip(headers, seqs)])
    with open(filename, "w") as w:
        w.write(text)


def mmseqs_cluster(
    headers: list[str],
    seqs: list[str],
    mmseqs_bin: str,
    name: str,
    flags: tuple = ("--cov-mode", "1", "--cluster-mode", "2", "--cluster-reassign"),
) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Cluster proteins

    Args:
        headers (list[str]): sequence headers
        seqs (list[str]): sequences to cluster
        mmseqs_bin (str): path to mmseqs binary
        name (str): name to organize sequences by
        flags (tuple, optional): Arguments for mmseqs easy-cluster.
            Defaults to ("--cov-mode", "1", "--cluster-mode", "2", "--cluster-reassign"), which are recommended for dealing with
        protein fragments in mmseqs documentation

    Returns:
        tuple[pl.DataFrame, pl.DataFrame]: df of representative sequences,
        df mapping representatives to their members
    """
    if len(headers) != len(seqs):
        raise ValueError("Number of headers and sequences do not match")
    reps: dict = {"ProteinId": [], "Group": []}
    with tempfile.TemporaryDirectory() as tmpdir:
        pop: str = os.getcwd()
        os.chdir(tmpdir)
        write_fasta(headers, seqs, "seqs.fasta")
        try:
            run(
                " ".join(
                    (mmseqs_bin, "easy-cluster", "seqs.fasta", "result", "tmp") + flags
                ),
                shell=True,
            ).check_returncode()
            mapping: pl.DataFrame = pl.read_csv(
                "result_cluster.tsv",
                separator="\t",
                new_columns=["representative", "member"],
            )
            reps: pl.DataFrame = mapping.unique("representative")
            reps = reps.with_columns(
                Group=pl.Series([f"{name}.{i}" for i in range(reps.shape[0])])
            ).drop("member")
        except CalledProcessError:
            print("Searching failed")
            print(f"\tGiven headers \n{headers}")
            print(f"\tGiven seqs \n{seqs}")
        finally:
            os.chdir(pop)
    return reps, mapping


def get_group_representatives(
    data: pl.DataFrame, mmseqs_bin: str
) -> tuple[pl.DataFrame, pl.DataFrame]:
    grouped: pl.DataFrame = (
        data.group_by(pl.col("Group"))
        .agg(pl.col("ProteinId"), pl.col("seq"), size=pl.len())
        .sort(pl.col("size"), descending=True)
    )
    reps: list = []
    groups: list = []
    mappings: list = []
    for group, ids, seqs, size in grouped.iter_rows():
        if size > 1 and len(set(seqs)) > 1:
            df, mapping = mmseqs_cluster(ids, seqs, mmseqs_bin, group)
            try:
                reps.extend(df["representative"])
                groups.extend(df["Group"])
                mappings.append(mapping)
            except pl.exceptions.ColumnNotFoundError as e:
                print(e)
                print("Column not found error")
                print(df)
                exit(1)
        else:
            reps.extend(ids)
            groups.extend([group] * len(ids))
            mappings.append(
                pl.DataFrame(
                    {"representative": pl.Series(ids), "member": pl.Series(ids)}
                )
            )
    map_df = pl.concat(mappings)
    rep_df = pl.DataFrame({"ProteinId": reps, "Subgroup": groups}).join(
        data, on="ProteinId"
    )
    data = data.join(map_df, left_on="ProteinId", right_on="member", how="left").rename(
        {"representative": "Group_representative"}
    )
    return data, rep_df


def parse_args():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-m", "--mmseqs")
    parser.add_argument("-i", "--input")
    parser.add_argument("-o", "--output")
    args = vars(parser.parse_args())
    return args


if __name__ == "__main__" and len(sys.argv) > 1:
    args = parse_args()
    df = pl.read_csv(args["input"], separator="\t", null_values="NA")
    if "Group_representative" in df.columns:
        df = df.drop("Group_representative")
    df, representatives = get_group_representatives(df, args["mmseqs"])
    df.write_csv(args["input"], separator="\t")
    representatives.write_csv(args["output"], separator="\t")
