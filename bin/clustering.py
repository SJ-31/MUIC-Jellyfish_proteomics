import numpy as np
import igraph as ig
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import dendrogram
import leidenalg as la

# Or perhaps run it on a knn version of the embeddings?


def saveDendogram(
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


def linkageMatrix(fitted_model):
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


def createGraph(nd_array, names) -> ig.Graph:
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


def partitionMetrics(partition) -> dict:
    metrics = {
        "n_clusters": len(partition),
        "n_elements": partition.n,
        "quality": partition.quality(),
        "modularity": partition.modularity,
    }
    return metrics


def partitionOptimise(graph: ig.Graph, p_type, it=10) -> la.VertexPartition:
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
