import numpy as np
import igraph as ig
import leidenalg as la


# Or perhaps run it on a knn version of the embeddings?


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


# m = pd.read_csv("/home/shannc/my_mat.txt", sep="\t")
# m = m.iloc[0:2000, 0:2000]
# nd = np.array(m)
#
# graph = createGraph(nd, m.index)
# md = la.find_partition(graph, partition_type=la.ModularityVertexPartition)
#
# rb = partitionOptimise(graph, la.RBERVertexPartition)
# nl = membershipWNames(rb)
# rb_sets = clusterSets(rb)
