import networkx as nx
import json
import sys
import argparse
import polars as pl
from collections import ChainMap
import obonet


class SubsetGO:
    """Class for subsetting the GO DAG based on the GO terms found in
    a protein sample.
    """

    def __init__(self, go_path: str, sample_path: str, metadata_path: str) -> None:
        metadata = pl.read_csv(metadata_path, separator="\t")
        self.roots = {"BP": "GO:0008150", "CC": "GO:0005575", "MF": "GO:0003674"}
        sample: pl.DataFrame = pl.read_csv(
            sample_path, separator="\t", null_values="NA"
        )
        self.sample_gos: set = {
            x
            for go_string in sample["GO_IDs"].drop_nulls()
            for x in go_string.split(";")
        }
        GO: nx.MultiDiGraph = obonet.read_obo(go_path)
        self.G: nx.MultiDiGraph = nx.MultiDiGraph()
        root_map: dict = dict(zip(metadata["GO_IDs"], metadata["ontology"]))
        GO = obonet.read_obo(go_path)
        self.G.add_nodes_from(self.roots.values())
        for go in self.sample_gos:
            if go in GO:
                paths: list = nx.all_simple_edge_paths(
                    GO, source=go, target=self.roots[root_map[go]]
                )
                is_a = relationPath(paths, "is_a")
                if is_a:
                    self.G.add_edges_from(is_a)
        self.metadata = self.__getNodeData().join(
            metadata.select("GO_IDs", "term", "definition", "ontology"), on="GO_IDs"
        )

    def __getNodeData(self):
        self.successors: dict = {}
        level_map: dict = {}
        from_sample: dict = {}
        # GO_IDs found in the sample are True, others (used to link GO terms back to their roots) are False
        for node in self.G.nodes():
            if node in self.sample_gos:
                from_sample[node] = True
            else:
                from_sample[node] = False
        nx.set_node_attributes(self.G, {"from_sample": from_sample})
        for root in self.roots.values():
            current = nx.bfs_tree(self.G, root)
            current.graph["root"] = root
            self.successors = ChainMap(self.successors, allSuccessors(current))
            level_map = ChainMap(level_map, levelMap(current))

        return (
            (
                pl.DataFrame(from_sample)
                .melt()
                .rename({"variable": "GO_IDs", "value": "in_sample"})
            )
            .with_columns(
                n_children=pl.col("GO_IDs").map_elements(
                    lambda x: len(self.successors[x]), return_dtype=pl.Int16
                ),
                level=pl.col("GO_IDs").map_elements(
                    lambda x: level_map[x], return_dtype=pl.Int16
                ),
            )
            .sort("n_children", descending=True)
        )

    def __checkRepresentative(self, to_check: str, representatives: list):
        """Make sure the roots do not contain each other"""
        return all([to_check not in self.successors[e] for e in representatives])

    def getRepresentatives(
        self, ontology: str, n: int = 18, show=True, min_depth=2, pre=tuple()
    ):
        """Select representative GO terms from the specified ontology that partition the ontology into `n` bins.
        Goal is to map all child terms below the representatives to each representative,
        making for more concise summarization
        representative terms are selected by the number of children they have, as
        well as the specified depth. Or can be pre-specified with the `pre` argument

        :return: A dictionary of the following
        `map`: Map of GO terms in the specified sub-ontology to the representatives above them
        `representatives`: Map of chosen representatives to their terms
        `unassigned`: GO terms in the sub-ontology graph that are not children of
        any representative. Happens when representatives have few children
        """
        children_per_cat = self.metadata.shape[0] / n
        filtered = self.metadata.filter(
            (
                (pl.col("ontology") == ontology)
                & (pl.col("n_children") <= children_per_cat)
                & (pl.col("level") >= min_depth)
            )
            | (pl.col("GO_IDs").is_in(pre))
        ).sort("n_children", descending=True)
        representatives = pre if pre else []
        for go_id, *_ in filtered.iter_rows():
            if self.__checkRepresentative(go_id, representatives):
                representatives.append(go_id)
            if len(representatives) == n:
                break
        representative_df = filtered.filter(pl.col("GO_IDs").is_in(representatives))
        if show:
            print(representative_df)
        child_to_parent: dict = {}
        unassigned: set = set(filtered["GO_IDs"])
        for representative in representatives:
            for child in self.successors[representative]:
                child_to_parent[child] = representative
            child_to_parent[representative] = representative
        unassigned = unassigned - child_to_parent.keys()
        return {
            "map": child_to_parent,
            "representatives": dict(
                zip(representative_df["GO_IDs"], representative_df["term"])
            ),
            "unassigned": unassigned,
        }


def relationPath(paths: list[tuple], relation: str) -> list:
    """Find the path from a list of paths (which are edge lists)
    that consists only of `relation`
    """
    result = []
    for path in paths:
        if all(map(lambda x: x[2] == relation, path)):
            for p in path:
                result.append((*p[:2][::-1], p[2]))
        # Note: This changes the direction of the obonet GO graph so that successors are children and predecessors are ancestors
        break
    return result


def joinDict(df: pl.DataFrame, dct: dict, by: str, colname: str):
    if isinstance(dct, ChainMap):
        dct = dict(dct)
    temp_df = pl.DataFrame(dct).melt().rename({"variable": by, "value": colname})
    return df.join(temp_df, on=by)


def allSuccessors(G: nx.DiGraph) -> dict:
    """Return a dictionary mapping the nodes of G to ALL of their
    children/successors (unlike an adjacency list, which
    only returns the immediate) children
    """
    successors = {}
    for node in G.nodes:
        tree = nx.bfs_tree(G, node)  # Creating a bfs tree for each
        # node restricts the view of G to `node` and everything below it
        children = list(tree.nodes)
        children.remove(node)
        successors[node] = children
    return successors


def levelMap(G: nx.DiGraph, root=None) -> dict:
    """Return a dictionary mapping the nodes of G to their levels in G"""
    level_map: dict = {}
    if not (root := G.graph.get("root", root)):
        raise ValueError("Root must be specified!")
    for i, layer in enumerate(nx.bfs_layers(G, root)):
        for node in layer:
            level_map[node] = i
    return level_map


# Group and find the lowest common ancestors???
# current = GoSubDag(["GO:0031577"], path)
# current.rcntobj.go2ancestors["GO:0031577"]
def parseArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument("-g", "--go_path")
    parser.add_argument("-s", "--sample_path")
    parser.add_argument("-p", "--go_info_path")
    parser.add_argument("-n", "--n_groups")
    parser.add_argument("--representatives")
    parser.add_argument("-o", "--outdir")
    args = vars(parser.parse_args())
    return args


def main(args):
    representatives = {"CC": None, "BP": None, "MF": None}
    S = SubsetGO(args["go_path"], args["sample_path"], args["go_info_path"])
    for sub in representatives.keys:
        rep = S.getRepresentatives()
        representatives[sub] = {"map": rep["map"], "unassigned": rep["unassigned"]}
    with open(f'{args["outdir"]}/go_representatives.json', "w") as j:
        json.dump(representatives, j)


if __name__ == "__main__" and not "ipykernel" in sys.argv[0]:
    args = parseArgs()
    main(args)
