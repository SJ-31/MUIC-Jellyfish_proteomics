#!/usr/bin/env python
import sys
import copy
import networkx as nx
import json
import tomllib
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
        self.G.add_nodes_from(self.roots.values())
        for go in self.sample_gos:
            if go in GO:
                paths: list = nx.all_simple_edge_paths(
                    GO, source=go, target=self.roots[root_map[go]]
                )
                is_a = relation_only_paths(paths, "is_a")
                if is_a:
                    for i in is_a:
                        self.G.add_edges_from(i)
        self.metadata = self.get_node_data().join(
            metadata.select("GO_IDs", "term", "definition", "ontology"), on="GO_IDs"
        )

    def save(self, path: str) -> None:
        nx.write_gml(self.G, path)

    def get_node_data(self):
        self.successors: dict = {}  # Map of GO_IDs to list of child terms
        level_map: dict = {}
        from_sample: dict = {}
        # GO_IDs found in the sample are True, others (used to link GO terms back to their roots) are False
        # Produces a data frame that maps GO terms to the number of children they have
        for node in self.G.nodes():
            if node in self.sample_gos:
                from_sample[node] = True
            else:
                from_sample[node] = False
        nx.set_node_attributes(self.G, {"from_sample": from_sample})
        for root in self.roots.values():
            current = nx.bfs_tree(self.G, root)
            current.graph["root"] = root
            self.successors = ChainMap(self.successors, get_successors(current))
            level_map = ChainMap(level_map, get_level_map(current))

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

    def __check_parent(self, to_check: str, ancestors: list):
        """Make sure the parent terms do not contain each other"""
        return all([to_check not in self.successors[e] for e in ancestors])

    def get_predefined(self, path: str, ontology: str) -> dict:
        with open(path, "r") as r:
            custom = json.load(r)
        id2group = {}
        reference = self.metadata.filter(pl.col("ontology") == ontology)
        for group, members in custom.items():
            for member in members:
                if member in reference["GO_IDs"]:
                    id2group[member] = group
        return id2group

    def get_parents(self, ontology: str, n: int = 18, show=True, min_depth=2, pre=""):
        """Select higher-level parent GO terms from the specified ontology that partition the ontology into `n` bins.
        Goal is to map all child terms to some higher-level term to make for more concise summarization
        higher-level terms are selected by the number of children they have, as well as the specified depth. They can also be pre-specified with the `pre` argument

        :return: A dictionary of the following
        `map`: Map of GO terms in the specified sub-ontology to their assigned higher-level terms
        `parents`: Map of chosen parents to their terms
        `unassigned`: GO terms in the sub-ontology graph that are not children of any of the chosen parents. Happens when chosen parents have few children
        `pre`: path to a json file containing pre-defined groups (mapping a GO id or group name to specific terms) that will override other mappings
        """
        children_per_cat = self.metadata.shape[0] / n
        filtered = self.metadata.filter(
            (
                (pl.col("ontology") == ontology)
                & (pl.col("n_children") <= children_per_cat)
                & (pl.col("level") >= min_depth)
            )
        ).sort("n_children", descending=True)
        parents: list = []
        for go_id, *_ in filtered.iter_rows():
            if self.__check_parent(go_id, parents):
                parents.append(go_id)
            if len(parents) == n:
                break
        parent_df = filtered.filter(pl.col("GO_IDs").is_in(parents))
        if show:
            print(parent_df)
        child_to_parent: dict = self.get_predefined(pre, ontology) if pre else {}
        unassigned: set = set(filtered["GO_IDs"])
        for parent in parents:
            for child in self.successors[parent]:
                if child not in child_to_parent:
                    child_to_parent[child] = parent
            child_to_parent[parent] = parent
        unassigned = unassigned - child_to_parent.keys()
        return {
            "map": child_to_parent,
            "parents": dict(zip(parent_df["GO_IDs"], parent_df["term"])),
            "unassigned": unassigned,
        }


def relation_only_paths(paths: list[tuple], relation: str) -> list[list]:
    """Find the path(s) from a list of paths (which are edge lists)
    that consists only of `relation`
    """
    all_results = []
    for path in paths:
        result = []
        if all(map(lambda x: x[2] == relation, path)):
            for p in path:
                result.append((*p[:2][::-1], p[2]))
        # Note: This changes the direction of the obonet GO graph so that successors are children and predecessors are ancestors
        if result:
            all_results.append(result)
    return all_results


def join_dict(df: pl.DataFrame, dct: dict, by: str, colname: str):
    if isinstance(dct, ChainMap):
        dct = dict(dct)
    temp_df = pl.DataFrame(dct).melt().rename({"variable": by, "value": colname})
    return df.join(temp_df, on=by)


def get_successors_nested(complete_GO: nx.MultiDiGraph, term) -> dict:
    """
    Return the complete child hierarchy of a GO term as a nested dictionary of dictionaries
    :param: complete_GO the Networkx object holding the complete GO graph
    :param: term the GO term whose succcessors we want to get
    """
    result = {}
    termlist = nx.bfs_tree(complete_GO, term).nodes
    d = nx.to_dict_of_dicts(complete_GO, nodelist=termlist)
    for adj_term in d.keys():
        if not d[adj_term]:
            d[adj_term] = {}
            continue
        for node in d[adj_term].keys():
            d[adj_term][node] = {}

    result[term] = d[term]

    def rec_helper(cur_dct, root):
        """Recursively reconstruct the nested structure in `result`"""
        current_members: dict = d[root]
        if not current_members:
            return
        for node in list(current_members.keys()):
            if node in d:
                cur_dct[node] = copy.copy(d[node])
                rec_helper(cur_dct[node], node)  # Descends into next dictionary

    for node, dct in result.items():
        rec_helper(dct, node)
    return result


def get_successors(
    G: nx.DiGraph,
    ignore: list = None,
    with_term: bool = False,
    from_node: str = "",
    flatten=False,
) -> dict | set:
    """Return a dictionary mapping the nodes of G to ALL of their
    children/successors (unlike an adjacency list, which
    only returns the immediate) children
    :param: from_node indicates that G is the main GO graph, and the function should
    start finding succesors from the GO term specified by `from_node`.
    :param: flatten Flatten the output entirely to return a set of all the successors
    of G rather than a map
    """
    successors = {}
    if from_node:
        G_old = G
        G = nx.bfs_tree(G, from_node)
        if with_term:
            term_dict = dict(G_old.nodes(data="term", default="NA"))
            nx.set_node_attributes(G, term_dict, name="term")
    to_ignore = set()
    if ignore:
        for ignored_node in ignore:
            try:
                ignored_subtree = nx.bfs_tree(G, ignored_node)
                to_ignore |= set(ignored_subtree.nodes)
            except nx.NetworkXError:
                continue
    if flatten:
        nodes = set(G.nodes) - to_ignore
        if with_term:
            nodes = {f"{c} {G.nodes[c].get('term', 'NA')}" for c in nodes}
        return nodes
    for node in G.nodes:
        if ignore and node in to_ignore:
            continue
        tree = nx.bfs_tree(G, node)  # Creating a bfs tree for each
        # node restricts the view of G to `node` and everything below it
        children = list(tree.nodes)
        if with_term:
            children = [f"{c} {G.nodes[c].get('term', 'NA')}" for c in children]
            node = f"{node} {G.nodes[node].get('term', 'NA')}"
        children.remove(node)
        successors[node] = children
    return successors


def get_level_map(G: nx.DiGraph, root=None) -> dict:
    """Return a dictionary mapping the nodes of G to their levels in G"""
    level_map: dict = {}
    if not (root := G.graph.get("root", root)):
        raise ValueError("Root must be specified!")
    for i, layer in enumerate(nx.bfs_layers(G, root)):
        for node in layer:
            level_map[node] = i
    return level_map


def write_go_metadata(go_obo_path: str, output: str) -> None:
    GO: nx.MultiDiGraph = obonet.read_obo(go_obo_path)
    ontologies = {
        "biological_process": "BP",
        "cellular_component": "CC",
        "molecular_function": "MF",
    }
    meta = {"GO_IDs": [], "term": [], "definition": [], "ontology": []}
    for node, data in GO.nodes.items():
        meta["GO_IDs"].append(node)
        meta["term"].append(data.get("name"))
        meta["definition"].append(data.get("def"))
        meta["ontology"].append(ontologies.get(data.get("namespace")))
    mdf = pl.DataFrame(meta)
    mdf.write_csv(output, separator="\t", null_value="NA")


def to_json(
    go_path: str,
    sample_path: str,
    go_info_path: str,
    outdir: str,
    n_groups: int,
    predefined: str = None,
):
    """Write the GO parent map to a json file"""
    data = {"CC": None, "BP": None, "MF": None}
    S = SubsetGO(go_path, sample_path, go_info_path)
    for sub in data.keys():
        rep = S.get_parents(sub, n=n_groups, pre=predefined)
        data[sub] = {
            "map": rep["map"],
            "unassigned": list(rep["unassigned"]),
        }
    with open(f"{outdir}/go_parents.json", "w") as j:
        json.dump(data, j)


def into_ontology(gos: list, mapping: dict) -> dict:
    partitioned: dict = {"CC": [], "BP": [], "MF": []}
    for go in gos:
        lookup: str = mapping[go]
        if lookup != "NA":  # Check if term is obsolete
            partitioned[lookup].append(go)
    return partitioned


def parent_gos(
    source: str,
    go_terms: list,
    rep_map: dict,
    term_map: dict,
    ontology_map: dict,
    priority: list = None,
) -> dict:
    """Obtain higher level parent GO terms from a term list (i.e. terms of a protein)

    Args:
        source (string): where do these ids come from?
        go_terms (list): list of GO ids to get higher-level terms from
        rep_map (dict): mapping of GO ids to their representative parent terms
        ontology_map (dict): map of GO ids to their sub-ontology
        ontology_map (dict): map of GO ids to their English terms
    """
    grouped = into_ontology(go_terms, ontology_map)
    terms = pl.DataFrame()
    missing = {"CC": 0, "BP": 0, "MF": 0}
    found_special = False
    for o in missing.keys():
        cur_group = grouped[o]
        reduced = {}
        cur_map = rep_map[o]
        for id in cur_group:
            found = cur_map["map"].get(id)
            if found and priority and found in priority:
                reduced[found] = sys.maxsize
                found_special = True
            elif found:
                reduced[found] = reduced.get(found, 0) + 1
            else:
                reduced[id] = 1
        if not reduced:
            continue
        temp = (
            pl.DataFrame(reduced)
            .melt()
            .rename({"variable": "GO_IDs", "value": "count"})
            .with_columns(ontology=pl.lit(o))
            .sample(n=len(reduced), shuffle=True)
            .sort("count", descending=True)
        )
        # Shuffling (using sample) is just a precaution in case all the entries have the same "count"
        terms = pl.concat([terms, temp], how="vertical")
    if terms.is_empty():
        # print(f"Warning: obsolete terms {list(go_terms)}")
        return {"CC": "unknown", "BP": "unknown", "MF": "unknown"}, missing
    aggregated = terms.group_by("ontology", maintain_order=True).agg(
        pl.col("GO_IDs").first()
    )
    result = dict(zip(aggregated["ontology"], aggregated["GO_IDs"]))
    if found_special:
        print(result)
    for k in missing.keys():
        v = result.get(k)
        if not v:
            result[k] = "unknown"
            missing[k] = 1
            continue
        cur_map = rep_map[k]
        if (v not in cur_map["map"] or v in cur_map["unassigned"]) and (
            not priority or v not in priority
        ):
            result[k] = "other"
        else:
            if ":" in v:  # Check that the key is actually a GO id, and
                # not one of the user-defined groups
                result[k] = term_map[v]
            else:
                print(f"from end loop {v}")
                result[k] = v
    return result, missing


def find_go_parents_auto(
    combined_results: str, go_info_path: str, parents_path: str, priority_path: str = ""
) -> tuple[pl.DataFrame, pl.DataFrame]:
    """Find higher-level GO terms of identified proteins in a results file

    Args:
        combined_results (str): path to results file
        go_info_path (str): path to file containing GO term metadata
        parents_path (str): path to json file mapping GO ids to their assigned parent terms
        output (str): output file name
    """
    info: pl.DataFrame = pl.read_csv(go_info_path, separator="\t")
    id2ontology = dict(zip(info["GO_IDs"], info["ontology"]))
    id2term = dict(zip(info["GO_IDs"], info["term"]))
    data = pl.read_csv(combined_results, separator="\t", null_values="NA")
    shape_before: tuple = data.shape
    has_go = data.filter(pl.col("GO_IDs").is_not_null())
    no_go = data.filter(pl.col("GO_IDs").is_null())
    if priority_path:
        with open(priority_path, "r") as j:
            priority = list(json.load(j).keys())
    else:
        priority = []
    with open(parents_path, "r") as j:
        rmap: dict = json.load(j)
    id_go: pl.DataFrame = (
        has_go.select("ProteinId", "GO_IDs", "GO_counts")
        .with_columns(
            GO_IDs=pl.col("GO_IDs").map_elements(
                lambda x: x.split(";"), return_dtype=pl.List(pl.Utf8)
            )
        )
        .sort("GO_counts", descending=True)
    )
    cc, bp, mf = [], [], []
    missing_counts = {"CC": 0, "BP": 0, "MF": 0}
    for id, gos in zip(id_go["ProteinId"], id_go["GO_IDs"]):
        parents, missing = parent_gos(
            id, gos, rmap, id2term, id2ontology, priority=priority
        )
        cc.append(parents["CC"])
        bp.append(parents["BP"])
        mf.append(parents["MF"])
        for k, v in missing.items():
            missing_counts[k] += v
    added = has_go.with_columns(
        GO_category_CC=pl.Series(cc),
        GO_category_MF=pl.Series(mf),
        GO_category_BP=pl.Series(bp),
    )
    no_go = no_go.with_columns(
        GO_category_CC=pl.lit("unknown"),
        GO_category_MF=pl.lit("unknown"),
        GO_category_BP=pl.lit("unknown"),
    )
    result = pl.concat([added, no_go], how="vertical")
    missing_df = (
        pl.DataFrame(missing_counts)
        .melt()
        .rename({"variable": "ontology", "value": "n_missing"})
    )
    if not result.shape[0] == shape_before[0]:
        raise ValueError(
            f"""
                        nrows before and after do not match!\n
                        rows before: {shape_before[0]}\n
                        rows current: {result.shape[0]}\n
                         """
        )
    return result.to_pandas(), missing_df.to_pandas()


class CompleteGO(SubsetGO):
    def __init__(
        self,
        go_path: str,
        metadata_path: str,
        get_metadata: bool = False,
    ) -> None:
        metadata = pl.read_csv(metadata_path, separator="\t")
        self.roots = {"BP": "GO:0008150", "CC": "GO:0005575", "MF": "GO:0003674"}
        self.root_map: dict = dict(zip(metadata["GO_IDs"], metadata["ontology"]))
        self.term_map: dict = dict(zip(metadata["GO_IDs"], metadata["term"]))
        if ".obo" in go_path:
            GO: nx.MultiDiGraph = obonet.read_obo(go_path)
            self.G: nx.MultiDiGraph = nx.MultiDiGraph()
            self.G.add_nodes_from(self.roots.values())
            for go in metadata["GO_IDs"]:
                if go in GO:
                    paths: list = nx.all_simple_edge_paths(
                        GO, source=go, target=self.roots[self.root_map[go]]
                    )
                    is_a = relation_only_paths(paths, "is_a")
                    if is_a:
                        for i in is_a:
                            self.G.add_edges_from(i)
            nx.set_node_attributes(self.G, self.term_map, name="term")
        elif ".gml" in go_path:
            self.G = nx.read_gml(go_path)
        self.definition_map = dict(zip(metadata["GO_IDs"], metadata["definition"]))
        self.metadata = metadata
        self.successors = None
        self.GL, self.name_map = self.w_term_label()
        if get_metadata:
            self.get_metadata()

    def get_metadata(self):
        self.metadata = self.get_node_data().join(
            self.metadata.select("GO_IDs", "term", "definition", "ontology"),
            on="GO_IDs",
        )

    def w_term_label(self) -> tuple[nx.MultiDiGraph, dict]:
        label_map = {}
        for t in self.term_map:
            label_map[t] = f"{t} {self.term_map[t]}"
        return nx.relabel_nodes(self.G, label_map), label_map

    def browse(
        self,
        regex: str,
        as_df: bool = False,
        ontology: str = None,
        sort: bool = True,
        with_label: bool = True,
    ) -> dict | pl.DataFrame:
        filtered = self.metadata.filter(
            pl.any_horizontal(pl.col(["term", "definition"]).str.contains(regex))
        )
        if ontology and ontology in {"BP", "MF", "CC"}:
            filtered = filtered.filter(pl.col("ontology") == ontology)
        if sort:
            filtered = filtered.sort(
                by=["n_children", "level"], descending=[True, False]
            )
        if as_df:
            return filtered
        result = {}
        if "level" in filtered.columns and "n_children" in filtered.columns:
            for row in filtered.iter_rows(named=True):
                result[row["GO_IDs"]] = {
                    "ontology": row["ontology"],
                    "term": row["term"],
                    "definition": row["definition"],
                    "level": row["level"],
                    "n_children": row["n_children"],
                }
        else:
            for row in filtered.iter_rows():
                result[row[0]] = {
                    "ontology": row[3],
                    "term": row[1],
                    "definition": row[2],
                }
        if sort:
            ordered = []
            for g in filtered["GO_IDs"]:
                metadata = {"M": result[g]}
                if with_label:
                    o = result[g].get("ontology")
                    del result[g]["term"]
                    metadata["I"] = f'"{g}|{o} {self.term_map.get(g)}"'
                else:
                    metadata["id"] = f'"{g}"'
                ordered.append(metadata)
            return ordered
        else:
            return result

    def nested_successors(self, term: str, with_label=True):
        if with_label:
            return get_successors_nested(self.GL, self.name_map[term])
        return get_successors_nested(self.G, term)

    def get_parents(self, term: str, with_label: bool = True):
        def parse_edge_list(edges: list):
            result = []
            for index, edge in enumerate(edges):
                if with_label:
                    t = f"{edge[0]} {self.term_map.get(edge[0])}"
                else:
                    t = edge[0]
                t = f"{' ' * index} {t}"
                result.append(t)
            if with_label:
                result.append(f"{' ' * len(edges)} {term} {self.term_map.get(term)}")
            else:
                result.append(f"{' ' * len(edges)} {term}")
            return result

        root = self.roots.get(self.root_map.get(term))
        if root:
            paths = nx.all_simple_edge_paths(self.G, source=root, target=term)
            return [parse_edge_list(p) for p in paths]

    def get_node_data(self):
        self.successors: dict = {}
        level_map: dict = {}
        for root in self.roots.values():
            current = nx.bfs_tree(self.G, root)
            current.graph["root"] = root
            self.successors = ChainMap(self.successors, get_successors(current))
            level_map = ChainMap(level_map, get_level_map(current))

        return (
            (pl.DataFrame({"GO_IDs": list(self.G.nodes)}))
            .with_columns(
                n_children=pl.col("GO_IDs").map_elements(
                    lambda x: len(get_successors(self.G, from_node=x)),
                    return_dtype=pl.Int16,
                ),
                level=pl.col("GO_IDs").map_elements(
                    lambda x: level_map[x], return_dtype=pl.Int16
                ),
            )
            .sort("n_children", descending=True)
        )


# TODO: Check that this works
def create_group_mapping(
    go_path: str,
    metadata_path: str,
    group_mapping_file: str,
    output: str,
    group_name=None,
) -> None:
    """Find all the relevant child GO terms to the terms specified in `group_mapping_file`
    :param: group_mapping_file path to a toml file that maps user-specified groups to a list of `seed` GO terms. For each term, all of its children will be added to the list, UNLESS that GO term has already been specified as a `seed` term in the original file
    :param: group_name the header of the toml table containing the desired group
    :return: None
    """
    C = CompleteGO(go_path=go_path, metadata_path=metadata_path)
    data_dict = {"GO_IDs": [], "Group": []}
    with open(group_mapping_file, "rb") as f:
        mapping = tomllib.load(f)
    if group_name:
        mapping = mapping[group_name]
    for group, parents in mapping.items():
        for parent in parents:
            data_dict["GO_IDs"].append(parent)
            data_dict["Group"].append(group)
            for child in C.successors.get(parent, []):
                data_dict["GO_IDs"].append(child)
                data_dict["GO_IDs"].append(group)
    pl.DataFrame(data_dict).write_csv(output, separator="\t")
