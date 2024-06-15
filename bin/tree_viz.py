from collections import ChainMap
from collections.abc import Callable
import polars as pl
import ete4 as et
from matplotlib.colors import rgb2hex
from matplotlib import colormaps
import ete4.treeview.faces as tvf
import ete4.treeview as tv
from enum import Enum

FONT = "Arial"


class Colormaps(Enum):
    LEAVES = colormaps.get_cmap("Oranges")
    PHYLUM = colormaps.get_cmap("Purples")
    RANKS = colormaps.get_cmap("Set1")


def dfToDict(df: pl.DataFrame, keys: str = None, values: str = None) -> dict:
    """Convert polars df into a dictionary. By default, takes the first and
    second columns as keys and values respectively
    """
    if not keys and not values:
        cols = df.columns
        return dict(zip(df[cols[0]], df[cols[1]]))
    return dict(zip(df[keys], df[values]))


RANKS: dict = {
    "kingdom": 0,
    "phylum": 1,
    "class": 2,
    "order": 3,
    "family": 4,
    "genus": 5,
}


def fillNodeStyle(node: et.PhyloTree, **kwargs):
    for k, v in kwargs.items():
        node.img_style[k] = v


def getRankColor(rank: str):
    return rgb2hex(Colormaps.RANKS.value(RANKS[rank]))


def getRankSize(rank: str):
    return RANKS[rank] + 10 - RANKS[rank]


def rankStyle(node: et.PhyloTree, rank: str):
    rgb: tuple = getRankColor(rank)
    color: str = rgb2hex(rgb)
    size = getRankSize(rank)
    fillNodeStyle(node, size=size, fgcolor=color, shape="square")


class TaxaTree:
    def __init__(self, taxa: pl.DataFrame | str, taxdump_path: str) -> None:
        if isinstance(taxa, str):
            taxa = pl.read_csv(taxa, separator="\t", null_values="NA")
        self.NCBI = et.NCBITaxa(taxdump_file=taxdump_path)
        counts = []
        percents = []
        for r in RANKS:
            r = r.capitalize()
            count_df = taxa[r].value_counts().filter(pl.col(r).is_not_null())
            percent_df = count_df.with_columns(
                count=pl.col("count") / pl.sum("count")
            ).rename({"count": "percent"})
            counts.append(dfToDict(count_df))
            percents.append(dfToDict(percent_df))
        all_rank_percents: ChainMap = ChainMap(*percents)
        all_rank_counts: ChainMap = ChainMap(*counts)

        complete_tree, name2taxid = self.tree(
            all_rank_percents.keys(), "Species", prune_others=True
        )
        all_nodes: set = {d.name for d in complete_tree.traverse()}
        taxid2name = {v: k for k, v in name2taxid.items()}
        for name, percent in all_rank_percents.items():
            taxid = name2taxid.get(name)
            count = all_rank_counts[name]
            if taxid and taxid in all_nodes:
                complete_tree[taxid].add_prop("weight", percent)
                complete_tree[taxid].add_prop("n", count)
                complete_tree[taxid].del_prop("taxid")
        self.T = complete_tree
        self.taxid2name = taxid2name

    def tree(
        self, taxa: set[str], rank_limit: str = "Order", prune_others: bool = False
    ) -> tuple[et.PhyloTree, dict]:
        id_map = {k: str(v[0]) for k, v in self.NCBI.get_name_translator(taxa).items()}
        T = self.NCBI.get_topology(
            id_map.values(),
            rank_limit=rank_limit.lower(),
            annotate=True,
        )
        if prune_others:
            others = [
                n for n in T.traverse() if n.props["rank"].capitalize() not in RANKS
            ]
            map(lambda x: x.delete(), others)
        return T, id_map

    def getSubtree(
        self, rank="", taxid: int = None, sci_name: str = ""
    ) -> et.PhyloTree:
        T = self.T.copy()
        if taxid or sci_name:
            for node in T.traverse():
                if node.props["sci_name"] == sci_name or node.props["name"] == str(
                    taxid
                ):
                    T = node.detach()
        if rank:
            keepRanksAbove(T, rank)
        return T

    def findTaxon(self, rank="genus", **props) -> et.PhyloTree:
        subtree: et.PhyloTree = list(self.T.search_nodes(**props))[0].copy()
        if rank:
            keepRanksAbove(subtree, rank)
        return subtree


def keepRanksAbove(T: et.PhyloTree, rank: str):
    """Keep all taxonomic ranks in T that are above `rank`.
    Removes all ranks lower than or equal to `rank`"""
    rank_num: int = RANKS[rank]
    kept_ranks = [k for k, v in RANKS.items() if v <= rank_num]
    kept_nodes = {n for n in T.traverse() if n.props.get("rank") in kept_ranks}
    T.prune(kept_nodes)


def defaultLayout(node):
    def addWeight(cm):
        if weight := props.get("weight"):
            node.img_style["bgcolor"] = rgb2hex(cm(weight))
            return weight
        return 0

    props: dict = node.props
    if props.get("legend"):
        return
    sci_name = props["sci_name"]
    if (count := props.get("n")) and node.is_leaf:
        label_text: str = f"{sci_name}, n = {count}"
    else:
        label_text: str = props["sci_name"]
    if (rank := props.get("rank")) and rank != "no rank":
        if rank in RANKS:
            rankStyle(node, rank)
        else:
            node.img_style["shape"] = "sphere"
            node.img_style["fgcolor"] = "black"
            other_clade = f"Rank: {rank}"
            tvf.add_face_to_node(
                TextFace(
                    other_clade, penwidth=100, margin_left=5, margin_right=3, bold=True
                ),
                node,
                column=0,
            )
    label_face = TextFace(
        label_text, fsize=10, ftype=FONT, margin_left=5, margin_right=3
    )
    if node.is_leaf:
        weight = addWeight(Colormaps.LEAVES.value)
        if weight > 0.7:
            label_face.fgcolor = "white"
    tvf.add_face_to_node(label_face, node, 0)


def show(
    tree: et.PhyloTree,
    circular=False,
    layout=defaultLayout,
    legendFun: Callable[[tv.FaceContainer], None] = None,
    arc_start=0,
    arc_span=0,
    save_to="",
    save_params: dict = None,
) -> None:
    style = tv.TreeStyle()
    if legendFun:
        legendFun(style.title)  # Fits better at this position
    style.layout_fn = layout
    style.show_scale = False
    style.show_leaf_name = False
    style.draw_guiding_lines = True
    style.draw_aligned_faces_as_table = True
    style.show_border = True
    if circular:
        style.mode = "c"
        if arc_start:
            style.arc_start = arc_start
        if arc_span:
            style.arc_span = arc_span
    else:
        style.branch_vertical_margin = 5
    if save_to:
        tree.render(save_to, tree_style=style, **save_params)
    else:
        tree.show(tree_style=style)


def TextFace(text: str, **kwargs) -> tvf.TextFace:
    face = tvf.TextFace(text)
    for k, v in kwargs.items():
        setattr(face, k, v)
    return face


def rankLegend(legend: tv.FaceContainer, max_rank: str) -> None:
    left_margin, top_margin = 20, 15
    header = TextFace(
        "Key to main taxonomic ranks",
        fsize=15,
        margin_left=left_margin,
        margin_top=top_margin,
        bold=True,
    )
    legend.add_face(header, column=1)
    for rank in RANKS:
        size = getRankSize(rank)
        item = TextFace(
            f"{rank.capitalize()}",
            fsize=size,
            fgcolor=getRankColor(rank),
            ftype=FONT,
            margin_left=left_margin + 5,
        )
        legend.add_face(item, column=1)
        if rank == max_rank:
            break
    return
