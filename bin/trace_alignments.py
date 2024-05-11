import sys
import pandas as pd
import polars as pl
import intervaltree as it
from thefuzz import process


def parentDir(file_str, parent_level=-1):
    if parent_level > 0:
        raise ValueError("The parent level must be negative!")
    dirs = file_str.split("/")
    return "/".join(dirs[:parent_level])


sys.path.append(parentDir(__file__))
from view_alignments import findChar


class AlignmentTracer:
    """Class for retracing how the peptides from a given engine
    were aligned onto the full-length protein. Methods returns the interval list that the engines' peptides cover on the protein
    """

    def __init__(
        self,
        proteins: pd.DataFrame,
        alignments: pd.DataFrame,
        matched_peptides: pd.DataFrame,
    ) -> None:
        """
        :param proteins: dataframe (as provided from R) containing protein identifications from percolator with an "engine" column specifying which engine made each identification. Must
        be filtered to be for a single ProteinId
        :param alignments: dataframe containing alignments (output of
        the `coverage_calc` process, filtered for the given ProteinId
        """
        self.proteins: pl.DataFrame = pl.from_pandas(proteins)
        self.alignments = pl.DataFrame = pl.from_pandas(alignments)
        self.alignments = self.alignments.with_columns(
            interval=pl.col("alignment").map_elements(
                lambda x: (findChar(x), findChar(x, True)),
                return_dtype=pl.List(pl.Int64),
            ),
        )
        self.alignments = self.alignments.with_columns(
            length=pl.col("interval").map_elements(
                lambda x: x[1] - x[0], return_dtype=pl.Int16
            ),
            seq=pl.struct("alignment", "interval").map_elements(
                lambda x: x["alignment"][x["interval"][0] : x["interval"][1]],
                return_dtype=pl.String,
            ),
        )
        self.matched_peptides = (
            pl.from_pandas(matched_peptides)
            .filter(pl.col("MatchedPeptideIds").is_not_null())
            .with_columns(
                MatchedPeptideIds=pl.col("MatchedPeptideIds").map_elements(
                    lambda x: set(x.split(";")), return_dtype=pl.Object
                )
            )
        )
        self.cur_protein: pl.DataFrame | None = None
        self.cur_alignments: pl.DataFrame | None = None

    def __getAlignmentInterval(self, sequence: str) -> list:
        found = None
        for match in process.extract(sequence, self.cur_alignments["seq"]):
            if match[1] >= 90:
                found = match[0]
                break
        if found:
            return self.cur_alignments.filter(pl.col("seq") == found)["interval"][
                0
            ].to_list()
        return [None, None]

    def __filterProteinId(self, protein_id: str) -> None:
        self.cur_alignments = self.alignments.filter(pl.col("ProteinId") == protein_id)
        find_matched = self.matched_peptides.filter(pl.col("ProteinId") == protein_id)
        if find_matched.is_empty():
            self.cur_protein = self.proteins.filter(pl.col("ProteinId") == protein_id)
        else:
            found = find_matched[0].select("MatchedPeptideIds").item()
            self.cur_protein = self.proteins.filter(
                (pl.col("ProteinId") == protein_id) | (pl.col("ProteinId").is_in(found))
            )

    def engineAlignmentIntervals(self, protein_id: str, engine: str) -> list[list[int]]:
        self.__filterProteinId(protein_id)
        temp_prot = self.cur_protein.filter(pl.col("engine") == engine)
        tree = it.IntervalTree()
        for seq in temp_prot["peptideIds"]:
            begin, end = self.__getAlignmentInterval(seq)
            if begin is not None:
                tree[begin : end + 1] = seq
        tree.merge_overlaps()
        return [list(v)[:-1] for v in tree]


ENGINE_LIST = [
    "comet",
    "identipy",
    "metamorpheus",
    "msfragger",
    "msgf",
    "tide",
    "metamorpheusGTPMD",
    "msfraggerGPTMD",
    "msfraggerGlyco",
]


def getIntervalGroups(interval_dict, id):
    tt = it.IntervalTree()
    for engine, interval_list in interval_dict.items():
        if interval_list:
            for interval in interval_list:
                tt[interval[0] : interval[1]] = engine
    tt.split_overlaps()
    tt.merge_equals(data_reducer=lambda x, y: f"{x};{y}")
    tt = sorted(tt)
    info = {e: [] for e in ENGINE_LIST}
    for index, interval in enumerate(tt):
        for engine in ENGINE_LIST:
            if engine in set(interval.data.split(";")):
                info[engine].append(str(index + 1))
    for key in info.keys():
        if not info[key]:
            info[key] = pd.NA
        else:
            info[key] = ";".join(info[key])
    info["n_unique"] = len(tt)
    return pd.DataFrame(info, index=[id])
