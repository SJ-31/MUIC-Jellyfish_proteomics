import sys
import re
import pandas as pd
import polars as pl
import intervaltree as it
from thefuzz import process


def group_peptide_origin(id_series: pl.Series) -> pl.DataFrame:
    def get_row(peptide_ids: set):
        denovo, transcriptome, unmatched = set(), set(), set()
        if peptide_ids:
            for p in peptide_ids:
                if "D" in p:
                    denovo.add(p)
                elif "T" in p:
                    transcriptome.add(p)
                else:
                    unmatched.add(p)
        return (denovo, transcriptome, unmatched)

    tmp = {
        "MatchedDenovoIds": [],
        "MatchedTranscriptomeIds": [],
        "MatchedUPeptideIds": [],
    }
    for r in map(get_row, id_series):
        tmp["MatchedDenovoIds"].append(r[0])
        tmp["MatchedTranscriptomeIds"].append(r[1])
        tmp["MatchedUPeptideIds"].append(r[2])
    return pl.DataFrame(tmp)


def clean_peptide(peptide):
    if len(set(peptide) & {"]"}) > 0:
        return "".join(re.findall("[A-Z]+", peptide))
    return peptide


def parent_dir(file_str, parent_level=-1):
    if parent_level > 0:
        raise ValueError("The parent level must be negative!")
    dirs = file_str.split("/")
    return "/".join(dirs[:parent_level])


sys.path.append(parent_dir(__file__))
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
        self.alignments: pl.DataFrame = pl.from_pandas(alignments)
        self.alignments = self.alignments.with_columns(
            interval=pl.col("alignment").map_elements(
                lambda x: (findChar(x), findChar(x, True)),
                return_dtype=pl.List(pl.Int64),
            ),
        )
        # Retrieve sequence of alignments i.e. ---ACTAGAT--- -> ACTAGAT
        self.alignments = self.alignments.with_columns(
            length=pl.col("interval").map_elements(
                lambda x: x[1] - x[0], return_dtype=pl.Int16
            ),
            seq=pl.struct("alignment", "interval").map_elements(
                lambda x: x["alignment"][x["interval"][0] : x["interval"][1]],
                return_dtype=pl.String,
            ),
        )
        # Get matched peptides from combined results
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

    def __get_alignment_interval(self, sequence: str) -> list:
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

    def filter_protein_id(self, protein_id: str) -> None:
        """
        Filter main dataframe onto given protein
        """
        self.cur_alignments = self.alignments.filter(pl.col("ProteinId") == protein_id)
        find_matched = self.matched_peptides.filter(pl.col("ProteinId") == protein_id)
        # Necessary when using aligning proteins from engines,
        # as the engines' matched proteins could have been aligned back into the given
        # protein
        if find_matched.is_empty():
            self.cur_protein = self.proteins.filter(pl.col("ProteinId") == protein_id)
        else:
            found = find_matched[0].select("MatchedPeptideIds").item()
            self.cur_protein = self.proteins.filter(
                (pl.col("ProteinId") == protein_id) | (pl.col("ProteinId").is_in(found))
            )

    def engine_alignment_intervals(
        self, protein_id: str, engine: str
    ) -> list[list[int]]:
        self.filter_protein_id(protein_id)
        temp_prot = self.cur_protein.filter(pl.col("engine") == engine)
        tree = it.IntervalTree()
        for seq in temp_prot["peptideIds"]:
            begin, end = self.__get_alignment_interval(seq)
            if begin is not None:
                tree[begin : end + 1] = seq
        tree.merge_overlaps()
        return [list(v)[:-1] for v in tree]


class DenovoAlignmentTracer(AlignmentTracer):
    def __init__(
        self,
        proteins: pd.DataFrame,
        alignments: pd.DataFrame,
        matched_peptides: pd.DataFrame,
        seq_map_path=None,
        unmatched_path=None,
    ) -> None:
        super().__init__(proteins, alignments, matched_peptides)
        grouped = group_peptide_origin(self.matched_peptides["MatchedPeptideIds"])
        self.matched_peptides = pl.concat(
            [
                self.matched_peptides.select("ProteinId", "MatchedPeptideIds"),
                grouped,
            ],
            how="horizontal",
        )
        split_list = list(
            map(
                lambda x: x.split(";"), self.proteins["MatchedPeptideIds"].drop_nulls()
            ),
        )
        wanted = list({s for split in split_list for s in split})
        self.seq_map: pl.DataFrame = get_seq_map(seq_map_path, unmatched_path, wanted)
        self.seq_dict: dict = dict(zip(self.seq_map["seq"], self.seq_map["id"]))
        self.seqs = list(self.seq_dict.keys())
        self.alignments = self.alignments.join(self.seq_map, on="seq", how="left")
        id_df = self.alignments.select(
            pl.struct("id", "seq").map_elements(
                lambda x: self.find_closest_seq(x["id"], x["seq"]),
                return_dtype=pl.String,
            )
        )
        # TODO: Takes ages, so save to some hidden file
        self.alignments = self.alignments.with_columns(id=id_df["id"])

    def find_closest_seq(self, id, seq):
        if id:
            return id
        for match in process.extract(seq, self.seqs):
            if match[1] >= 90:
                return self.seq_dict[match[0]]
        return pl.Null

    # def denovo_alignment_intervals(self, protein_id: str):
    #     m = self.matched_peptides.filter(pl.col("ProteinId") == protein_id)
    #     coverage_dct = {"Denovo": 0, "Transcriptome": 0, "UPeptide": 0, "undetermined": 0}
    #     for origin in ["Denovo", "Transcriptome", "UPeptide"]:
    #         for id in m[f"Matched{origin}Ids"]:
    #             for
    # TODO: finish this
    # return coverage_dct


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


def get_interval_groups(interval_dict, id):
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


def substitution_type(type_str: str, target: str = "conservative") -> int:
    splits = type_str.split("->")
    if splits[0] == splits[1] and target == "conservative":
        return 1
    if splits[0] != splits[1] and target == "non_conservative":
        return 1
    return 0


def aggregate_mismatches(
    df: pl.DataFrame, grouping_col: str = "ProteinId", average=False
) -> pl.DataFrame:
    """Aggregate useful stats for mismatches"""

    def get_average(expr):
        return expr / pl.col("ProteinId").unique().len()

    exprs = {
        "mismatches": pl.len(),
        "disagreements": pl.col("index").is_duplicated().sum() / 2,
        "cons": pl.sum("conservative"),
        "n_cons": pl.sum("non_conservative"),
    }
    if average:
        for k, v in exprs.items():
            exprs[k] = get_average(v)
    mismatch_metrics = (
        df.group_by(grouping_col)
        .agg(
            n_mismatches=exprs["mismatches"],
            n_disagreements=exprs["disagreements"],
            n_conservative=exprs["cons"],
            n_non_conservative=exprs["n_cons"],
        )
        .with_columns(
            nc_c_ratio=pl.col("n_non_conservative") / pl.col("n_conservative")
        )
    )
    return mismatch_metrics


def classify_mismatches(df: pl.DataFrame) -> pl.DataFrame:
    """Classify amino acid mismatches as conservative or non-conservative"""
    if isinstance(df, pd.DataFrame):
        df = pl.from_pandas(df)
    return df.with_columns(
        conservative=pl.col("type").map_elements(
            lambda x: substitution_type(x, "conservative"),
            return_dtype=pl.Int16,
        )
    ).with_columns(
        non_conservative=pl.col("conservative").map_elements(
            lambda x: 1 if x == 0 else 0,
            return_dtype=pl.Int16,
        )
    )


def seq_from_alignment(alignment: str) -> str:
    start, stop = findChar(alignment), findChar(alignment, True)
    return alignment[start : stop + 1]


def denovo_mismatch_metrics(
    ids_to_keep: list,
    ids_to_keep_denovo: list,
    seq_map_path: str,
    unmatched_path: str,
    aligned_peptides: pd.DataFrame,
    mismatches: pd.DataFrame,
):
    """Get metrics for the amino acid mismatches occuring on de novo, transcriptome and unknown proteins
    :param: ids_to_keep List of ProteinIds for proteins that were mapped by de novo peptides
    :param: ids_to_keep_denovo List of ProteinIds for the de novo proteins mapping to the main proteins
    :return: A dictionary with two entries:
        "map": Mapping of denovo peptides to their protein ids for convenience
        "metrics": metrics of denovo peptides' amino acid mismatches
    """
    ids_to_keep, ids_to_keep_denovo = pl.Series(ids_to_keep), pl.Series(
        ids_to_keep_denovo
    )
    mismatches = pl.from_pandas(mismatches)
    id_map = get_seq_map(seq_map_path, unmatched_path, ids_to_keep_denovo)
    peptides = (
        pl.from_pandas(aligned_peptides)
        .filter(pl.col("ProteinId").is_in(ids_to_keep))
        .with_columns(
            sequence=pl.col("alignment").map_elements(
                seq_from_alignment, return_dtype=pl.String
            ),
            start=pl.col("alignment").map_elements(
                lambda x: findChar(x), return_dtype=pl.Int32
            ),
            stop=pl.col("alignment").map_elements(
                lambda x: findChar(x, True), return_dtype=pl.Int32
            ),
        )  # Obtain sequence from alignment i.e. ---AGCMNAR--- AGCMNAR
        # And try to map this with original
    ).select(pl.col("*").exclude("UniProtKB_ID", "alignment"))
    denovo_map = peptides.join(id_map, left_on="sequence", right_on="seq", how="inner")

    metrics = (
        denovo_map.join(mismatches, on="ProteinId")
        .filter(
            (pl.col("start") <= pl.col("index")) & (pl.col("index") <= pl.col("stop"))
        )
        .pipe(classify_mismatches)
        .pipe(lambda x: aggregate_mismatches(x, "id", True))
    )
    return {"mapping": denovo_map, "metrics": metrics}


def get_seq_map(file_path: str, unmatched_path: str, ids_to_keep: list) -> pl.DataFrame:
    map_to = pl.Series(ids_to_keep)
    id_map = (
        pl.read_csv(file_path, separator="\t")
        .filter(pl.col("id").is_in(map_to))
        .select(pl.col("*").exclude("header", "mass", "length"))
    )
    unmatched = (
        (
            pl.read_csv(unmatched_path, separator="\t", null_values="NA")
            .with_columns(
                seq=pl.col("peptideIds").map_elements(
                    clean_peptide, return_dtype=pl.String
                )
            )
            .rename({"ProteinId": "id"})
        )
        .select("id", "seq")
        .filter(pl.col("id").is_in(map_to))
    )
    return pl.concat([id_map, unmatched])
