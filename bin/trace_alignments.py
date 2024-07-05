import sys
from Bio.motifs.matrix import GenericPositionMatrix
from Bio.SeqRecord import SeqRecord
from Bio import Align
import math
from Bio.Seq import Seq
from Bio import motifs
import numpy as np
import re
import pandas as pd
import polars as pl
import intervaltree as it

AMINO_ACIDS = "ACDEFGHIKLMNPQRSTVWYUO"


def parent_dir(file_str, parent_level=-1):
    if parent_level > 0:
        raise ValueError("The parent level must be negative!")
    dirs = file_str.split("/")
    return "/".join(dirs[:parent_level])


sys.path.append(parent_dir(__file__))
import view_alignments as va


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
from view_alignments import find_char


def categorize_by_id(id):
    if "T" in id:
        return "transcriptome"
    elif "P" in id:
        return "database"
    elif "D" in id:
        return "denovo"
    return "unmatched_peptide"


class AlignmentTracer:
    """Class for retracing how the peptides from a given engine
    were aligned onto the full-length protein. Methods returns the interval list that the engines' peptides cover on the protein
    """

    def __init__(
        self,
        alignment_path: str,
        peptide_map_path: str,
    ) -> None:
        """
        :param alignment_path: path to output of the `coverage_calc` process
        :param peptide_map_path: path to `percolator_peptide_map` file, file mapping
        identified peptides in the pipeline to the engines which identified them
        """
        alignments: pl.DataFrame = pl.read_csv(
            alignment_path, separator="\t", null_values="NA"
        )
        unknown = alignments.filter(pl.col("id") == "unknown")
        print(
            f"Proportion of unidentified alignments: {unknown.shape[0]/alignments.shape[0] * 100}"
        )
        self.alignments = alignments.with_columns(
            length=pl.col("end") - pl.col("start"),
            interval=pl.Series(zip(alignments["start"], alignments["end"])),
        )
        self.peptide_map: pl.DataFrame = pl.read_csv(
            peptide_map_path, separator="\t", null_values="NA"
        )
        temp = (
            self.peptide_map.filter(pl.col("engine").is_not_null())
            .unique(subset=["ProteinId", "engine"])
            .group_by("ProteinId")
            .agg(pl.col("engine"))
        )
        self.id2engines = dict(
            zip(list(temp["ProteinId"]), [list(set(x)) for x in temp["engine"]])
        )
        self.alignment_types = [
            "denovo",
            "database",
            "transcriptome",
            "unmatched_peptide",
        ]
        self.engines: list = list(self.peptide_map["engine"].unique().drop_nulls())
        self.cur_protein: pl.DataFrame | None = None
        self.cur_alignments: pl.DataFrame | None = None
        self.id2metadata: dict = {}
        self.temp_data: dict = {
            "ProteinId": [],
            "unidentified_coverage": [],
            "unidentified_count": [],
            "total_count": [],
        }
        for lst in [self.alignment_types, self.engines]:
            for val in lst:
                for metric in ["count", "coverage"]:
                    self.temp_data[f"{val}_{metric}"] = []

    def get_id2seq(self, data_path: str) -> None:
        df = pl.read_csv(data_path, separator="\t", null_values="NA")
        meta = [(s, h) for s, h in zip(df["seq"], df["header"])]
        self.id2metadata = dict(zip(df["ProteinId"], meta))

    def filter_protein_id(self, protein_id: str) -> pl.DataFrame:
        """
        Filter main dataframe onto given protein
        """
        df = self.alignments.filter(pl.col("ProteinId") == protein_id)
        df = df.with_columns(
            engine_matches=pl.col("id").map_elements(
                lambda x: self.id2engines.get(x, []),
                return_dtype=pl.List(pl.String),
            ),
            alignment_type=pl.col("id").map_elements(
                lambda x: categorize_by_id(x), return_dtype=pl.String
            ),
        )
        # TODO Find ids and unmatched peptides better
        # unmatched_peptides =

        return df

    def get_coverage_contributions(self, protein_id: str, df: pl.DataFrame) -> dict:
        """
        Retrieves the coverage (and aligned peptide count) contributions of the different engines and aligned peptide types for the given protein
        :param: df Dataframe resulting from calling `filter_protein_id` on
        `protein_id`
        """
        metrics: dict = {}
        for engine in self.engines:
            current = df.filter(pl.col("engine_matches").list.contains(engine))
            metrics[f"{engine}_coverage"] = coverage_helper(current)
            metrics[f"{engine}_count"] = current.shape[0]
        for type in self.alignment_types:
            current = df.filter(pl.col("alignment_type") == type)
            metrics[f"{type}_coverage"] = coverage_helper(current)
            metrics[f"{type}_count"] = current.shape[0]
        unidentified = df.filter(
            (pl.col("engine_matches").list.len() == 0) | (pl.col("id") == "unknown")
        )
        metrics["unidentified_coverage"] = coverage_helper(unidentified)
        metrics["unidentified_count"] = unidentified.shape[0]
        for k in self.temp_data:
            if k == "ProteinId":
                self.temp_data[k].append(protein_id)
            elif k == "total_count":
                self.temp_data[k].append(df.shape[0])
            else:
                self.temp_data[k].append(metrics[k])
        return metrics

    def run(self) -> pl.DataFrame:
        for id in self.alignments["ProteinId"].unique():
            df = self.filter_protein_id(id)
            _ = self.get_coverage_contributions(id, df)
        return pl.DataFrame(self.temp_data)

    def plot_engines_alignment(self, protein_id: str, outdir: str) -> None:
        """Obtain the consensus sequnces of each engines' peptides
        for `protein_id`, creating a visualization of the peptides
        aligned onto the sequence of `protein_id`
        """
        if not self.id2metadata:
            raise ValueError("Must initialize id2seq mapping first!")
        filtered_df: pl.DataFrame = self.filter_protein_id(protein_id)
        seq_list = [str(get_engine_consensus(filtered_df, e)) for e in self.engines]
        seqs = into_msa(dict(zip(self.engines, seq_list)))
        seq, header = self.id2metadata[protein_id]
        cur_seq: SeqRecord = SeqRecord(seq=seq, id=header)
        msa = va.PeptideViz(
            seqs,
            wrap_length=90,
            show_consensus=True,
            aligned_to=cur_seq,
        )
        msa.set_custom_color_scheme(va.COLOR_SCHEME)
        # fig = msa.plotfig()
        outfile = f"{outdir}/{protein_id}.png"
        msa.savefig(outfile)


def get_engine_consensus(df, engine) -> Seq:
    """Convert all the alignments produced by the given engine into `filename`"""
    cur_engine = df.filter(pl.col("engine_matches").list.contains(engine))
    motif = motifs.create(cur_engine["alignment"], alphabet=AMINO_ACIDS)
    return get_consensus(motif.counts)


def into_msa(header2seq: dict[str, str]) -> Align.MultipleSeqAlignment:
    seqs = [SeqRecord(seq=v, id=k) for k, v in header2seq.items()]
    return Align.MultipleSeqAlignment(seqs)


def get_consensus(m: GenericPositionMatrix):
    if set(m.alphabet).union("ACGTUN-") == set("ACGTUN-"):
        undefined = "N"
    else:
        undefined = "X"
    sequence = ""
    for i in range(m.length):
        maximum = -math.inf
        current_letter = ""
        for letter in m.alphabet:
            count = m[letter][i]
            if count > maximum and count != 0:
                maximum = count
                current_letter = letter
        if not current_letter:
            sequence += undefined
        else:
            sequence += current_letter
    return Seq(sequence)


def coverage_helper(df: pl.DataFrame) -> float:
    """
    Helper function for computing coverage for the alignments of the given
    dataframe
    """

    def total_coverage(intervals: list[list]) -> float:
        return np.sum(list(map(lambda x: x[1] - x[0], intervals))) / length

    if df.is_empty():
        return np.float64(0)
    tree = it.IntervalTree()
    length = df["seq_length"][0]
    for seq, interval in zip(df["original"], df["interval"]):
        start, end = interval
        tree[start:end] = seq
    tree.merge_overlaps()
    lst = [list(v)[:-1] for v in tree]
    return total_coverage(lst)


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
    start, stop = find_char(alignment), find_char(alignment, True)
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
        pl.from_pandas(aligned_peptides).filter(pl.col("ProteinId").is_in(ids_to_keep))
    ).select(pl.col("*").exclude("UniProtKB_ID", "alignment"))
    denovo_map = peptides.join(
        id_map, left_on="sequence", right_on="original", how="inner"
    )

    metrics = (
        denovo_map.join(mismatches, on="ProteinId")
        .filter(
            (pl.col("start") <= pl.col("index")) & (pl.col("index") <= pl.col("end"))
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


def get_engine_counts(percolator_path: str, data: pd.DataFrame) -> pl.DataFrame:
    def get_engine_counts_row(id_list):
        filtered = (
            PERCOLATOR_ALL.filter(pl.col("ProteinId").is_in(id_list))
            .group_by("engine")
            .agg(num_peps=pl.col("num_peps").sum())
        )
        engine2num_peps = dict(zip(filtered["engine"], filtered["num_peps"]))
        for e in ENGINES:
            if e in engine2num_peps:
                continue
            else:
                engine2num_peps[e] = 0
        df = pl.DataFrame(engine2num_peps).with_columns(ProteinId=pl.lit(id_list[0]))
        return df.select(sorted(df.columns))

    PERCOLATOR_ALL = pl.read_csv(percolator_path, separator="\t").with_columns(
        num_peps=pl.col("peptideIds").str.count_matches(" ") + 1
    )

    ids_to_get = (
        (
            pl.from_pandas(data)
            .fill_null(value="")
            .with_columns(
                pl.col("MatchedPeptideIds").str.split(";"),
                pl.col("ProteinId").map_elements(
                    lambda x: [x], return_dtype=pl.List(pl.String)
                ),
            )
            .with_columns(
                ids_to_get=pl.col("ProteinId").list.concat(pl.col("MatchedPeptideIds"))
            )
        )
        .select("ids_to_get")
        .to_series()
    )
    ENGINES = PERCOLATOR_ALL["engine"].unique()
    rows = map(get_engine_counts_row, ids_to_get)
    result = pl.concat(rows, how="vertical", rechunk=True)
    return result
