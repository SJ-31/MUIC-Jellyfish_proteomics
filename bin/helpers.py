#!/usr/bin/env python

from collections import Counter
import kegg_pull.pull as kpp
from concurrent.futures import ProcessPoolExecutor
import sys
from Bio import SeqIO
import pandas as pd
from scipy.cluster.hierarchy import DisjointSet
import numpy as np
import thefuzz
import polars.selectors as cs
from rapidfuzz import distance as di
import polars as pl
import functools
import re

import thefuzz.process

AA_LOOKUP = {
    "A": "Ala",  # Alanine
    "R": "Arg",  # Arginine
    "N": "Asn",  # Asparagine
    "D": "Asp",  # Aspartic acid
    "C": "Cys",  # Cysteine
    "E": "Glu",  # Glutamic acid
    "Q": "Gln",  # Glutamine
    "G": "Gly",  # Glycine
    "H": "His",  # Histidine
    "I": "Ile",  # Isoleucine
    "L": "Leu",  # Leucine
    "K": "Lys",  # Lysine
    "M": "Met",  # Methionine
    "F": "Phe",  # Phenylalanine
    "P": "Pro",  # Proline
    "S": "Ser",  # Serine
    "T": "Thr",  # Threonine
    "W": "Trp",  # Tryptophan
    "Y": "Tyr",  # Tyrosine
    "V": "Val",  # Valine
    "n": "Nterm",
}


REGEXES = {
    "pfam": re.compile("(PF.*?)_"),
    "uniprot": re.compile("\\S*\\|\\S*\\|\\S* (.*) OS="),
    "ncbi": re.compile("\\S* (.*) \\["),
}


# Split elements in these columns and keep only the unique ones
COLS = {
    "split_keep_unique": [
        "ProteinGroupId",
        "PANTHER",
        "COG_category",
        "KEGG_Genes",
        "entry_name",
        "KEGG_ko",
        "KEGG_Pathway",
        "KEGG_Module",
        "KEGG_Reaction",
        "interpro_accession",
        "interpro_description",
        "interpro_pathways",
        "interpro_db",
        "GO",
        "GO_evidence",
        "PFAMs",
        "eggNOG_OGs",
        "BRITE",
        "Description",
        "Preferred_name",
        "MatchedPeptideIds",
    ],
    "drop_nulls_get_first": [
        "NCBI_ID",
        "UniProtKB_ID",
        "organism",
        "lineage",
        "ID_method",
        "inferred_by",
        "seed_ortholog",
    ],
    "get_first":  # These columns will already be the same,
    # (or the choice is arbitrary)
    # so its fine to just get the first
    ["ProteinId", "mass", "length"],
    "concat": [
        "header",
        "is_blast_best",
        "is_blast_one_hit",
        "q.value",
        "posterior_error_prob",
        "q_adjust",
        "pep_adjust",
        "peptideIds",
        "IdsFromDupes",
    ],
}


def resolve_duplicate_seq(data: pd.DataFrame, as_pandas=True):
    data = pl.from_pandas(data).with_columns(
        entry_name=pl.col("header").map_elements(
            entry_name_from_header, return_dtype=pl.String
        )
    )

    for c in COLS:
        COLS[c] = list(filter(lambda x: x in data.columns, COLS[c]))

    exprs = {
        "split_keep_unique": (
            pl.col(COLS["split_keep_unique"])
            .list.join(";")
            .str.split(";")
            .list.unique()
            .list.join(";")
        ),
        "drop_nulls_get_first": pl.col(COLS["drop_nulls_get_first"])
        .list.drop_nulls()
        .list.first(),
        "get_first": pl.col(COLS["get_first"]).list.first(),
        "concat": pl.col(COLS["concat"]).list.join(";"),
    }

    all_cols = [_ for col in COLS.values() for _ in col]
    joined: pl.DataFrame = (
        data.group_by("seq")
        .agg(pl.col(all_cols), pl.col("ProteinId").alias("IdsFromDupes"))
        .with_columns(*exprs.values())
        .select(data.columns)
    )
    if as_pandas:
        return joined.to_pandas()
    return joined


def entry_name_from_header(header: str) -> str:
    if not header:
        return header
    if REGEXES["uniprot"].match(header):
        return REGEXES["uniprot"].findall(header)[0]
    elif REGEXES["ncbi"].match(header):
        return REGEXES["ncbi"].findall(header)[0]
    return header


def parse_named(named_mod, count, sep="|"):
    s = named_mod.split(":")
    info = s[1].split("_")
    return f"{info[2]}{sep}{info[0]}{sep}{count}"


def get_mods(peptide, sep="|"):
    if not "[" in peptide:
        return "NA"
    mass_changes: dict = Counter(re.findall(r"(.)\[([0-9\.]+)\]", peptide))
    mods = [f"{AA_LOOKUP[k[0]]}{sep}{k[1]}{sep}{v}" for k, v in mass_changes.items()]
    if "_" in peptide:
        named_mods = Counter(re.findall(r"\[([A-Za-z:_]+)\]", peptide))
        mods.extend([parse_named(k, v, sep) for k, v in named_mods.items()])
    return ";".join(mods)


def resolve_matches(dlfq: pl.DataFrame, df: pl.DataFrame):
    """
    Map UPs (e.g. denovo, transcriptome peptides) that
    matched to full-length proteins in `df` via BLAST to the UP intensities originally
    recorded in `dlfq`.
    Allows the protein intensities to properly account for the UPs, which are effectively
    treated as ions of the protein
    """
    mp: str = "MatchedPeptideIds"
    has_matched = (
        df.filter(pl.col(mp).is_not_null())
        .with_columns(pl.col(mp).str.split(";"))
        .explode(mp)
    )
    up_matches = dlfq.filter(pl.col("protein").is_in(has_matched[mp]))
    return (
        up_matches.join(
            has_matched.select(["ProteinId", mp]),
            left_on="protein",
            right_on=mp,
        )
        .drop("protein")
        .rename({"ProteinId": "protein"})
        .select(dlfq.columns)
    )


def py_cat(lines: list[str], filename: str, append: bool = False):
    text = "\n".join([str(l) for l in lines])
    if not append:
        with open(filename, "w") as f:
            f.write(text)
    else:
        with open(filename, "a") as f:
            f.write(text)


def find_matches(queries, targets) -> pd.DataFrame:
    matches: dict = {"query": [], "best_hit": [], "similarity": []}
    for q in queries:
        find = thefuzz.process.extract(
            q, targets, scorer=di.Levenshtein.normalized_similarity
        )
        if find:
            matches["best_hit"].append(find[0][0])
            matches["similarity"].append(find[0][1])
        else:
            matches["best_hit"].append("None")
            matches["similarity"].append(0)
        matches["query"].append(q)
    return pd.DataFrame(matches)


def find_matches_par(queries, targets) -> pd.DataFrame:
    matches: dict = {"query": [], "best_hit": [], "similarity": []}

    def helper(q):
        result = thefuzz.process.extract(
            q, targets, scorer=di.Levenshtein.normalized_similarity
        )
        if result:
            return (q, result[0][0], result[0][1])
        return (q, "None", 0)

    with ProcessPoolExecutor() as exec:
        matched = exec.map(helper, queries)

    for m in matched:
        matches["query"].append(m[0])
        matches["best_hit"].append(m[1])
        matches["similarity"].append(m[2])

    return pd.DataFrame(matches)


def get_top3(dlfq_path: str, df: pd.DataFrame):
    """
    Estimate absolute protein abundance with Top3 method
    Reports protein intensity as the average of the protein's top 3 most intense peptides
    """
    df = pl.from_pandas(df)
    dlfq = (
        pl.read_csv(dlfq_path, separator="\t", infer_schema_length=None)
        .with_columns(pl.col("protein").str.split(";"))
        .explode("protein")
    )
    full_prot = dlfq.filter(pl.col("protein").str.contains("P"))
    full_prot = pl.concat([full_prot, resolve_matches(dlfq, df)])

    n_peptides_matched = full_prot.group_by("protein").len()
    at_least_3 = n_peptides_matched.filter(pl.col("len") >= 3)

    sample_files = dlfq.select(cs.numeric()).columns
    ranking_exprs = [
        pl.col(s).rank(descending=True).over("protein").alias(f"rank_{s}")
        for s in sample_files
    ]

    ranked = dlfq.filter(pl.col("protein").is_in(at_least_3["protein"])).with_columns(
        ranking_exprs
    )

    all_top3: list[pl.DataFrame] = []
    for s in sample_files:
        t3 = (
            ranked.filter(pl.col(f"rank_{s}") < 4)
            .group_by("protein")
            .agg(pl.col(s).mean().alias(f"top3-{s}"))
        )
        all_top3.append(t3)

    top3: pl.DataFrame = (
        functools.reduce(
            lambda x, y: x.join(y, on="protein", how="full", coalesce=True),
            all_top3,
        )
        .with_columns(cs.numeric().log(base=2))
        .with_columns(cs.numeric().replace(-np.inf, 0))
    )
    top3 = add_mean_median(top3, "top3")
    return top3


def write_new_dlfq(dlfq_path: str, db_file: str, output: str):
    df = pl.from_pandas(pd.read_csv(db_file, sep="\t"))
    dlfq = (
        pl.read_csv(dlfq_path, separator="\t", infer_schema_length=None)
        .with_columns(pl.col("protein").str.split(";"))
        .explode("protein")
    )
    matches = resolve_matches(dlfq, df)
    has_prot = dlfq.filter(pl.col("protein").str.contains("P"))
    new_dlfq = (
        pl.concat([matches, has_prot])
        .group_by("ion")
        .agg(pl.col("protein"), cs.numeric().first())
        .with_columns(pl.col("protein").list.unique().list.join(separator=";"))
        .select(matches.columns)
    )
    new_dlfq.write_csv(output, separator="\t")


def write_fasta(headers: list[str], seqs: list[str], filename: str) -> None:
    text = "\n".join([f">{h}\n{s}" for h, s in zip(headers, seqs)])
    with open(filename, "w") as w:
        w.write(text)


def add_mean_median(df: pl.DataFrame, prefix: str):
    numeric_only = df.select(cs.numeric())
    return df.with_columns(
        numeric_only.mean_horizontal(ignore_nulls=True).alias(f"{prefix}_mean"),
        pl.Series(np.nanmedian(numeric_only, axis=1)).alias(f"{prefix}_median"),
    )


def prefix_numeric_cols(df: pl.DataFrame, prefix: str, with_mean_median: bool = True):
    numeric_cols = df.select(cs.numeric()).columns
    name_mapping = {c: f"{prefix}-{c}" for c in numeric_cols}
    if with_mean_median:
        df = add_mean_median(df, prefix)
    return df.rename(name_mapping)


def read_dlfq_prot(dlfq_path: str) -> pl.DataFrame:
    result = pl.read_csv(dlfq_path, separator="\t")
    named_cols = result.columns
    named_cols.remove("")
    result = result.select(named_cols)
    result = prefix_numeric_cols(result, "directlfq")
    exploded = (
        (
            result.with_columns(pl.col("protein").str.split(";"))
            .explode("protein")
            .group_by("protein")
            .agg(cs.numeric().mean())
        )
        .with_columns(cs.numeric().log(base=2))
        .with_columns(cs.numeric().replace(-np.inf, 0))
    )
    return exploded


def flatten_by(lst, by=";"):
    return [i for string in lst for i in string.split(by)]


def fasta2df(fasta: str) -> pl.DataFrame:
    tmp = {"header": [], "seq": []}
    for entry in SeqIO.parse(fasta, format="fasta"):
        tmp["header"].append(entry.id)
        tmp["seq"].append(entry.seq)
    return pl.DataFrame(tmp)


def get_queries(query_path: str, header_mapping: str) -> tuple[pl.DataFrame, dict]:
    queries = fasta2df(query_path)
    header_map = pl.read_csv(header_mapping, separator="\t")
    queries = (
        queries.join(header_map, left_on="header", right_on="id")
        .rename({"header_right": "H"})
        .drop("mass", "seq_right", "length")
    )
    header2id: dict = dict(zip(queries["H"], queries["header"]))
    return queries, header2id


def retrieve_saved_eggnog(
    query_path: str,
    saved_path: str,
    header_mapping: str,
) -> tuple[pl.DataFrame, pl.DataFrame]:
    queries, header2id = get_queries(query_path, header_mapping)
    filenames = ["annotations", "orthologs", "seed_orthologs"]
    qcols = ["#query", "#query", "#qseqid"]
    dfs = [
        pl.read_csv(
            f"{saved_path}/eggnog_{t}.tsv",
            separator="\t",
            comment_prefix="##",
            null_values="NA",
        )
        for t in filenames
    ]
    # Retrieve annotations that were already found for "queries"
    filtered = [
        df.filter(pl.col("header").is_in(queries["H"]))
        .with_columns(
            pl.col("header")
            .map_elements(lambda x: header2id[x], return_dtype=pl.String)
            .alias(qcol)
        )
        .select([qcol] + df.columns)
        .drop("header")
        for df, qcol in zip(dfs, qcols)
    ]
    results = {m: df for m, df in zip(filenames, filtered)}
    queries = queries.filter(~pl.col("H").is_in(dfs[0]["header"]))  # Remove all
    # queries that were previously found
    return results, queries


def retrieve_saved_interpro(
    query_path: str,
    saved_path: str,
    header_mapping: str | None = None,
) -> tuple[pl.DataFrame, pl.DataFrame]:
    queries, header2id = get_queries(query_path, header_mapping)
    interpro = pl.read_csv(saved_path, separator="\t")
    filtered = (
        interpro.filter(pl.col("header").is_in(queries["H"]))
        .with_columns(
            pl.col("header")
            .map_elements(lambda x: header2id[x], return_dtype=pl.String)
            .alias("query")
        )
        .select(["query"] + interpro.columns)
        .drop("header")
    )
    queries = queries.filter(~pl.col("H").is_in(interpro["header"]))
    return filtered, queries


def get_denovo_stats(denovo_file: str, identifications: str, header_map_path: str):
    data: pl.DataFrame = pl.read_csv(identifications, separator="\t", null_values="NA")
    header_map = pl.read_csv(header_map_path, separator="\t")
    matched: list = flatten_by(
        data.filter(pl.col("MatchedPeptideIds").is_not_null())["MatchedPeptideIds"]
    ) + list(data.filter(pl.col("ProteinId").str.contains("D"))["ProteinId"])
    matched = list(set(filter(lambda x: "D" in x, matched)))
    denovo: pl.DataFrame = fasta2df(denovo_file).join(header_map, on="header")
    print(f"Number of denovo peptides given: {denovo.shape[0]}")
    print(f"Number of denovo peptides matched: {len(matched)}")
    was_matched = denovo.filter(pl.col("id").is_in(matched))
    percent_id = was_matched.shape[0] / denovo.shape[0]
    print(f"Percent identified: {percent_id * 100}")


def merge_subsets(items: dict) -> DisjointSet:
    """
    Given a bunch of sets that may or may not be subsets of one another,
    resolve the sets such that no sets left are subsets of one another (i.e. finding the largest sets that contain the others)

    If a set A is a subset of B AND C, and B and C are disjoint, whichever comes first
    will take on A
    """
    DS: DisjointSet = DisjointSet(items.keys())
    candidates = set(items.keys())
    while candidates:
        were_subsets: set = set()
        for x in candidates:
            for y in candidates:
                if x == y:
                    continue
                x_vals, y_vals = items.get(x), items.get(y)
                if x_vals <= y_vals and x not in were_subsets:
                    DS.merge(x, y)
                    were_subsets.add(x)
                elif y_vals <= x_vals and y not in were_subsets:
                    DS.merge(x, y)
                    were_subsets.add(y)
        if not were_subsets:
            break
        candidates -= were_subsets
    return DS


def subsets2df(DS: DisjointSet, items: dict) -> pl.DataFrame:
    """The representative in this case is the set that contains
    all others in the subset
    """
    size_key: dict = {k: len(v) for k, v in items.items()}
    results: dict = {"item": [], "set": [], "representative": []}
    for i, subset in enumerate(DS.subsets()):
        ordering: dict = {}
        for e in subset:
            ordering[e] = size_key.get(e)
            results["item"].append(e)
            results["set"].append(i)
        ordered: list = sorted(ordering.items(), key=lambda x: x[1], reverse=True)
        largest = ordered[0][0]
        results["representative"].extend([largest] * len(subset))
    return pl.DataFrame(results)


def clean_peptide(peptide):
    if re.search("[a-z]", peptide):
        mod_regex = re.compile(r"\[[A-Za-z_]+\:[_A-Za-z]+\]")
        peptide = re.subn(mod_regex, "", peptide)[0]
    peptide = peptide.replace("X", "")
    peptide = re.sub("^n", "", peptide)
    return "".join(re.findall("[A-Z]+", peptide))


def clean_peptide_joined(peptide_str) -> str:
    cleaned = ""

    if peptide_str:
        splits = set(peptide_str.split(";"))
        cleaned = ";".join([clean_peptide(p) for p in splits])
    return cleaned


def get_unique_peptides_py(data: pd.DataFrame, filename: str) -> pd.DataFrame:
    df = pl.from_pandas(data)
    key = (
        df.unique("ProteinId")
        .select("ProteinId", "peptideIds")
        .with_columns(
            unique_peptides=pl.col("peptideIds").map_elements(
                clean_peptide_joined, return_dtype=pl.String
            )
        )
    )
    df = df.join(key.select("ProteinId", "unique_peptides"), on="unique_peptides")
    df.write_csv(filename, separator="\t", null_value="NA")
    return df.to_pandas()


def group_by_subsets(data: pd.DataFrame) -> pd.DataFrame:
    df: pl.DataFrame = pl.from_pandas(data).with_columns(
        pl.col("unique_peptides").str.split(";").alias("split_peps_temp")
    )
    pep_dict = {k: set(v) for k, v in zip(df["ProteinId"], df["split_peps_temp"])}
    merged = merge_subsets(pep_dict)
    result = subsets2df(merged, pep_dict).rename(
        {"set": "GroupSB", "representative": "sb_rep"}
    )
    return (
        df.join(result, left_on="ProteinId", right_on="item")
        .drop("split_peps_temp")
        .to_pandas()
    )


def get_unique_peptides_py(filename: str) -> None:
    df = pl.read_csv(filename, separator="\t", null_values="NA").drop(
        cs.contains("unique_peptides")
    )
    key = (
        df.unique("ProteinId")
        .select("ProteinId", "peptideIds")
        .with_columns(
            unique_peptides=pl.col("peptideIds").map_elements(
                clean_peptide_joined, return_dtype=pl.String
            )
        )
    )
    df = df.join(key.select("ProteinId", "unique_peptides"), on="ProteinId")
    df.write_csv(filename, separator="\t", null_value="NA")


def is_kegg_header(query: str) -> bool:
    return query in {
        "NAME",
        "CLASS",
        "DESCRIPTION",
        "COMPOUND",
        "REFERENCE",
        "AUTHORS",
        "TITLE",
        "JOURNAL",
        "DOI",
        "ENTRY",
        "CLASS",
        "PATHWAY_MAP",
        "DBLINKS",
        "ORTHOLOGY",
        "///",
    }


def get_splits(line: str):
    return list(filter(lambda x: x != "", line.split(" ")))


class KeggParser:

    def __init__(self, kegg_text: str, fill_dict: dict = None) -> None:
        self.lines = kegg_text.splitlines()
        self.index = 0
        self.length = len(self.lines)
        if fill_dict:
            self.fill_into = True
            self.parsed = fill_dict
        else:
            self.fill_into = False
            self.parsed: dict = {}

    def get_entry_list(self, start) -> list:
        lst = []
        lst.append(start[1])
        l = self.lines[self.index + 1]
        split_tmp = get_splits(l)
        while not is_kegg_header(split_tmp[0]) and self.index < self.length:
            lst.append(split_tmp[0])
            self.index += 1
            if self.index >= self.length:
                return lst
            l = self.lines[self.index]
            split_tmp = get_splits(l)
        if len(lst) > 1:
            self.index -= 1
        return lst

    def fill(self, key, value) -> None:
        if self.fill_into:
            self.parsed[key].append(value)
        else:
            self.parsed[key] = value

    def __call__(self) -> dict:
        to_join_lines: dict = {
            "NAME": "name",
            "CLASS": "class",
            "DESCRIPTION": "description",
        }
        multiline_entries: dict = {
            "ORTHOLOGY": "orthologs",
            "COMPOUND": "compounds",
            "REL_PATHWAY": "related",
        }

        not_filled: set = {
            "name",
            "class",
            "description",
            "orthologs",
            "compounds",
            "related",
        }
        while self.index < self.length:
            l: str = self.lines[self.index]
            split: list = get_splits(l)
            first: str = split[0]
            if first == "ENTRY":
                self.fill("entry", split[1])
                self.fill("db", split[2])
            elif first in to_join_lines:
                self.fill(to_join_lines[first], " ".join(split[1:]))
                not_filled.remove(to_join_lines[first])
            elif first in multiline_entries:
                joined = ";".join(self.get_entry_list(split))
                self.fill(multiline_entries[first], joined)
                not_filled.remove(multiline_entries[first])
            self.index += 1
        if not_filled:
            for f in not_filled:
                self.fill(f, None)
        return self.parsed


def get_kegg_metadata(kegg_ids: list) -> pl.DataFrame:
    single_pull = kpp.SinglePull()
    data = {
        "entry": [],
        "db": [],
        "name": [],
        "class": [],
        "description": [],
        "orthologs": [],
        "compounds": [],
        "related": [],
    }
    for k in kegg_ids:
        find = single_pull.pull_dict([k])
        if k in find[0].successful_entry_ids:
            KeggParser(find[1][k], fill_dict=data)()
    return pl.DataFrame(data)


def parse_args():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-t", "--task")
    parser.add_argument("-i", "--input")
    parser.add_argument("--save_type")
    parser.add_argument("-m", "--maxlfq")
    parser.add_argument("-p", "--top3")
    parser.add_argument("-s", "--saved")
    parser.add_argument("--seq_header_mapping")
    parser.add_argument("-d", "--dlfq_input")
    parser.add_argument("-o", "--output")
    args = vars(parser.parse_args())
    return args


if __name__ == "__main__" and len(sys.argv) > 1 and not "radian" in sys.argv[0]:
    args = parse_args()
    if args["task"] == "write_dlfq":
        write_new_dlfq(args["dlfq_input"], args["input"], args["output"])
    elif args["task"] == "top3":
        df = pd.read_csv(args["input"], sep="\t")
        result = get_top3(args["dlfq_input"], df)
        result.write_csv(args["output"], separator="\t")
    elif args["task"] == "merge":
        dlfq = read_dlfq_prot(args["dlfq_input"])
        maxlfq = pl.read_csv(args["maxlfq"], separator="\t", null_values="NA").rename(
            {"ProteinId": "protein"}
        )
        maxlfq = add_mean_median(maxlfq, "maxlfq")
        top3 = pl.read_csv(args["top3"], separator="\t")
        merged: pl.DataFrame = functools.reduce(
            lambda x, y: x.join(y, on="protein", how="full", coalesce=True),
            [dlfq, maxlfq, top3],
        ).rename({"protein": "ProteinId"})
        sm = pl.read_csv(args["seq_header_mapping"], separator="\t").select(
            "id", "header"
        )
        merged = merged.join(sm, left_on="ProteinId", right_on="id")
        merged.write_csv(args["output"], separator="\t", null_value="NA")
    elif args["task"] == "get_saved" and args["save_type"] == "eggnog":
        eggnog, queries = retrieve_saved_eggnog(
            query_path=args["input"],
            saved_path=args["saved"],
            header_mapping=args["seq_header_mapping"],
        )
        write_fasta(queries["header"], queries["seq"], filename="new_query.fasta")
        for t, df in eggnog.items():
            df.write_csv(f"saved_{t}.tsv", separator="\t", include_header=False)
    elif args["task"] == "get_saved" and args["save_type"] == "interpro":
        interpro, queries = retrieve_saved_interpro(
            query_path=args["input"],
            saved_path=args["saved"],
            header_mapping=args["seq_header_mapping"],
        )
        interpro.write_csv("saved_interpro.tsv", separator="\t", include_header=False)
        write_fasta(queries["header"], queries["seq"], filename="new_query.fasta")
    elif args["task"] == "denovo_stats":
        get_denovo_stats(args["input"], args["saved"], args["seq_header_mapping"])
