from collections import Counter
import pandas as pd
import numpy as np
import polars.selectors as cs
import polars as pl
import functools
import re

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
            has_matched.select(["ProteinId", mp]), left_on="protein", right_on=mp
        )
        .drop("protein")
        .rename({"ProteinId": "protein"})
        .select(dlfq.columns)
    )


def get_top3(dlfq_path: str, df: pd.DataFrame):
    """
    Estimate absolute protein abundance with Top3 method
    Reports protein intensity as the average of the protein's top 3 most intense peptides
    """
    df = pl.from_pandas(df)
    dlfq = (
        pl.read_csv(dlfq_path, separator="\t")
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
            lambda x, y: x.join(y, on="protein", how="full", coalesce=True), all_top3
        )
        .with_columns(cs.numeric().log())
        .with_columns(cs.numeric().replace(-np.inf, 0))
    )
    cols_only = top3.select(cs.numeric())
    top3 = top3.with_columns(
        top3_mean=cols_only.mean_horizontal(ignore_nulls=True),
        top3_median=np.nanmedian(cols_only, axis=1),
    )
    return df.join(top3, how="left", left_on="ProteinId", right_on="protein")
