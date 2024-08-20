import pandas as pd
import re
import sys
import polars as pl
from pathlib import Path

if "Bio_SDD" in str(Path().absolute()):
    wd = "/home/shannc/Bio_SDD/MUIC_senior_project/workflow"
else:
    wd = "/home/shannc/workflow"


outdir = f"{wd}/docs/figures/denovo_matches"
python_source = f"{wd}/bin"
sys.path.append(python_source)
import helpers as hh

COLS = ["JOIN", "Peptide", "Proteins", "n_full", "n", "engine"]

standard_search_engines = [
    "comet",
    "identipy",
    "metamorpheus",
    "msfragger",
    "msgf",
    "tide",
]

passes = ["1-First_pass", "2-Second_pass"]
prefixes = {"ND": "ND_C_indra", "default": "C_indra"}
results = f"{wd}/results"


ENGINE_FILES = dict(
    identipy="Identipy/identipy_all_pins.temp",
    comet="Comet/comet_all_pins.temp",
    msfragger="MsFragger/fragger_all_pins.temp",
    msgf="MSGF",
    metamorpheus="Metamorpheus/metamorpheus_AllPSMs_FormattedForPercolator.tab",
    tide="Tide/tide_search.target.txt",
)


def get_engine(path, engine) -> pl.DataFrame:
    engine_path = f"{path}/{ENGINE_FILES[engine]}"
    pin: pl.DataFrame
    if engine == "tide":
        pin = hh.read_tide(engine_path)
    elif engine == "msgf":
        pin_files = Path(engine_path).glob("*pin")
        pin = pl.concat([hh.read_pin(p) for p in pin_files])
    else:
        pin = hh.read_pin(engine_path)
        if engine == "identipy":
            pin = pin.with_columns(pl.col("Proteins").str.replace_all("\t", ";"))
    pin = pin.filter(pl.col("Label") == 1).with_columns(engine=pl.lit(engine))
    return pin


def format_pin(pin: pl.DataFrame) -> pl.DataFrame:
    return (
        pin.with_columns(
            pl.col("SpecId").map_elements(
                lambda x: re.sub(".*/", "", x), return_dtype=pl.String
            )
        )
        .with_columns(
            JOIN=pl.concat_str(["SpecId", "ScanNr"]),
            n_full=pl.col("Proteins").str.count_matches("T|P"),
            n=pl.col("Proteins").str.count_matches(";") + 1,
        )
        .select(COLS)
    )


outdir = f"{wd}/docs/figures/denovo_matches"
dpath = f"{results}/C_indra"
ndpath = f"{results}/ND_C_indra"
d_engines = {}
nd_engines = {}
joined = {}
for e in standard_search_engines:
    if e in {"metamorpheus", "tide"}:
        p = f"{dpath}/{passes[0]}/Engines"
        nd_p = f"{ndpath}/{passes[0]}/Engines"
    else:
        p = f"{dpath}/{passes[1]}/Engines"
        nd_p = f"{ndpath}/{passes[1]}/Engines"
    cur = format_pin(get_engine(p, e)).with_columns(
        n_denovo=pl.col("Proteins").str.count_matches("D")
    )
    nd_cur = format_pin(get_engine(nd_p, e))
    d_engines[e] = cur
    nd_engines[e] = nd_cur
    if e != "metamorpheus":
        j = (
            cur.join(nd_cur, on="JOIN", suffix="_nd")
            .with_columns(JOIN_nd=pl.lit("NA"))
            .unique("JOIN")
        )
    else:
        j = (
            cur.join(nd_cur, on="Peptide", suffix="_nd")
            .with_columns(Peptide_nd=pl.lit("NA"))
            .unique("Peptide")
        )
    joined[e] = j.select(sorted(j.columns))

all_d: pl.DataFrame = pl.concat(d_engines.values())
all_joined: pl.DataFrame = pl.concat(joined.values())

all_d.write_csv(f"{outdir}/default_prot_all.tsv", separator="\t", null_value="NA")
all_joined.write_csv(f"{outdir}/joined_prot_all.tsv", separator="\t", null_value="NA")
