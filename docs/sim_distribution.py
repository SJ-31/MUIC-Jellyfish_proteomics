import sys
import itertools
import matplotlib.pyplot as plt

sys.path.append("/home/shannc/Bio_SDD/MUIC_senior_project/workflow/bin")
import null_distribution as nd
import polars as pl
import h5py
import concurrent.futures as futures
import numpy as np


results_file = "/home/shannc/Bio_SDD/MUIC_senior_project/workflow/results/C_indra/Analysis/C_indra_all_wcategory.tsv"
data = pl.read_csv(results_file, separator="\t", null_values="NA")
sem_file = "/home/shannc/Bio_SDD/MUIC_senior_project/workflow/results/C_indra/Analysis/sem_matrices.hdf5"
prottrans_file = "/home/shannc/Bio_SDD/MUIC_senior_project/workflow/results/C_indra/Analysis/Embeddings_prottrans/distances.hdf5"
outdir = "/home/shannc/Bio_SDD/MUIC_senior_project/workflow/.cache/"
file = h5py.File(sem_file)
pt_file = h5py.File(prottrans_file)

mf = nd.Dist(file, ontology="MF")
bp = nd.Dist(file, ontology="BP")
cc = nd.Dist(file, ontology="CC")
prottrans = nd.Dist(pt_file, metric="cosine")

try:
    sem_distribution = np.loadtxt(f"{outdir}/sem_distribution")
except FileNotFoundError:
    value = 1_000_000
    with futures.ProcessPoolExecutor() as exec:
        result: list[futures.Future] = list(
            exec.map(nd.null_distribution, [cc, mf, bp], [value, value, value])
        )
    sem_distribution = np.mean(result, axis=0)
    sem_distribution.tofile(f"{outdir}/sem_distribution", sep=" ")

try:
    pt_distribution = np.loadtxt(f"{outdir}/pt_distribution")
except FileNotFoundError:
    pt_distribution = nd.null_distribution(prottrans, value)
    pt_distribution.tofile(f"{outdir}/pt_distribution", sep=" ")

sem_pdf = nd.PdfEstimate(sem_distribution, 7, (0, 1))
pt_pdf = nd.PdfEstimate(pt_distribution, 6, (0, 1))

prot = data.select("ProteinId", "Group")
grouped = (
    prot.group_by("Group")
    .agg(pl.col("ProteinId"), pl.len())
    .sort("len", descending=True)
)
# %%

sem_pdf.plot(with_hist=False)

grouped = grouped.with_columns(
    average_dist=pl.col("ProteinId").map_elements(
        lambda x: nd.average_dist(x, prottrans), return_dtype=pl.Float64
    )
)
grouped = grouped.with_columns(
    p_value=pl.col("average_dist").map_elements(
        lambda x: pt_pdf.get_pr(x), return_dtype=pl.Float64
    )
)
# grouped["p_value"]
# grouped.drop("ProteinId").write_csv("validated.tsv", separator="\t")

grouped

blah = pl.Series(list(range(grouped.shape[0])))


def add_dynamic(series, name, df):
    return df.with_columns(series.alias(name))


add_dynamic(blah, "ba", grouped)
