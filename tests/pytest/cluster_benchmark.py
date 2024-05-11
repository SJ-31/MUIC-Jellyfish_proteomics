import os
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns
import umap

# from context import get_distances as gd
import h5py


if "Bio_SDD" in __file__:
    wd = "/home/shannc/Bio_SDD/MUIC_senior_project/workflow"
else:
    wd = "/home/shannc/workflow"
os.chdir(f"{wd}/tests/pytest")


OUTDIR = "./output/umap_tests"
if not os.path.isdir(OUTDIR):
    os.makedirs(OUTDIR)

embeddings = h5py.File(
    f"{wd}/tests/nf-test-out/C_indra_prottrans_embeddings/distances.hdf5"
)
combined_results = pd.read_csv(
    f"{wd}/results/C_indra_A/1-First_pass/C_indra_all_wcoverage.tsv", sep="\t"
)
color = "category"
dist = pd.DataFrame(embeddings["metric/cosine"][:])
dist.insert(0, "ProteinId", embeddings["names"][:].astype("str"))
dist = pd.merge(combined_results[["ProteinId", color]], dist, on="ProteinId")
min_distances = np.arange(0.001, 0.05, 0.005)
n_neighbors = np.arange(5, 50, 10)


fig, ax = plt.subplots(nrows=len(min_distances), ncols=len(n_neighbors))
for index, md in enumerate(min_distances):
    for index2, nn in enumerate(n_neighbors):
        reducer = umap.UMAP(
            min_dist=md, n_components=2, n_neighbors=nn, metric="precomputed"
        )
        current_ax = ax[index, index2]
        result = reducer.fit_transform(dist.iloc[:, 2:])
        result = pd.concat([pd.DataFrame(result), dist[color]], axis=1)
        sns.scatterplot(
            x=result[0],
            y=result[1],
            hue=result["category"],
            ax=current_ax,
            legend=False,
        )
        current_ax.set_title(f"md: {md}, nn: {nn}")
fig.set_size_inches(50, 60)
fig.savefig(f"{OUTDIR}/umap_compare.png")
