#!/usr/bin/env python
import matplotlib.pyplot as plt
import numpy as np
import numpy.polynomial.polynomial as poly
import scipy.integrate as integrate
import itertools
import polars as pl
import h5py
import concurrent.futures as futures
import scipy.optimize as op
import h5py

ITER = 1_000_000


class PdfEstimate:
    def __init__(
        self, data: np.array, degree: int, domain=tuple[float, float], bins="auto"
    ) -> None:
        self.data = data
        ys, xs = np.histogram(data, bins=bins, density=True)
        self.pdf = poly.Polynomial = poly.Polynomial.fit(
            xs[:-1], ys, deg=degree, domain=domain
        )
        self.domain = domain

    def critical_value(self, alpha: float, left_tailed=True):
        start, stop = self.domain
        if left_tailed:
            f = lambda x: integrate.quad(self.pdf, start, x)[0] - alpha
        else:
            f = lambda x: integrate.quad(self.pdf, x, stop)[0] - alpha
        guess = np.linspace(start, stop)[11]
        solve = op.root(f, x0=[guess])
        if solve.success:
            return solve.x[0]
        return None

    def get_pr(self, x: float, left_tailed=True) -> float:
        if left_tailed:
            return integrate.quad(self.pdf, self.domain[0], x)[0]
        return integrate.quad(self.pdf, x, self.domain[1])[0]

    def plot(self, critical_value=None, with_hist=True) -> None:
        start, stop = self.domain
        xs, ys = plot_in_range(self.pdf, start, stop, num=1000)
        if with_hist:
            plt.hist(self.data, bins="auto")
        if critical_value:
            cv = self.critical_value(critical_value)
            plt.vlines(cv, ymin=start, ymax=self.pdf(cv), colors="green")
            plt.fill_between(xs, ys, where=(xs < cv), color="green")
            plt.axis([start, stop, 0, max(ys) + np.mean(ys) / 2])


class Dist:
    def __init__(
        self, file: h5py.File, ontology: str = None, metric: str = None
    ) -> None:
        if not ontology and not metric:
            raise ValueError("Type of distance matrix not specified")
        if ontology:
            self.mat: np.array = file[f"matrices/{ontology}"][:]
            names = file[f"names/{ontology}"][:]
        elif metric:
            self.mat: np.array = file[f"metric/{metric}"][:]
            names = file["names"][:]
        if not (self.mat.T == self.mat).all():
            raise ValueError("Not a symmetric distance matrix")
        self.name_map: dict[str, int] = {v.decode(): i for i, v in enumerate(names)}

    def distance_name(self, ab: tuple[str, str]):
        a, b = ab
        a_index, b_index = self.name_map.get(a), self.name_map.get(b)
        if not a_index or not b_index:
            return np.nan
        return self.mat[a_index, b_index]

    def distance_index(self, ab: tuple[int, int]):
        a, b = ab
        shape = self.mat.shape
        if a < shape[0] and b < shape[1]:
            return self.mat[a, b]
        raise IndexError("a or b are out of the matrix")


MF: Dist
CC: Dist
BP: Dist
EMBEDDINGS: Dist


def plot_in_range(f, left_bound, right_bound, num=None, ax=None, **kwargs):
    if num:
        xs = np.linspace(start=left_bound, stop=right_bound, num=num)
    else:
        xs = np.linspace(start=left_bound, stop=right_bound)
    ys = list(map(f, xs))
    if not ax:
        plt.plot(xs, ys, **kwargs)
    else:
        ax.pl
    return xs, ys


def null_distribution(D: Dist, n: int) -> np.array:
    shape = D.mat.shape

    def new_pair():
        pair = np.random.randint(low=0, high=shape[0], size=2)
        while pair[0] == pair[1]:
            pair = np.random.randint(low=0, high=shape[0], size=2)
        return pair

    return np.array([D.distance_index(new_pair()) for _ in range(n)])


def average_dist(protein_ids, dist: Dist):
    if len(protein_ids) < 2:
        return np.float64(0)
    distances = list(
        filter(
            lambda x: not np.isnan(x),
            [dist.distance_name(p) for p in itertools.pairwise(protein_ids)],
        )
    )
    return np.mean(distances)


def get_sem_dist(cache: str) -> np.array:
    try:
        sem_distribution = np.loadtxt(f"{cache}/sem_distribution")
    except FileNotFoundError:
        with futures.ProcessPoolExecutor() as exec:
            result: list[futures.Future] = list(
                exec.map(null_distribution, [CC, MF, BP], [ITER, ITER, ITER])
            )
        sem_distribution = np.mean(result, axis=0)
        sem_distribution.tofile(f"{cache}/sem_distribution", sep=" ")
    return sem_distribution


def get_embedding_dist(cache: str) -> np.array:
    try:
        e_distribution = np.loadtxt(f"{cache}/pt_distribution")
    except FileNotFoundError:
        e_distribution = null_distribution(EMBEDDINGS, ITER)
        e_distribution.tofile(f"{cache}/pt_distribution", sep=" ")
    return e_distribution


def parse_args():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input")
    parser.add_argument("-c", "--cache")
    parser.add_argument("-p", "--protein_embedding_distances")
    parser.add_argument("-s", "--semantic_similarity")
    parser.add_argument("-g", "--grouping_col")
    parser.add_argument("-o", "--output")
    args = vars(parser.parse_args())
    return args


def average_dist_sem(protein_ids: list):
    if len(protein_ids) < 2:
        return np.float64(0)
    nums = list(
        filter(
            lambda x: not np.isnan(x),
            np.array(
                [average_dist(protein_ids, ont) for ont in [CC, BP, MF]]
            ).flatten(),
        )
    )
    return np.mean(nums)


def add_p_values(grouped_df: pl.DataFrame, mode: str, pdf: PdfEstimate) -> pl.DataFrame:
    if mode == "embedding":
        dist_fun = lambda x: average_dist(x, EMBEDDINGS)
        dist_col = "avg_embedding_distance"
        pr_fun = lambda x: pdf.get_pr(x)
    elif mode == "semantic":
        dist_fun = average_dist_sem
        dist_col = "avg_semantic_distance"
        pr_fun = lambda x: pdf.get_pr(x, left_tailed=False)
    else:
        raise ValueError("mode must be one of 'embedding' or 'semantic'")
    grouped = grouped_df.with_columns(
        pl.col("ProteinId")
        .map_elements(dist_fun, return_dtype=pl.Float64)
        .alias(dist_col)
    )
    grouped = grouped.with_columns(
        pl.col(dist_col)
        .map_elements(pr_fun, return_dtype=pl.Float64)
        .alias(f"p_value_{mode}")
    )
    return grouped


if __name__ == "__main__":
    args = parse_args()
    data = pl.read_csv(args["input"], separator="\t", null_values="NA")
    sem_file = h5py.File(args["semantic_similarity"])
    embedding_file = h5py.File(args["protein_embedding_distances"])
    MF, BP, CC = (
        Dist(sem_file, ontology="MF"),
        Dist(sem_file, ontology="BP"),
        Dist(sem_file, ontology="CC"),
    )
    EMBEDDINGS = Dist(embedding_file, metric="cosine")

    sem_distribution = get_sem_dist(args["cache"])
    embedding_distribution = get_embedding_dist(args["cache"])
    sem_pdf = PdfEstimate(sem_distribution, 7, (0, 1))
    embedding_pdf = PdfEstimate(embedding_distribution, 6, (0, 1))

    prot = data.select("ProteinId", args["grouping_col"])
    grouped = (
        prot.group_by(args["grouping_col"])
        .agg(pl.col("ProteinId"), pl.len())
        .sort("len", descending=True)
    )
    grouped = add_p_values(grouped, "semantic", sem_pdf)
    grouped = add_p_values(grouped, "embedding", embedding_pdf)
    grouped.drop("ProteinId").write_csv(args["output"], separator="\t")
