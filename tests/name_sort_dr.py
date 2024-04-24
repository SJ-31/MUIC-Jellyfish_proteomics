#!/usr/bin/env python
import pathlib
import shutil
from subprocess import Popen
from pathlib import Path
from natsort import natsorted
import pandas as pd
import re

DIST = ["euclidean", "cosine"]
TECHNIQUES = ["tsne", "pcoa", "umap"]


def rreplace(string: str, old: str, new: str, count: int = 1) -> str:
    """
    Replace from the right of a string
    """
    rev = string[::-1]
    replaced = rev.replace(old[::-1], new[::-1], count)[::-1]
    return replaced


def fileSuffix(filename, suffix):
    """
    Apply `suffix` to filename while being aware of any file extensions
    """
    ext = filename[filename.rfind(".") + 1 :]
    if not ext:
        raise ValueError("Filename has no extension!")
    output = f'{rreplace(filename, f".{ext}", "")}{suffix}.{ext}'
    return output


def deletePaths(paths: list[str]) -> None:
    """
    :param paths: list of pathnames to delete
    """
    for f in paths:
        current = Path(f)
        if current.exists() and current.is_dir():
            shutil.rmtree(current.absolute())
        elif current.exists():
            current.unlink()


def addLabel(
    image: str,
    label: str,
    coord: tuple,
    gravity: str = "None",
    color: str = "black",
    output: str = "",
    size: int = 50,
) -> str:
    if not output:
        output = fileSuffix(image, "-LABELED")
    command = [
        "magick",
        image,
        "-gravity",
        gravity,
        "-pointsize",
        str(size),
        "-fill",
        color,
        "-annotate",
        f"+{coord[0]}+{coord[1]}",
        label,
        output,
    ]
    with Popen(command) as m:
        m.communicate()
    return output


def arrangeImages(
    images: list, output: str, mode: str, remove: bool = False
) -> str:
    com = ["convert"]
    if mode == "horizontal":
        com.append("+append")
    elif mode == "vertical":
        com.append("-append")
    com.extend(images)
    com.append(output)
    with Popen(com) as m:
        m.communicate()
        if m.returncode != 0:
            print(com)
    if remove:
        deletePaths(images)
    return output


# Arrange pca and pcoa
def getPcaPcoa(input_path, output_path):
    pca = []
    pcoa = []
    for tech, lst in zip(["pca", "pcoa"], [pca, pcoa]):
        for d in input_path.glob(f"{tech}*"):
            if d.is_dir():
                model = re.sub(".*_", "", d.name)
                img = [png for png in d.rglob("*.png")][0]
                lbl = f"Model: {model}"
                labeled = addLabel(
                    str(img), lbl, (300, 50), gravity="NorthEast", size=70
                )
                lst.append(labeled)
    pca = sorted(pca)
    pcoa = sorted(pcoa)
    p1 = arrangeImages(pca, f"{input_path}/pca_all.png", "horizontal")
    p2 = arrangeImages(pcoa, f"{input_path}/pcoa_all.png", "horizontal")
    arrangeImages([p1, p2], f"{output_path}/pca_pcoa.png", "vertical", True)


def cropSide(image: str, by: int, side: str, output: str = "") -> str:
    if not output:
        output = fileSuffix(image, "-CROPPED")
    match side:
        case "t":
            crop_str = f"+0+{by}"
        case "b":
            crop_str = f"+0-{by}"
        case "l":
            crop_str = f"+{by}+0"
        case "r":
            crop_str = f"-{by}+0"
        case other:
            raise ValueError(
                "Side must be one of [t]op, [b]ottom, [l]eft or [r]ight"
            )
    command = ["magick", image, "-crop", crop_str, "+repage", output]
    with Popen(command) as m:
        m.communicate()
    return output


# For arranging tsne & umap
def arrangeToModelMetric(
    technique: str, model_path: Path, dmetric: str, outdir: str = ""
) -> None:
    model = re.sub(".*_", "", str(model_path))
    pics = [str(png) for png in model_path.rglob(f"*{dmetric}*.png")]
    pics = natsorted(pics)
    formatted = []
    for index, img in enumerate(pics):
        if index == 0:
            labeled = addLabel(
                pics[0],
                f"Model: {model} ",
                (200, 30),
                size=70,
                gravity="NorthEast",
            )
            formatted.append(labeled)
        else:
            formatted.append(cropSide(img, 600, "r"))
    ncols = 5
    row_pics = []
    count = 0
    while count <= len(pics):
        sliced = formatted[count : count + ncols]
        count += ncols
        if len(sliced) <= 1:
            continue
        row_pics.append(
            arrangeImages(
                sliced,
                f"{model_path}/r{count}.png",
                "horizontal",
            )
        )
    if not outdir:
        outdir = model_path
    arrangeImages(
        row_pics,
        f"{outdir}/{technique}_{model}_{dmetric}_all.png",
        "vertical",
        True,
    )
    deletePaths(formatted)


def arrangeMultiple(technique, input_path, output_path) -> None:
    paths = list(input_path.glob(f"{technique}*"))
    for t in paths:
        for m in DIST:
            arrangeToModelMetric(technique, t, m, outdir=output_path)


def getTrustworthinessModel(technique, model_path):
    model = re.sub(".*_", "", str(model_path))
    files = natsorted(list(model_path.rglob("*trustworthiness.txt")))
    tw_dict = {
        "trustworthiness": [],
        "technique": [],
        "model": [],
        "distance_metric": [],
        "extra_params": [],
    }
    for f in files:
        if technique in {"tsne", "umap"}:
            extra = f.parent.stem.split("_")
            tw_dict["extra_params"].append(extra[0])
            tw_dict["distance_metric"].append(extra[1])
        else:
            tw_dict["extra_params"].append("-")
            if technique == "pcoa":
                tw_dict["distance_metric"].append("cosine")
            else:
                tw_dict["distance_metric"].append("euclidean")
        with open(f, "r") as r:
            trust = r.read()
            if trust:
                tw_dict["trustworthiness"].append(trust)
            else:
                tw_dict["trustworthiness"].append(pd.NA)
        tw_dict["technique"].append(technique)
        tw_dict["model"].append(model)
    df = pd.DataFrame(tw_dict)
    df = df.sort_values("distance_metric", kind="mergesort")
    return df


def getTrustWorthiness(input_path, output_path) -> None:
    t_df = pd.DataFrame()
    for t in ["pca", "pcoa", "umap", "tsne"]:
        paths = list(input_path.glob(f"{t}*"))
        for p in paths:
            t_df = pd.concat([t_df, getTrustworthinessModel(t, p)])
    t_df.to_csv(f"{output_path}/trustworthiness.tsv", sep="\t", index=False)


def listFiles(directory: str | pathlib.PosixPath, pattern: str):
    if isinstance(directory, pathlib.PosixPath):
        return [str(f) for f in Path(directory).rglob(pattern)]
    return [str(f) for f in Path(directory).rglob(pattern)]


if "Bio_SDD" in str(Path(".").absolute()):
    wd = "/home/shannc/Bio_SDD/MUIC_senior_project/workflow/tests"
else:
    wd = "/home/shannc/workflow/tests"


COMPARE_PROTEIN_EMBEDDINGS = False
if COMPARE_PROTEIN_EMBEDDINGS:
    # Arrange protein embedding DR
    outdir = Path(f"{wd}/testthat/output/protein_dr")
    all_arranged = outdir.joinpath("all_compared")
    getTrustWorthiness(outdir, all_arranged)
    getPcaPcoa(outdir, all_arranged)
    arrangeMultiple("tsne", outdir, all_arranged)
    arrangeMultiple("umap", outdir, all_arranged)


def compareSsA2v(outdir):
    arranged = []
    for v in ["semantic", "a2v"]:
        figures = listFiles(outdir.joinpath(v), "*png")
        joined = arrangeImages(
            figures, f"{outdir}/{v}_joined.png", "horizontal"
        )
        arranged.append(joined)
    arrangeImages(arranged, f"{outdir}/comparison.png", "vertical", True)


def whichTechnique(file_path):
    for t in TECHNIQUES:
        if t in file_path:
            return t


def ssA2vGetTrustWorthiness(outdir):
    comparison_dict = {"type": [], "trustworthiness": [], "technique": []}
    for v in ["semantic", "a2v"]:
        ss = outdir.joinpath(v)
        files = listFiles(ss, "*trustworthiness*")
        for file in files:
            comparison_dict["type"].append(v)
            comparison_dict["technique"].append(whichTechnique(file))
            with open(file, "r") as f:
                trst = f.read()
            comparison_dict["trustworthiness"].append(trst)
    pd.DataFrame(comparison_dict).to_csv(
        f"{outdir}/trustworthiness.tsv", sep="\t", index=False
    )


COMPARE_SS_A2V = True
if COMPARE_SS_A2V:
    # Arrange images for ss and a2v comparison
    # outdir = Path(f"{wd}/testthat/output/go_dr_protein_level")
    outdir = Path(f"{wd}/testthat/output/go_dr_go_level")
    compareSsA2v(outdir)
    ssA2vGetTrustWorthiness(outdir)
