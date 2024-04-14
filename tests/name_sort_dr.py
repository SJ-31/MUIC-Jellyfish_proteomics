#!/usr/bin/env python

from subprocess import Popen
from pathlib import Path
import re


def addLabel(
    image: str,
    label: str,
    position: tuple,
    color: str = "black",
    output: str = "",
    size: int = 50,
) -> str:
    if not output:
        ext = re.findall("(png|jpeg|jpg)", image)
        if not ext:
            raise ValueError("Image file has no extension!")
        ext = ext[0]
        output = f"{image.replace(f".{ext}", "")}-labeled.{ext}"
    command = [
        "convert",
        image,
        "-pointsize",
        str(size),
        "-fill",
        color,
        "-annotate",
        f"+{position[0]}+{position[1]}",
        label,
        output,
    ]
    with Popen(command) as m:
        m.communicate()
    return output


def arrangeImages(
    images: list, output: str, mode: str, remove: bool = False
) -> str:
    if mode == "horizontal":
        hori = ["convert", "+append"]
    elif mode == "vertical":
        hori = ["convert", "-append"]
    hori.extend(images)
    hori.append(output)
    with Popen(hori) as m:
        m.communicate()
    if remove:
        for i in images:
            Path(i).unlink()
    return output


if "Bio_SDD" in str(Path(".").absolute()):
    wd = "/home/shannc/Bio_SDD/MUIC_senior_project/workflow/tests"
else:
    wd = "/home/shannc/workflow/tests"
outdir = Path(f"{wd}/testthat/output/protein_dr")
temp = outdir.joinpath("temp")

# Arrange pca and pcoa
pca = []
pcoa = []
for tech, lst in zip(["pca", "pcoa"], [pca, pcoa]):
    for d in outdir.glob(f"{tech}*"):
        if d.is_dir():
            model = re.sub(".*_", "", d.name)
            img = [png for png in d.rglob("*.png")][0]
            lbl = f"Model: {model}"
            labeled = addLabel(str(img), lbl, (600, 60))
            lst.append(labeled)
pca = sorted(pca)
pcoa = sorted(pcoa)

p1 = arrangeImages(pca, f"{outdir}/pca_all.png", "horizontal")
p2 = arrangeImages(pcoa, f"{outdir}/pcoa_all.png", "horizontal")
arrangeImages([p1, p2], f"{outdir}/pca_pcoa.png", "vertical", True)

# Arrange tsne
