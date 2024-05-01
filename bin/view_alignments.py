import os
import string
import re
from subprocess import Popen, PIPE, CalledProcessError
from Bio import SeqIO
import dna_features_viewer as dv
import pandas as pd
import matplotlib as ml

font = {
    "family": "Source Code Pro",
    "weight": "medium",
}
ml.rc("font", **font)


class ReadTemporary:
    """
    Initializes with a command that creates some temporary file, reads from it,
    and deletes it afterward
    """

    def __init__(self, command: list, temp_file: str = "temp"):
        self.temp_file = temp_file
        self.command = command

    def __enter__(self):
        with Popen(self.command, stderr=PIPE, stdout=PIPE) as c:
            result = c.communicate()
            returncode = c.returncode
        print(result[0].decode())
        print(result[1].decode())
        if returncode == 0:
            return self.temp_file
        raise CalledProcessError(returncode=returncode, cmd=self.command)

    def __exit__(self, exc_type, exc_val, exc_tb):
        os.remove(self.temp_file)


def findChar(given: str, right: bool = False) -> int:
    length: int = len(given)
    start, end, step = (-1, -length - 1, -1) if right else (0, length, 1)
    found_index = -1
    for i in range(start, end, step):
        if given[i] in string.ascii_letters:
            found_index = i
            break
    if right:
        found_index = length + found_index
    return found_index


def fasta2Dict(filename) -> dict[str, str]:
    dct = {}
    for record in SeqIO.parse(filename, "fasta"):
        dct[record.id] = str(record.seq)
    return dct


def writeFasta(filename, seqs, headers=None) -> dict[str, str]:
    fasta: list = []
    fdict: dict = {}
    if not headers:
        count = 0
        for seq in seqs:
            header = f"S{count}"
            fasta.append(f">{header}\n{seq}")
            fdict[header] = seq
            count += 1
    else:
        for seq, header in zip(seqs, headers):
            fasta.append(f">{header}\n{seq}")
            fdict[header] = seq
    with open(filename, "w") as w:
        w.write("\n".join(fasta))
    return fdict


def recordAlignments(alignments, ids=None) -> tuple[pd.DataFrame, int]:
    a_dict = {"full_seq": [], "start": [], "stop": [], "seq": [], "id": []}
    if ids and not len(ids) != len(alignments):
        raise ValueError("Uneven number of ids and alignments!")
    for index, a in enumerate(alignments):
        a_dict["start"].append((start := findChar(a)))
        a_dict["stop"].append((stop := findChar(a, right=True)))
        a_dict["full_seq"].append(a)
        a_dict["seq"].append(a[start : stop + 1])
        if not ids:
            a_dict["id"].append(f"S{index}")
        else:
            a_dict["id"].append(ids[index])

    return pd.DataFrame(a_dict), len(alignments)


def getFeatureList(entry_map) -> list[dv.GraphicFeature]:
    feature_list = []
    covered_sequences = {}

    def isWithin() -> bool:
        return False

    seen_seqs = set()
    for i in range(len(entry_map)):
        start = entry_map["start"][i]
        stop = entry_map["stop"][i]
        seq = entry_map["seq"][i]
        if seq in seen_seqs:
            continue
        else:
            seen_seqs.add(seq)
        feature = dv.GraphicFeature(start=start, end=stop, strand=0, label=seq)
        feature_list.append(feature)
    return feature_list


alignments = []
with open("/home/shannc/Dropbox/testzone/alignments.txt", "r") as t:
    seq = next(t).strip()
    for l in t:
        alignments.append(l.strip())

OUTDIR = "/home/shannc/Downloads/thesis_testzone"
records, num_alignments = recordAlignments(alignments)


fd = writeFasta(f"{OUTDIR}/alignments.fasta", records["seq"], list(records["id"]))

cdhit_path = "/home/shannc/Bio_SDD/tools/cd-hit-v4.8.1-2019-0228/cd-hit"
with ReadTemporary(
    command=[cdhit_path, "-i", f"{OUTDIR}/alignments.fasta", "-o", "temp"],
) as cdhit:
    clustered = fasta2Dict(cdhit)


representative_seqs = records[records["id"].isin(clustered.keys())].reset_index()
alignment_features = getFeatureList(representative_seqs)
graphic = dv.GraphicRecord(sequence=seq, features=alignment_features)

fig, axs = graphic.plot_on_multiple_lines(
    nucl_per_line=100, plot_sequence=True, figure_width=20
)
fig.savefig(f"{OUTDIR}/sequence_and_translation.svg", bbox_inches="tight")
