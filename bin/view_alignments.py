import string
import intervaltree as it
import polars as pl
import polars.selectors as cs
import sys
import io
import requests
from requests.adapters import Retry, HTTPAdapter
from Bio import SeqIO
import matplotlib as ml
import dna_features_viewer as dv


def parentDir(file_str, parent_level=-1):
    if parent_level > 0:
        raise ValueError("The parent level must be negative!")
    dirs = file_str.split("/")
    return "/".join(dirs[:parent_level])


sys.path.append(parentDir(__file__))
from dna_features_viewer_modified import GraphicRecordCustom

font = {
    "family": "monospace",
    "weight": "medium",
}
ml.rc("font", **font)

POLLING_INTERVAL = 3
API_URL = "https://rest.uniprot.org"
retries = Retry(total=5, backoff_factor=0.25, status_forcelist=[500, 502, 503, 504])
SESSION = requests.Session()
SESSION.mount("https://", HTTPAdapter(max_retries=retries))
FIELDS = "ft_var_seq%2Cft_variant%2Cft_non_cons%2Cft_non_std%2Cft_non_ter%2Cft_conflict%2Cft_unsure%2Cft_act_site%2Cft_binding%2Cft_dna_bind%2Cft_site%2Cft_mutagen%2Cft_intramem%2Cft_topo_dom%2Cft_transmem%2Cft_chain%2Cft_crosslnk%2Cft_disulfid%2Cft_carbohyd%2Cft_init_met%2Cft_lipid%2Cft_mod_res%2Cft_peptide%2Cft_propep%2Cft_signal%2Cft_transit%2Cft_strand%2Cft_helix%2Cft_turn%2Cft_coiled%2Cft_compbias%2Cft_domain%2Cft_motif%2Cft_region%2Cft_repeat%2Cft_zn_fing"


def checkResponse(response):
    try:
        response.raise_for_status()
    except requests.HTTPError:
        print(response)
        raise


def uniprotGetFeatures(uniprot_id: str, file_format: str = "gff") -> str | None:
    url = f"{API_URL}/uniprotkb/{uniprot_id}.{file_format}?fields={FIELDS}"
    request = SESSION.get(url)
    checkResponse(request)
    response = request.content.decode()
    if len(split := response.splitlines()) <= 2 and len(split[1]) <= 3:
        return None
    return response


def readGff(filename: str | io.StringIO) -> pl.DataFrame:
    columns = [
        "seqname",
        "source",
        "feature",
        "start",
        "end",
        "score",
        "strand",
        "frame",
        "attribute",
    ]
    df = pl.read_csv(
        filename,
        separator="\t",
        comment_prefix="#",
        new_columns=columns,
        has_header=False,
    )
    return df.select(cs.by_name(columns))


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


def mergeIntervalData(x, y) -> tuple[list, list]:
    return (x[0] + y[0], x[1] + y[1])


def recordAlignments(alignments, ids=None) -> tuple[pl.DataFrame, it.IntervalTree]:
    ittree = it.IntervalTree()
    a_dict = {
        "full_seq": [],
        "start": [],
        "stop": [],
        "seq": [],
        "id": [],
        "length": [],
    }
    if ids and not len(ids) != len(alignments):
        raise ValueError("Uneven number of ids and alignments!")
    for index, a in enumerate(alignments):
        a_dict["start"].append((start := findChar(a)))
        a_dict["stop"].append((stop := findChar(a, right=True)))
        a_dict["full_seq"].append(a)
        seq = a[start : stop + 1]
        a_dict["seq"].append(seq)
        a_dict["length"].append(len(seq))
        if not ids:
            a_dict["id"].append(f"S{index}")
        else:
            a_dict["id"].append(ids[index])
        ittree[start : stop + 1] = ([a_dict["id"][index]], [seq])

    ittree.merge_overlaps(data_reducer=mergeIntervalData)
    return pl.DataFrame(a_dict), ittree


def mergeContaining(ittree: it.IntervalTree) -> it.IntervalTree:
    removed: set = set()
    merged_tree = it.IntervalTree()
    intervals = ittree.items()
    for i in intervals:
        if i in removed:
            continue
        data = i.data
        enveloped_by = ittree.envelop(i.begin, i.end)
        if not enveloped_by:
            continue
        for e in enveloped_by:
            removed.add(e)
            if e.begin == i.begin and e.end == i.end:
                continue
            data = mergeIntervalData(data, e.data)
        merged_tree[i.begin : i.end] = data
    return merged_tree


def alignmentFeatureList(
    ittree: it.IntervalTree, entry_map: pl.DataFrame
) -> list[dv.GraphicFeature]:
    feature_list = []
    entry_map = entry_map.sort(by="length", descending=True)

    font = {"color": "#cdd6f4", "weight": "bold"}

    for i in ittree:
        start = i.begin
        stop = i.end - 1
        ids, seqs = i.data
        feature = dv.GraphicFeature(
            start=start,
            end=stop,
            strand=0,
            label=f"ALIGNMENTS  n: {len(ids)}",
            fontdict=font,
            color="#1e1e2e",
        )
        feature_list.append(feature)
    return feature_list


def gffFeatureList(gff: pl.DataFrame, colors=list) -> list[dv.GraphicFeature]:
    """
    Create a list of GraphicFeatures from a gff file
    If `colors` is supplied, each unique feature in the gff file is randomly
    given a color from the list
    """
    features: list = []
    color_map = {}
    if colors:
        color_set = set(colors)
    for row in gff.iter_rows():
        gff_feature = row[2]
        label = f"{row[1]}: {gff_feature}"
        if colors and gff_feature not in color_map:
            color_map[gff_feature] = color_set.pop()
            if not color_set:
                color_set = set(colors)
        feature = dv.GraphicFeature(
            start=row[3],
            end=row[4],
            label=label,
            color=color_map.get(gff_feature, "#74c7ec"),
        )
        features.append(feature)
    return features


def generateVisual(
    protein_id: str,
    alignment_df: pl.DataFrame = pl.DataFrame(),
    alignments: list = None,
    sequence_df: pl.DataFrame = None,
    sequence: str = None,
    uniprot_id=None,
):
    if not alignment_df.is_empty():
        found = alignment_df.filter(pl.col("ProteinId") == protein_id)
        alignments = found["alignment"].to_list()
        uniprot_id = found["UniProtKB_ID"][0]
    elif not alignments or not uniprot_id:
        raise ValueError(
            "If no alignment dataframe is provided, the alignments and UniProt Id must be given!"
        )
    if not sequence_df.is_empty():
        sequence = sequence_df.filter(pl.col("ProteinId") == protein_id)["seq"][0]
    elif not sequence:
        raise ValueError(
            "If no sequence dataframe is provided, the sequence must be given!"
        )
    records, interval_tree = recordAlignments(alignments)
    features = alignmentFeatureList(interval_tree, records)
    gff_df = pl.DataFrame()
    if uniprot_id != "NA" and uniprot_id:
        uniprot_request = uniprotGetFeatures(uniprot_id)
        if uniprot_request:
            gff_df = readGff(io.StringIO(uniprot_request))
    if not gff_df.is_empty():
        feature_colors = [
            "#eed49f",
            "#f4dbd6",
            "#fab387",
            "#eba0ac",
            "#f0c6c6",
            "#c6a0f6",
            "#b8c0e0",
            "#cad3f5",
            "#f5a97f",
            "#a6da95",
        ]
        features = gffFeatureList(gff_df, feature_colors) + features
    graphic = GraphicRecordCustom(sequence=sequence, features=features)
    non_polar = {a: "#a6adc8" for a in ["A", "G", "L", "I", "M", "P", "F", "W", "V"]}
    polar_neutral = {a: "#89b4fa" for a in ["T", "C", "N", "Q", "S", "Y"]}
    polar_acidic = {a: "#f38ba8" for a in ["E", "D", "U"]}
    polar_basic = {a: "#cba6f7" for a in ["R", "K", "H", "O"]}
    bg = {**non_polar, **polar_neutral, **polar_acidic, **polar_basic}
    bg["*"] = "#a6e3a1"

    fig, axs = graphic.plot_on_multiple_lines(
        nucl_per_line=100,
        plot_sequence=True,
        figure_width=20,
        sequence_params={"background_map": bg},
    )
    return fig, axs


def parseArgs():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--results_file")
    parser.add_argument("-c", "--coverage_threshold")
    parser.add_argument("-a", "--alignment_file")  # "aligned_peptides.tsv" file
    parser.add_argument("-o", "--outdir")
    args = vars(parser.parse_args())
    return args


def main(args):
    seq_df: pl.DataFrame = (
        pl.read_csv(args["results_file"], separator="\t", null_values="NA")
        .filter(pl.col("pcoverage_align") > args["coverage_threshold"])
        .select(cs.by_name(["ProteinId", "seq"]))
    )
    align_df = pl.read_csv(args["alignment_file"], separator="\t", null_values="NA")
    ids = seq_df["ProteinId"].to_list()
    for id in ids:
        fig, ax = generateVisual(id, alignment_df=align_df, sequence_df=seq_df)
        fig.savefig(f"{args['outdir']}/{id}_alignment.svg", bbox_inches="tight")


if __name__ == "__main__" and not "ipykernel" in sys.argv[0]:
    args = parseArgs()
    main(args)
