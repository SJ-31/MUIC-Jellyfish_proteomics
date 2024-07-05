#!/usr/bin/env python

from matplotlib.collections import PatchCollection
from matplotlib.axes import Axes
from Bio import Align
from pathlib import Path
import pymsaviz as pv
from matplotlib.patches import Rectangle
import string
from Bio.SeqRecord import SeqRecord
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

non_polar = {a: "#a6adc8" for a in ["A", "G", "L", "I", "M", "P", "F", "W", "V"]}
polar_neutral = {a: "#89b4fa" for a in ["T", "C", "N", "Q", "S", "Y"]}
polar_acidic = {a: "#f38ba8" for a in ["E", "D", "U"]}
polar_basic = {a: "#cba6f7" for a in ["R", "K", "H", "O"]}
COLOR_SCHEME = {**non_polar, **polar_neutral, **polar_acidic, **polar_basic}


def parent_dir(file_str, parent_level=-1):
    if parent_level > 0:
        raise ValueError("The parent level must be negative!")
    dirs = file_str.split("/")
    return "/".join(dirs[:parent_level])


sys.path.append(parent_dir(__file__))
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


def check_response(response):
    try:
        response.raise_for_status()
    except requests.HTTPError:
        print(response)
        raise


def get_uniprot_features(uniprot_id: str, file_format: str = "gff") -> str | None:
    url = f"{API_URL}/uniprotkb/{uniprot_id}.{file_format}?fields={FIELDS}"
    request = SESSION.get(url)
    check_response(request)
    response = request.content.decode()
    if len(split := response.splitlines()) <= 2 and len(split[1]) <= 3:
        return None
    return response


def read_gff(filename: str | io.StringIO) -> pl.DataFrame:
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


def find_char(given: str, right: bool = False) -> int:
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


def fasta2dict(filename) -> dict[str, str]:
    dct = {}
    for record in SeqIO.parse(filename, "fasta"):
        dct[record.id] = str(record.seq)
    return dct


def write_fasta(filename, seqs, headers=None) -> dict[str, str]:
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


def merge_interval_data(x, y) -> tuple[list, list]:
    return (x[0] + y[0], x[1] + y[1])


def record_alignments(alignments, ids=None) -> tuple[pl.DataFrame, it.IntervalTree]:
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
        a_dict["start"].append((start := find_char(a)))
        a_dict["stop"].append((stop := find_char(a, right=True)))
        a_dict["full_seq"].append(a)
        seq = a[start : stop + 1]
        a_dict["seq"].append(seq)
        a_dict["length"].append(len(seq))
        if not ids:
            a_dict["id"].append(f"S{index}")
        else:
            a_dict["id"].append(ids[index])
        ittree[start : stop + 1] = ([a_dict["id"][index]], [seq])

    ittree.merge_overlaps(data_reducer=merge_interval_data)
    return pl.DataFrame(a_dict), ittree


def merge_containing(ittree: it.IntervalTree) -> it.IntervalTree:
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
            data = merge_interval_data(data, e.data)
        merged_tree[i.begin : i.end] = data
    return merged_tree


def alignment_feature_list(
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


def gff_feature_list(gff: pl.DataFrame, colors=list) -> list[dv.GraphicFeature]:
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


def add_uniprot_features(uniprot_id, features: list) -> list:
    gff_df = pl.DataFrame()
    if uniprot_id != "NA" and uniprot_id:
        uniprot_request = get_uniprot_features(uniprot_id)
        if uniprot_request:
            gff_df = read_gff(io.StringIO(uniprot_request))
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
        features = gff_feature_list(gff_df, feature_colors) + features
    return features


def generate_visual(
    protein_id: str,
    alignment_df: pl.DataFrame = pl.DataFrame(),
    alignments: list = None,
    sequence_df: pl.DataFrame = None,
    sequence: str = None,
    uniprot_id=None,
    add_uniprot=False,
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
    records, interval_tree = record_alignments(alignments)
    features = alignment_feature_list(interval_tree, records)
    if (uniprot_id != "NA" and uniprot_id) and add_uniprot:
        features = add_uniprot_features(uniprot_id, features)
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


class PeptideViz(pv.MsaViz):
    """Class for representing peptides aligned to a central protein sequence
    Code adapted from https://github.com/moshi4/pyMSAviz
    """

    def __init__(
        self,
        msa: str | Path | Align.MultipleSeqAlignment,
        *,
        aligned_to: SeqRecord = None,
        format: str = "fasta",
        color_scheme: str | None = None,
        start: int = 1,
        end: int | None = None,
        wrap_length: int | None = None,
        wrap_space_size: float = 3.0,
        show_label: bool = True,
        label_type: str = "id",
        show_seq_char: bool = True,
        show_grid: bool = False,
        show_count: bool = False,
        unknown_char="-",
        show_consensus: bool = False,
        consensus_color: str = "#1f77b4",
        consensus_size: float = 2.0,
        sort: bool = False,
    ):
        super().__init__(
            msa,
            format=format,
            color_scheme=color_scheme,
            start=start,
            end=end,
            wrap_length=wrap_length,
            wrap_space_size=wrap_space_size,
            show_label=show_label,
            label_type=label_type,
            show_seq_char=show_seq_char,
            show_grid=show_grid,
            show_count=show_count,
            show_consensus=show_consensus,
            consensus_color=consensus_color,
            consensus_size=consensus_size,
            sort=sort,
        )
        self.unknown_char = unknown_char
        self.aligned_to: SeqRecord = aligned_to
        self.matched_positions: set = set()
        for i in range(msa.get_alignment_length()):
            if set(msa[:, i]) | {"X", "-"} != {"X", "-"}:
                self.matched_positions.add(i)

    def _plot_consensus(
        self, ax: Axes, start: int | None = None, end: int | None = None
    ) -> None:
        """Plot the sequence of the complete protein that the peptides
        were aligned to
        """
        # Set xlim, ylim
        start = 0 if start is None else start
        end = self.alignment_length if end is None else end
        ax.set_xlim(start, start + self._wrap_length)
        ax.set_ylim(0, 80)  # 0 - 100 [%]

        # Plot label text
        y_lower = 50
        if self._show_label and self._consensus_size != 0:
            ax.text(
                start - 1, y_lower, self.aligned_to.id, ha="right", va="center", size=10
            )

        # Set spines & tick params
        for pos in ("left", "right", "top", "bottom"):
            ax.spines[pos].set_visible(False)
        ax.tick_params(bottom=False, left=False, labelleft=False, pad=0)
        ax.axis("off")

        plot_patches = []
        y_center = y_lower + 0.5
        seq = self.aligned_to.seq
        for x_left in range(start, end):
            seq_char = seq[x_left]
            x_center = x_left + 0.5
            ax.text(
                x_center,
                y_center,
                seq_char,
                ha="center",
                va="center",
                size=10,
            )
            if x_left in self.matched_positions:
                # Highlight residues that have at least one peptide residue aligned
                rect_prop: dict = dict(
                    xy=(x_left, 32), width=1, height=40, color="none", lw=0
                )
                highlight_positions = self._highlight_positions
                if highlight_positions is None or x_left in highlight_positions:
                    color = self.color_scheme.get(seq_char, "#FFFFFF")
                    if self._color_scheme_name == "Identity":
                        color = self._get_identity_color(seq_char, x_left)
                    if self._custom_color_func is not None:
                        custom_color = self._custom_color_func(
                            0, x_left, seq_char, self.msa
                        )
                        color = color if custom_color is None else custom_color
                    rect_prop.update(**dict(color=color, lw=0, fill=True))
                plot_patches.append(Rectangle(**rect_prop))

        # Plot colored rectangle patch collection (Use collection for speedup)
        collection = PatchCollection(plot_patches, match_original=True, clip_on=False)
        ax.add_collection(collection)  # type: ignore

    def _plot_msa(
        self, ax: Axes, start: int | None = None, end: int | None = None
    ) -> None:
        """Plot MSA

        Parameters
        ----------
        ax : Axes
            Matplotlib axes to be plotted
        start : int | None, optional
            Start position. If None, `0` is set.
        end : int | None, optional
            End position. If None, `alignment_length` is set.
        """
        # Set xlim, ylim
        start = 0 if start is None else start
        end = self.alignment_length if end is None else end
        ax.set_xlim(start, start + self._wrap_length)
        ax.set_ylim(0, self.msa_count)

        # Set spines & tick params (Only show bottom ticklables)
        for pos in ("left", "right", "top", "bottom"):
            ax.spines[pos].set_visible(False)
        ax.tick_params(left=False, labelleft=False)

        # Plot alignment position every 10 chars on xticks
        ticks_interval = self._ticks_interval
        if ticks_interval is None:
            ax.tick_params(bottom=False, labelbottom=False)
        else:
            tick_ranges = range(start + 1, end + 1)
            xticklabels = list(filter(lambda n: n % ticks_interval == 0, tick_ranges))
            xticks = [n - 0.5 for n in xticklabels]
            ax.set_xticks(xticks, xticklabels, size=8)  # type: ignore

        plot_patches = []
        for cnt in range(self.msa_count):
            msa_seq = self.seq_list[cnt]
            y_lower = self.msa_count - (cnt + 1)
            y_center = y_lower + 0.5
            # Plot label text
            if self._show_label:
                if self._label_type == "id":
                    label = self.id_list[cnt]
                elif self._label_type == "description":
                    label = self.desc_list[cnt]
                else:
                    err_msg = f"{self._label_type=} is invalid (`id`|`description`)"
                    raise ValueError(err_msg)
                ax.text(
                    start - 1,
                    y_center,
                    label,
                    ha="right",
                    va="center",
                    size=10,
                )
            # Plot count text
            if self._show_count:
                scale = end - self._start - msa_seq[self._start : end].count("-")
                ax.text(
                    end + 1,
                    y_center,
                    str(scale),
                    ha="left",
                    va="center",
                    size=10,
                )
            for x_left in range(start, end):
                # Add colored rectangle patch
                seq_char = msa_seq[x_left]
                rect_prop: dict = dict(
                    xy=(x_left, y_lower), width=1, height=1, color="none", lw=0
                )
                highlight_positions = self._highlight_positions
                if highlight_positions is None or x_left in highlight_positions:
                    color = self.color_scheme.get(seq_char, "#FFFFFF")
                    if self._color_scheme_name == "Identity":
                        color = self._get_identity_color(seq_char, x_left)
                    if self._custom_color_func is not None:
                        custom_color = self._custom_color_func(
                            cnt, x_left, seq_char, self.msa
                        )
                        color = color if custom_color is None else custom_color
                    rect_prop.update(**dict(color=color, lw=0, fill=True))
                if self._show_grid:
                    rect_prop.update(**dict(ec=self._grid_color, lw=0.5))
                plot_patches.append(Rectangle(**rect_prop))

                # Plot seq char text
                if seq_char == "X" or seq_char == "-":
                    seq_char = self.unknown_char
                x_center = x_left + 0.5
                if self._show_seq_char:
                    ax.text(
                        x_center,
                        y_center,
                        seq_char,
                        ha="center",
                        va="center",
                        size=10,
                    )
                # Plot marker
                if cnt == 0 and x_left in self._pos2marker_kws:
                    marker_kws = self._pos2marker_kws[x_left]
                    ax.plot(x_center, y_center + 1, **marker_kws)
                # Plot text annotation
                if cnt == 0 and x_left in self._pos2text_kws:
                    text_kws = self._pos2text_kws[x_left]
                    ax.text(**text_kws)

        # Plot colored rectangle patch collection (Use collection for speedup)
        collection = PatchCollection(plot_patches, match_original=True, clip_on=False)
        ax.add_collection(collection)  # type: ignore


def main(args):
    seq_df: pl.DataFrame = (
        pl.read_csv(args["results_file"], separator="\t", null_values="NA")
        .filter(pl.col("pcoverage_align") > args["coverage_threshold"])
        .select(cs.by_name(["ProteinId", "seq"]))
    )
    align_df = pl.read_csv(args["alignment_file"], separator="\t", null_values="NA")
    ids = seq_df["ProteinId"].to_list()
    if args["mode"] == "engine_alignment":
        from trace_alignments import AlignmentTracer

        tracer = AlignmentTracer(args["alignment_file"], args["peptide_map_file"])
        tracer.get_id2metadata(args["results_file"])
        for id in ids:
            tracer.plot_engines_alignment(id, args["outdir"])
    elif args["mode"] == "with_uniprot":
        for id in ids:
            fig, ax = generate_visual(id, alignment_df=align_df, sequence_df=seq_df)
            fig.savefig(f"{args['outdir']}/{id}_alignment.svg", bbox_inches="tight")


def parse_args():
    import argparse

    parser = argparse.ArgumentParser()
    parser.add_argument("-r", "--results_file")
    parser.add_argument("-c", "--coverage_threshold", type=float)
    parser.add_argument("-a", "--alignment_file")  # "aligned_peptides.tsv" file
    parser.add_argument("-p", "--peptide_map_file")
    parser.add_argument("-m", "--mode")
    parser.add_argument("-o", "--outdir")
    args = vars(parser.parse_args())
    return args


if __name__ == "__main__" and not "ipykernel" in sys.argv[0]:
    args = parse_args()
    main(args)
