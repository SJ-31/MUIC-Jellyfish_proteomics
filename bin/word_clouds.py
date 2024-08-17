import wordcloud as wc
from matplotlib.patches import Patch
from matplotlib.colors import Colormap
from matplotlib.axes import Axes
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import polars as pl

import matplotlib

matplotlib.rcParams["font.family"] = "sans-serif"


def get_cmap(colormap: str):
    if "ch:" in colormap:
        return sns.color_palette(colormap, as_cmap=True)
    elif "#" in colormap:
        return sns.dark_palette(colormap, as_cmap=True, reverse=True)
    return mpl.colormaps.get(colormap)


class ColorMapper:
    def __init__(
        self,
        items: dict,
        colormap: str = "magma",
        category2colormap: dict = None,
        item2category: dict = None,
    ) -> None:
        """
        :param: a dictionary mapping items to values that determine the intensity of their color
        :param: colormap A valid sequential colormap
        :param: category_mapping dict of category->
        """
        if category2colormap and not item2category:
            raise ValueError(
                "If a category->color is given, then you must specify a mapping of item->category"
            )
        vmin, vmax = min(items.values()), max(items.values())
        self.normalizer = mpl.colors.Normalize(vmin, vmax, clip=False)
        self.items = items
        self.default_cm: str = colormap
        self.category2colormap = category2colormap
        self.all_colormaps: dict[str, Colormap] = {colormap: get_cmap(colormap)}
        if category2colormap:
            self.item2colormap: dict[str, str] = {}
            for i, c in item2category.items():
                cm = category2colormap[c]
                self.item2colormap[i] = cm
                if c not in self.all_colormaps:
                    self.all_colormaps[cm] = get_cmap(cm)
            print("Initiated the following colormaps:\n")
            print(self.all_colormaps)
        else:
            self.item2colormap = None

    def add_cmap_legend(self, ax: Axes, **kwargs):
        if not self.item2colormap:
            return
        for cat, cm_name in self.category2colormap.items():
            cm: Colormap = self.all_colormaps.get(cm_name)
            sm = mpl.cm.ScalarMappable(norm=self.normalizer, cmap=cm)
            bar = plt.colorbar(
                sm,
                ax=ax,
                location="right",
                ticks=None,
                label=cat,
                **kwargs,
            )

    def __call__(self, item: str, **kwargs) -> str:
        if self.item2colormap:
            cm: Colormap = self.all_colormaps.get(self.item2colormap[item])
        else:
            cm: Colormap = self.all_colormaps.get(self.default_cm)
        val = self.normalizer(self.items[item])
        return mpl.colors.to_hex(cm(val))


def cols2dict(df: pl.DataFrame, keys: str, vals: str) -> dict:
    return dict(zip(df[keys], df[vals]))


def add_abbrev_legend(
    ax: Axes, abbrevs: dict, size: int = 15, location=(0.2, 1)
) -> None:
    """
    Add abbreviation annotations to an empty axes
    """
    start = -(size + 5)
    header: mpl.text.Annotation = ax.annotate(
        "Abbreviations",
        location,
        xycoords="axes fraction",
        ha="center",
        weight="bold",
        size=size,
    )
    for k, v in abbrevs.items():
        anno_size = size - 2
        lab = ax.annotate(
            k,
            xy=(0.2, 0),
            xycoords=header,
            xytext=(0, start),
            textcoords="offset points",
            size=anno_size,
            weight="bold",
        )
        ax.annotate(
            f"= {v}",
            xy=(0.55, 0),
            xycoords=header,
            xytext=(0, start),
            textcoords="offset points",
            size=anno_size,
        )
        start -= size
    return lab


def word_cloud_main(tokens: dict, abbrevs: dict = None, params: dict = {}):
    """
    Generates a word cloud image from the given tokens and optionally includes a legend for abbreviations.

    Parameters
    ----------
    tokens : dict
        A dictionary where keys are words or phrases and values are their respective frequencies.
    abbrevs : dict
        A dictionary mapping abbreviations to their full forms for the legend. If empty, no legend is added.
    filename : str
        The path where the generated word cloud image will be saved.
    params : dict, optional
        A dictionary of optional parameters for customization. Supported keys include:

        - colormap: A colormap to use for coloring words.
        - category2colormap: A mapping from categories to colormaps.
        - item2category: A mapping from items to categories.
        - background_color: Background color of the word cloud image (default is "white").
        - width: Width of the word cloud image in pixels (default is 1000).
        - height: Height of the word cloud image in pixels (default is 1000).
        - min_font_size: Minimum font size for words in the word cloud (default is 10).
        - font: Path to the font file to be used for the word cloud (default is "/home/shannc/.fonts/FiraSans-Regular.ttf").
        - interpolation: The interpolation method used for displaying the image (default is "bicubic").
        - title: Title to be displayed above the word cloud.
        - title_size: Font size for the title (default is 20).
        - abbrev_size: Font size for the abbreviation legend (default is 15).
        - cb_fraction: Fraction of original axes for colorbar (default 0.10)
        - cb_shrink: Fraction by which to multiply size of colorbar (default 1.5)
        - fig_size: Figure size as a tuple (width, height) in inches (default is (15, 15)).
    """

    color_mapper: ColorMapper = ColorMapper(
        tokens,
        colormap=params.get("colormap", "magma"),
        category2colormap=params.get("category2colormap"),
        item2category=params.get("item2category"),
    )
    cloud: wc.WordCloud = wc.WordCloud(
        background_color=params.get("background_color", "white"),
        width=params.get("width", 1000),
        height=params.get("height", 1000),
        min_font_size=params.get("min_font_size", 10),
        font_path=params.get("font", "/home/shannc/.fonts/FiraSans-Regular.ttf"),
        color_func=color_mapper,
    )
    image = cloud.generate_from_frequencies(tokens)
    if abbrevs:
        fig, ax = plt.subplots(ncols=2, width_ratios=[3, 1], layout="constrained")
        ax[0].imshow(image, interpolation=params.get("interpolation", "bicubic"))
        if title := params.get("title"):
            ax[0].set_title(title, weight="bold", size=params.get("title_size", 20))
        ax[1].set_aspect("equal")
        add_abbrev_legend(ax[1], abbrevs, size=params.get("abbrev_size", 15))
        color_mapper.add_cmap_legend(
            ax[1],
            fraction=params.get("cb_fraction", 0.15),
            shrink=params.get("cb_shrink", 1),
        )
    else:
        fig, ax = plt.subplots(layout="constrained")
        color_mapper.add_cmap_legend(
            ax,
            fraction=params.get("cb_fraction", 0.10),
            shrink=params.get("cb_shrink", 1.5),
        )
    for a in ax:
        a.axis("off")
    fig.set_size_inches(params.get("fig_size", (15, 15)))
    return fig
