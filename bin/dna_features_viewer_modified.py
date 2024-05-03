import matplotlib.pyplot as plt
import numpy
import dna_features_viewer as dv


class GraphicRecordCustom(dv.GraphicRecord):
    def __init__(
        self,
        sequence_length=None,
        sequence=None,
        features=...,
        feature_level_height=1,
        first_index=0,
        plots_indexing="biopython",
        labels_spacing=8,
        ticks_resolution="auto",
    ):
        super().__init__(
            sequence_length,
            sequence,
            features,
            feature_level_height,
            first_index,
            plots_indexing,
            labels_spacing,
            ticks_resolution,
        )

    def plot_sequence(
        self,
        ax,
        location=None,
        y_offset=1,
        fontdict=None,
        guides_intensity=0,
        color_map: dict = None,
        background_map: dict = None,
    ):
        """Plot a sequence of nucleotides at the bottom of the plot.

        Parameters
        ----------

        ax
        Which axes the translation should be plotted to.

        location
        location of the segment to translate, either (start, end) or
        (start, end, strand).

        y_offset
        Number of text levels under the plot's base line where to draw the
        nucleotides. Should be 1 if the nucleotide sequence is to be plotted
        directly under the main line.

        fontdict
        Matplotlib fontdict for the text, e.g.
        ``{'size': 11, 'weight':'bold'}``

        color_map
        Dictionary of colors (following matplotlib specification)
        for each symbol e.g. {"A": "r", "T": "b", "*": "#0f0f0f80"}.
        "*" is the default color applied to symbols that aren't specified
        Leave to None for no background.

        background_map
        Dictionary of background colors (following matplotlib specification)
        for the background of each symbol e.g. {"A": "r", "T": "b", "*": "#0f0f0f80"}.
        "*" is the default color applied to symbols that aren't specified
        Leave to None for no background.

        guides_intensity
        Intensity of the vertical guides marking the different nucleotides
        (0 = no guides).
        """
        if color_map and "*" not in color_map:
            raise ValueError("Default color must be specified!")
        if background_map and "*" not in background_map:
            raise ValueError("Default background must be specified!")

        def get_symbol_color(symbol):
            color_dict = {}
            if color_map:
                color = color_map.get(symbol, color_map.get("*"))
                color_dict = {"color": color}
            if background_map:
                bg = background_map.get(symbol, background_map.get("*"))
                color_dict = {**color_dict, **{"backgroundcolor": bg}}
            return color_dict

        if self.sequence is None:
            raise ValueError("No sequence in the graphic record")
        if location is None:
            location = self.span
        location_start, location_end = location
        fontdict = {"size": 11, **(fontdict or {})}
        for i, nucleotide in enumerate(self.sequence):
            color_dict = get_symbol_color(nucleotide)
            index = i + location_start
            if location_start <= index <= location_end:
                ax.text(
                    index,
                    -0.7 * self.feature_level_height * y_offset,
                    nucleotide,
                    ha="center",
                    va="center",
                    fontdict=fontdict,
                    **color_dict,
                )
        if guides_intensity:
            color = (0, 0, 0, guides_intensity)
            for i in range(location_start, location_end + 1):
                ax.axvline(i - 0.5, linewidth=0.1, color=color, zorder=-10000)
        ymin = ax.get_ylim()[0]
        if ymin < -500:
            ymin = 0
        ax.set_ylim(bottom=min(ymin, -y_offset * self.feature_level_height))

    def crop(self, window):
        start, end = window
        first_index = self.first_index
        if (start < first_index) or (end > self.last_index):
            raise ValueError("out-of-bound cropping")
        new_features = []
        for f in self.features:
            cropped_feature = f.crop(window)
            if cropped_feature is not None:  # = has ovelap with the window
                new_features.append(cropped_feature)

        return GraphicRecordCustom(
            sequence=(
                self.sequence[start - first_index : end - first_index]
                if self.sequence is not None
                else None
            ),
            sequence_length=end - start,
            features=new_features,
            feature_level_height=self.feature_level_height,
            first_index=start,
            plots_indexing=self.plots_indexing,
            labels_spacing=self.labels_spacing,
            ticks_resolution=self.ticks_resolution,
        )

    # Override
    def plot_on_multiple_lines(
        self,
        n_lines=None,
        nucl_per_line=None,
        plot_sequence=False,
        figure_width="auto",
        translation_params=None,
        sequence_params=None,
        **plot_params,
    ):
        if n_lines is None:
            n_lines = int(numpy.ceil(self.sequence_length / nucl_per_line))
        else:
            nucl_per_line = self.sequence_length // n_lines + 1
        if figure_width == "auto":
            if plot_sequence:
                figure_width = 0.15 * nucl_per_line
            else:
                figure_width = 10

        figures_heights = []

        def plot_line(line_index, ax=None, translation_params=None):
            first, last = self.first_index, self.last_index
            line_start = first + line_index * nucl_per_line
            line_virtual_end = first + (line_index + 1) * nucl_per_line
            line_end = min(last, line_virtual_end)
            line_record = self.crop((line_start, line_end))
            line_ax, _ = line_record.plot(
                figure_width=figure_width,
                x_lim=(line_start, line_virtual_end),
                ax=ax,
                plot_sequence=plot_sequence,
                sequence_params=sequence_params,
                **plot_params,
            )
            if translation_params is not None:
                line_record.plot_translation(ax=line_ax, **translation_params)
            return line_ax

        for line_index in range(n_lines):
            line_ax = plot_line(line_index, translation_params=translation_params)
            figures_heights.append(line_ax.figure.get_figheight())
            plt.close(line_ax.figure)
        fig, axes = plt.subplots(
            n_lines,
            1,
            gridspec_kw={"height_ratios": figures_heights},
            figsize=(figure_width, 0.9 * sum(figures_heights)),
        )
        if n_lines == 1:
            axes = [axes]
        for line_index, ax in enumerate(axes):
            plot_line(line_index, ax=ax, translation_params=translation_params)
        fig.tight_layout()
        return fig, axes
