import plotly.graph_objects as go
import plotly.express as px
import pandas as pd
import matplotlib.pyplot as plt


def plotly_save(figure, metadata) -> None:
    file = metadata.get("file")
    if file and "html" in file:
        figure.write_html(file)
    elif file:
        figure.write_image(
            file,
            height=metadata.get("height", 1000),
            width=metadata.get("width", 1500),
            scale=3,
        )


def plotly_sunburst(sunburst, metadata: dict, layout_params: dict = None):
    fig = go.Figure(
        go.Sunburst(
            labels=sunburst["name"],
            parents=sunburst["parents"],
            ids=sunburst["ids"],
            values=sunburst["values"],
            branchvalues="total",
            outsidetextfont={"weight": "bold"},
            textinfo="label+value",
        )
    )
    if layout_params:
        fig.update_layout(**layout_params)
    plotly_save(fig, metadata)


def plotly_treemap(treemap, metadata, layout_params: dict = None):
    fig = go.Figure(
        go.Treemap(
            labels=treemap["name"],
            parents=treemap["parents"],
            ids=treemap["ids"],
            values=treemap["values"],
            branchvalues="total",
            outsidetextfont={"weight": "bold"},
            textinfo="label+value",
        )
    )
    if layout_params:
        fig.update_layout(**layout_params)
    plotly_save(fig, metadata)


def plotly_psm_comparisons(df, compare_col: str):
    left, right = f"{compare_col}.first", f"{compare_col}.sec"
    if compare_col == "e_value":
        color = "E-value lower in First pass"
        df[color] = df[left] < df[right]
    elif compare_col == "psm_score":
        color = "PSM score higher in First"
        df[color] = df[left] > df[right]
    fig = px.scatter(
        df, y=left, x=right, facet_row="engine", facet_col="param", color=color
    )
    fig.update_yaxes(matches=None)
    fig.update_xaxes(matches=None)
    fig.show()
    return fig
