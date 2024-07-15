import plotly.graph_objects as go


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
