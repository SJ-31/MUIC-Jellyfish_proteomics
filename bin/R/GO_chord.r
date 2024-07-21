library("gridGraphics")
library("gridBase")
library("ComplexHeatmap")
library("circlize")

save_legend <- function(
    color_map,
    filename, title, nrows) {
  lgd <- Legend(
    labels = color_map,
    legend_gp = list(fill = names(color_map)), title = title, nrow = nrows
  )
  height <- grobHeight(lgd@grob) |> convertHeight("inches", valueOnly = TRUE)
  width <- grobWidth(lgd@grob) |> convertWidth("inches", valueOnly = TRUE)
  if (str_detect(filename, ".svg$")) {
    saveFun <- \() svg(filename, width = width + 0.5, height = height)
  } else if (str_detect(filename, ".png$")) {
    saveFun <- \() png(filename, width = width + 0.5, height = height, units = "in")
  }
  saveFun()
  grid.draw(lgd)
  dev.off()
}

save_chord <- function(tc, filename, width = 10, height = 10) {
  if (!"from" %in% colnames(tc) || !"to" %in% colnames(tc)) {
    stop("`from` and `to` columns must be specified in given dataframe")
  }
  go_colors <- setNames(tc$from, sample(COLORS, length(tc$from)))

  # Use this setting to control colors based on intensity or smthng
  protein_colors <- setNames(tc$to, sample(COLORS, length(tc$to)))

  colors <- c(go_colors, protein_colors)
  if (str_detect(filename, ".svg$")) {
    ext <- "svg"
    saveFun <- \() svg(filename, width = width + 0.5, height = height)
  } else if (str_detect(filename, ".png$")) {
    ext <- "png"
    saveFun <- \() png(filename, width = width + 0.5, height = height, units = "in")
  }
  saveFun()
  plot_chord(tc, go_colors, TRUE, TRUE)
  dev.off()
  save_legend(
    go_colors,
    glue("{filename}-legend.{ext}"),
    "GO key",
    nrows = 6
  )
}

plot_chord <- function(tc, colors = NULL, show_to = TRUE, show_from = FALSE) {
  if (is.null(colors)) {
    chordDiagramFromDataFrame(as.data.frame(tc),
      directional = 1,
      link.target.prop = TRUE,
      big.gap = 20,
      direction.type = "arrows",
      link.arr.type = "big.arrow",
      preAllocateTracks = list(track.hieght = 5),
      annotationTrack = "grid",
      transparency = 0.8
    )
  } else {
    chordDiagramFromDataFrame(as.data.frame(tc),
      directional = 1,
      link.target.prop = TRUE,
      big.gap = 20,
      direction.type = "arrows",
      link.arr.type = "big.arrow",
      preAllocateTracks = list(track.hieght = 5),
      annotationTrack = "grid",
      grid.col = colors,
      transparency = 0.8
    )
  }
  circos.trackPlotRegion(
    track.index = 1,
    panel.fun = function(x, y) {
      xlim <- get.cell.meta.data("xlim")
      ylim <- get.cell.meta.data("ylim")
      sector.name <- get.cell.meta.data("sector.index")
      if ((sector.name %in% tc$to && show_to) || (sector.name %in% tc$from && show_from)) {
        circos.text(xlim, ylim, NULL)
        circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
        circos.axis(h = "top", labels.cex = 0.5, major.tick.length = 0.2, sector.index = sector.name, track.index = 2)
      }
    }, bg.border = NA
  )
}
