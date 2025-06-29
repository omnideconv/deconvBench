#' ggradar
#'
#' @param plot.data dataframe comprising one row per group
#' @param base.size text size
#' @param font.radar text font family
#' @param values.radar values to print at minimum, 'average', and maximum gridlines
#' @param axis.labels  names of axis labels if other than column names supplied via plot.data
#' @param grid.min value at which mininum grid line is plotted
#' @param grid.mid value at which 'average' grid line is plotted
#' @param grid.max value at which maximum grid line is plotted
#' @param centre.y value of y at centre of plot
#' @param plot.extent.x.sf controls relative size of plot horizontally
#' @param plot.extent.y.sf controls relative size of plot vertically
#' @param x.centre.range controls axis label alignment
#' @param label.centre.y whether value of y at centre of plot should be labelled
#' @param grid.line.width width of gridline
#' @param gridline.min.linetype line type of minimum gridline
#' @param gridline.mid.linetype line type of 'average' gridline
#' @param gridline.max.linetype line type of maximum gridline
#' @param gridline.min.colour colour of minimum gridline
#' @param gridline.mid.colour colour of 'average' gridline
#' @param gridline.max.colour colour of maximum gridline
#' @param grid.label.size text size of gridline label
#' @param gridline.label.offset displacement to left/right of central vertical axis
#' @param label.gridline.min whether or not to label the mininum gridline
#' @param label.gridline.mid whether or not to label the 'mininum'average' gridline
#' @param label.gridline.max whether or not to label the maximum gridline
#' @param axis.label.offset vertical displacement of axis labels from maximum grid line, measured relative to circle diameter
#' @param axis.label.size text size of axis label
#' @param axis.line.colour colour of axis line
#' @param group.line.width line width of group
#' @param group.point.size point size of group
#' @param group.colours colour of group
#' @param background.circle.colour colour of background circle/radar
#' @param background.circle.transparency transparency of background circle/radar
#' @param plot.legend whether to include a plot legend
#' @param legend.title title of legend
#' @param plot.title title of radar plot
#' @param legend.text.size text size in legend
#' @param legend.position position of legend, valid values are "top", "right", "bottom", "left"
#' @param fill whether to fill polygons
#' @param fill.alpha if filling polygons, transparency values
#' @param draw.points whether to draw points
#' @param point.alpha alpha for points, can be a single value or vector
#' @param line.alpha alpha for lines, can be a single value or vector
#'
#' @import ggplot2
#' @return a ggplot object
#'
#' @name ggradar-package
#'
#' @export
#'
#' @source
#' Most of the code is from \url{http://rstudio-pubs-static.s3.amazonaws.com/5795_e6e6411731bb4f1b9cc7eb49499c2082.html}.
#'
#' @examples
#' \dontrun{
#' library(ggradar)
#' library(dplyr)
#' library(scales)
#' library(tibble)
#'
#' mtcars_radar <- mtcars %>%
#'   as_tibble(rownames = "group") %>%
#'   mutate_at(vars(-group), rescale) %>%
#'   tail(4) %>%
#'   select(1:10)
#' mtcars_radar
#' ggradar(mtcars_radar)
#' }
ggradar_custom <- function(plot.data,
                           base.size = 15,
                           font.radar = "sans",
                           values.radar = c("0%", "50%", "80%","100%"),
                           axis.labels = colnames(plot.data)[-1],
                           grid.min = 0, # 10,
                           grid.mid = 0.5, # 50,
                           grid.add = 0.8, # 80,
                           grid.max = 1, # 100,
                           centre.y = grid.min - ((1 / 9) * (grid.max - grid.min)),
                           plot.extent.x.sf = 1,
                           plot.extent.y.sf = 1.2,
                           x.centre.range = 0.02 * (grid.max - centre.y),
                           label.centre.y = FALSE,
                           grid.line.width = 0.5,
                           gridline.min.linetype = "longdash",
                           gridline.mid.linetype = "longdash",
                           gridline.add.linetype = "longdash",
                           gridline.max.linetype = "longdash",
                           gridline.min.colour = "grey",
                           gridline.mid.colour = "#007A87",
                           gridline.add.colour = "#007A87",
                           gridline.max.colour = "grey",
                           grid.label.size = 6,
                           gridline.label.offset = -0.1 * (grid.max - centre.y),
                           label.gridline.min = TRUE,
                           label.gridline.mid = TRUE,
                           label.gridline.add = TRUE,
                           label.gridline.max = TRUE,
                           axis.label.offset = 1.15,
                           axis.label.size = 5,
                           axis.line.colour = "grey",
                           group.line.width = 1.5,
                           group.point.size = 6,
                           group.colours = NULL,
                           background.circle.colour = "#D7D6D1",
                           background.circle.transparency = 0.2,
                           plot.legend = if (nrow(plot.data) > 1) TRUE else FALSE,
                           legend.title = "",
                           plot.title = "",
                           legend.text.size = 14,
                           legend.position = "left",
                           fill = FALSE,
                           fill.alpha = 0.5,
                           draw.points = TRUE, # Whether to draw points
                           point.alpha = 1, # Alpha for points, can be a single value or vector
                           line.alpha = 1 # Alpha for lines, can be a single value or vector
) {
  plot.data <- as.data.frame(plot.data)
  # if there are several groups in the first column with differing values
  # on the dimensions, we should aggregate them by taking the mean, otherwise
  # only the first row is taken into account in the function CalculateGroupPath.
  plot.data <- aggregate(
    x = plot.data[, -1], 
    by = list(plot.data[, 1]), 
    FUN = "mean")
  
  if (!is.factor(plot.data[, 1])) {
    plot.data[, 1] <- as.factor(as.character(plot.data[, 1]))
  }
  
  var.names <- colnames(plot.data)[-1] # Short version of variable names
  # axis.labels [if supplied] is designed to hold 'long version' of variable names
  # with line-breaks indicated using \n
  
  # calculate total plot extent as radius of outer circle x a user-specifiable scaling factor
  plot.extent.x <- (grid.max + abs(centre.y)) * plot.extent.x.sf
  plot.extent.y <- (grid.max + abs(centre.y)) * plot.extent.y.sf
  
  # Check supplied data makes sense
  if (length(axis.labels) != ncol(plot.data) - 1) {
    stop("'axis.labels' contains the wrong number of axis labels", call. = FALSE)
  }
  if (min(plot.data[, -1]) < centre.y) {
    stop("plot.data' contains value(s) < centre.y", call. = FALSE)
  }
  
  if (max(plot.data[, -1]) > grid.max) {
    plot.data[, -1] <- (plot.data[, -1]/max(plot.data[, -1]))*grid.max
    warning("'plot.data' contains value(s) > grid.max, data scaled to grid.max", call. = FALSE)
  }
  
  
  ### Convert supplied data into plottable format
  # (a) add abs(centre.y) to supplied plot data
  # [creates plot centroid of 0,0 for internal use, regardless of min. value of y
  # in user-supplied data]
  plot.data.offset <- plot.data
  plot.data.offset[, 2:ncol(plot.data)] <- plot.data[, 2:ncol(plot.data)] + abs(centre.y)
  # print(plot.data.offset)
  # (b) convert into radial coords
  group <- NULL
  group$path <- ggradar:::CalculateGroupPath(plot.data.offset)
  
  # print(group$path)
  # (c) Calculate coordinates required to plot radial variable axes
  axis <- NULL
  axis$path <- ggradar:::CalculateAxisPath(var.names, grid.min + abs(centre.y), grid.max + abs(centre.y))
  # print(axis$path)
  # (d) Create file containing axis labels + associated plotting coordinates
  # Labels
  axis$label <- data.frame(
    text = axis.labels,
    x = NA,
    y = NA
  )
  # print(axis$label)
  # axis label coordinates
  n.vars <- length(var.names)
  angles <- seq(from = 0, to = 2 * pi, by = (2 * pi) / n.vars)
  axis$label$x <- sapply(1:n.vars, function(i, x) {
    ((grid.max + abs(centre.y)) * axis.label.offset) * sin(angles[i])
  })
  axis$label$y <- sapply(1:n.vars, function(i, x) {
    ((grid.max + abs(centre.y)) * axis.label.offset) * cos(angles[i])
  })
  # print(axis$label)
  # (e) Create Circular grid-lines + labels
  # caclulate the cooridinates required to plot circular grid-lines for three user-specified
  # y-axis values: min, mid and max [grid.min; grid.mid; grid.max]
  gridline <- NULL
  gridline$min$path <- ggradar:::funcCircleCoords(c(0, 0), grid.min + abs(centre.y), npoints = 360)
  gridline$mid$path <- ggradar:::funcCircleCoords(c(0, 0), grid.mid + abs(centre.y), npoints = 360)
  gridline$add$path <- ggradar:::funcCircleCoords(c(0, 0), grid.add + abs(centre.y), npoints = 360)
  gridline$max$path <- ggradar:::funcCircleCoords(c(0, 0), grid.max + abs(centre.y), npoints = 360)
  # print(head(gridline$max$path))
  # gridline labels
  gridline$min$label <- data.frame(
    x = gridline.label.offset, y = grid.min + abs(centre.y),
    text = as.character(grid.min)
  )
  gridline$max$label <- data.frame(
    x = gridline.label.offset, y = grid.max + abs(centre.y),
    text = as.character(grid.max)
  )
  gridline$mid$label <- data.frame(
    x = gridline.label.offset, y = grid.mid + abs(centre.y),
    text = as.character(grid.mid)
  )
  gridline$add$label <- data.frame(
    x = gridline.label.offset, y = grid.add + abs(centre.y),
    text = as.character(grid.add)
  )
  # print(gridline$min$label)
  # print(gridline$max$label)
  # print(gridline$mid$label)
  ### Start building up the radar plot
  
  # Declare 'theme_clear', with or without a plot legend as required by user
  # [default = no legend if only 1 group [path] being plotted]
  theme_clear <- theme_bw(base_size = base.size) +
    theme(
      axis.text.y = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      legend.key = element_rect(linetype = "blank")
    )
  
  if (plot.legend == FALSE) legend.position <- "none"
  
  # Base-layer = axis labels + plot extent
  # [need to declare plot extent as well, since the axis labels don't always
  # fit within the plot area automatically calculated by ggplot, even if all
  # included in first plot; and in any case the strategy followed here is to first
  # plot right-justified labels for axis labels to left of Y axis for x< (-x.centre.range)],
  # then centred labels for axis labels almost immediately above/below x= 0
  # [abs(x) < x.centre.range]; then left-justified axis labels to right of Y axis [x>0].
  # This building up the plot in layers doesn't allow ggplot to correctly
  # identify plot extent when plotting first (base) layer]
  
  # base layer = axis labels for axes to left of central y-axis [x< -(x.centre.range)]
  base <- ggplot(axis$label) +
    xlab(NULL) +
    ylab(NULL) +
    coord_equal() +
    geom_text(
      data = subset(axis$label, axis$label$x < (-x.centre.range)),
      aes(x = x, y = y, label = text), size = axis.label.size, hjust = 1, family = font.radar
    ) +
    scale_x_continuous(limits = c(-1.5 * plot.extent.x, 1.5 * plot.extent.x)) +
    scale_y_continuous(limits = c(-plot.extent.y, plot.extent.y))
  
  # ... + circular grid-lines at 'min', 'mid' and 'max' y-axis values
  base <- base + geom_path(
    data = gridline$min$path, aes(x = x, y = y),
    lty = gridline.min.linetype, colour = gridline.min.colour, linewidth = grid.line.width
  )
  base <- base + geom_path(
    data = gridline$mid$path, aes(x = x, y = y),
    lty = gridline.mid.linetype, colour = gridline.mid.colour, linewidth = grid.line.width
  )
  base <- base + geom_path(
    data = gridline$add$path, aes(x = x, y = y),
    lty = gridline.add.linetype, colour = gridline.add.colour, linewidth = grid.line.width
  )
  base <- base + geom_path(
    data = gridline$max$path, aes(x = x, y = y),
    lty = gridline.max.linetype, colour = gridline.max.colour, linewidth = grid.line.width
  )
  
  # + axis labels for any vertical axes [abs(x)<=x.centre.range]
  base <- base + geom_text(
    data = subset(axis$label, abs(axis$label$x) <= x.centre.range),
    aes(x = x, y = y, label = text), size = axis.label.size, hjust = 0.5, family = font.radar
  )
  # + axis labels for any vertical axes [x>x.centre.range]
  base <- base + geom_text(
    data = subset(axis$label, axis$label$x > x.centre.range),
    aes(x = x, y = y, label = text), size = axis.label.size, hjust = 0, family = font.radar
  )
  # + theme_clear [to remove grey plot background, grid lines, axis tick marks and axis text]
  base <- base + theme_clear
  #  + background circle against which to plot radar data
  base <- base + geom_polygon(
    data = gridline$max$path, aes(x, y),
    fill = background.circle.colour,
    alpha = background.circle.transparency
  )
  
  # + radial axes
  base <- base + geom_path(
    data = axis$path, aes(x = x, y = y, group = axis.no),
    colour = axis.line.colour
  )
  
  theGroupName <- names(group$path[1])
  
  # ... + group (cluster) 'paths'
  # base <- base + geom_path(
  #   data = group$path, aes(x = .data[["x"]], y = .data[["y"]], group = theGroupName, colour = theGroupName),
  #   size = group.line.width
  # )
  if (length(line.alpha) == 1) {
    base <- base + geom_path(data = group$path, aes(x = .data[["x"]], y = .data[["y"]], group = .data[[theGroupName]], colour = .data[[theGroupName]]), linewidth = group.line.width, alpha = line.alpha)
  } else {
    # Assuming line.alpha is a vector with the same length as the number of groups
    # This will apply different alpha values to each line
    base <- base + geom_path(data = group$path, aes(x = .data[["x"]], y = .data[["y"]], group = .data[[theGroupName]], colour = .data[[theGroupName]]), linewidth = group.line.width) +
      scale_alpha_manual(values = line.alpha)
  }
  
  # ... + group points (cluster data)
  # Modify point drawing logic based on draw.points
  if (draw.points) {
    # Check if point.alpha is a vector or single value
    if (length(point.alpha) == 1) {
      base <- base + geom_point(data = group$path, aes(x = .data[["x"]], y = .data[["y"]], group = .data[[theGroupName]], colour = .data[[theGroupName]]), size = group.point.size, alpha = point.alpha)
    } else {
      # Assuming point.alpha is a vector with the same length as the number of groups
      # This will apply different alpha values to each group
      base <- base + geom_point(data = group$path, aes(x = .data[["x"]], y = .data[["y"]], group = .data[[theGroupName]], colour = .data[[theGroupName]]), size = group.point.size) +
        scale_alpha_manual(values = point.alpha)
    }
  }
  
  # ... + group (cluster) fills
  if (fill == TRUE) {
    base <- base + geom_polygon(data = group$path, aes(x = .data[["x"]], y = .data[["y"]], group = .data[[theGroupName]], fill = .data[[theGroupName]]), alpha = fill.alpha)
  }
  
  
  # ... + amend Legend title
  if (plot.legend == TRUE) base <- base + labs(colour = legend.title, size = legend.text.size)
  
  # ... + grid-line labels (max; mid; min)
  if (label.gridline.min == TRUE) {
    base <- base + geom_text(aes(x = x, y = y, label = values.radar[1]), data = gridline$min$label, size = grid.label.size * 0.8, hjust = 1, family = font.radar)
  }
  if (label.gridline.mid == TRUE) {
    base <- base + geom_text(aes(x = x, y = y, label = values.radar[2]), data = gridline$mid$label, size = grid.label.size * 0.8, hjust = 1, family = font.radar)
  }
  if (label.gridline.add == TRUE) {
    base <- base + geom_text(aes(x = x, y = y, label = values.radar[3]), data = gridline$add$label, size = grid.label.size * 0.8, hjust = 1, family = font.radar)
  }
  if (label.gridline.max == TRUE) {
    base <- base + geom_text(aes(x = x, y = y, label = values.radar[4]), data = gridline$max$label, size = grid.label.size * 0.8, hjust = 1, family = font.radar)
  }
  # ... + centre.y label if required [i.e. value of y at centre of plot circle]
  if (label.centre.y == TRUE) {
    centre.y.label <- data.frame(x = 0, y = 0, text = as.character(centre.y))
    base <- base + geom_text(aes(x = x, y = y, label = text), data = centre.y.label, size = grid.label.size, hjust = 0.5, family = font.radar)
  }
  
  if (!is.null(group.colours)) {
    colour_values <- rep(group.colours, length(unique(plot.data[, 1])) / length(group.colours))
  } else {
    colour_values <- generate_color_values(length(unique(plot.data[, 1])))
  }
  
  base <- base +
    theme(
      legend.key.width = unit(3, "line"),
      text = element_text(
        size = 20,
        family = font.radar
      )
    ) +
    theme(legend.text = element_text(size = legend.text.size), legend.position = legend.position) +
    theme(legend.key.height = unit(2, "line")) +
    scale_colour_manual(values = colour_values) +
    theme(text = element_text(family = font.radar)) +
    theme(legend.title = element_blank())
  
  
  if (isTRUE(fill)) {
    base <- base +
      scale_fill_manual(values = colour_values, guide = "none")
  }
  
  if (legend.title != "") {
    base <- base + theme(legend.title = element_text())
  }
  
  if (plot.title != "") {
    base <- base + ggtitle(plot.title)
  }
  
  return(base)
}

