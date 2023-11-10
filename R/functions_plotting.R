#' Out plotting style.
#'
#'@param title The plot title.
#'@param col A vector of colours to use.
#'@param fill A vector of fill colours to use.
#'@param legend_title The legend title.
#'@param legend_position The legend position.
#'@param xlab The title of the x-axis.
#'@param ylab The title of the y-axis.
#'@return None, add as theme.
AddStyle = function(title=NULL, col=NULL, fill=NULL, legend_title=NULL, legend_position=NULL, xlab=NULL, ylab=NULL) {
  list(
    ggplot2::theme_light(11) + ggplot2::theme(panel.border = element_blank()), 
    if (!is.null(title)) ggtitle(title), 
    if (length(col) > 0) ggplot2::scale_colour_manual(values=col),
    if (length(fill) > 0)  ggplot2::scale_fill_manual(values=fill),
    if (!is.null(legend_title)) {
      labs(color=legend_title, fill=legend_title)
    } else {
      theme(legend.title = element_blank()) 
    },
    if (!is.null(legend_position)) ggplot2::theme(legend.position=legend_position),
    if (!is.null(xlab)) xlab(xlab),
    if (!is.null(ylab)) ylab(ylab)
  )
}
#' Transform a matrix cells (rows) x htos (cols) into a format that can be understood by feature_grid: cell class, name hto1, value hto1, name hto2, value hto2
#' 
#' @param x: A matrix cells (rows) x htos (cols).
#' @param cell_classification A vector of cell classifications.
DfAllColumnCombinations = function(x, cell_classification) {
  out = combn(x, 2, simplify=FALSE)
  out = lapply(out, function(o) {
    return(data.frame(cell_classification=unname(cell_classification[rownames(o)]), name1=colnames(o)[1], value1=o[, 1], name2=colnames(o)[2], value2=o[, 2]))
  }) %>% dplyr::bind_rows()
  
  # Define plot order so that the two levels of interest are always on top, then negatives, doublets, 
  #   and finally all other samples
  out$order = 0
  out[out$cell_classification==out$name1, "order"] = 3
  out[out$cell_classification==out$name2, "order"] = 3
  out[out$cell_classification=="Negative", "order"] = 2
  out[out$cell_classification=="Doublet", "order"] = 1
  out = out %>% dplyr::group_by(name1, name2) %>% dplyr::arrange(order)
  out$order = NULL # remove column again -> only needed to order data points
  
  return(out)
}

# Plot Relative log expression per cell 
PlotRLE = function(x, col, id) { 
  # x - data
  # id - cell identity (same order as cells)
  # col - colours for cell identities

  # Median of a gene across all cells
  genes.median = sapply(1:nrow(x), function(gene) median(x[gene,], na.rm=TRUE))
  
  # Subtract gene median from gene count
  y = x - genes.median
  
  # Get statistics
  y_stats = boxplot(y, plot=FALSE)
  y_outlier_x = y_stats$group
  y_outlier_y = y_stats$out
  y_outlier = cbind(cell=y_outlier_x, out=y_outlier_y) %>% as.data.frame()
  y_outlier$id = id[y_outlier$cell]
  y_stats = y_stats$stats %>% t() %>% as.data.frame()
  colnames(y_stats) = c("lowerWhisker", "q25", "med", "q75", "upperWhisker")
  y_stats$cell = 1:nrow(y_stats) # convert x-axis to numeric
  
  # Actual plotting
  p = ggplot(y_stats) +
    geom_ribbon(aes(x=cell, ymin=q75, ymax=upperWhisker), fill="darkgrey") + 
    geom_ribbon(aes(x=cell, ymin=med, ymax=q75), fill="lightgrey") + 
    geom_ribbon(aes(x=cell, ymin=q25, ymax=med), fill="lightgrey") + 
    geom_ribbon(aes(x=cell, ymin=lowerWhisker, ymax=q25), fill="darkgrey") + 
    geom_line(aes(x=cell, y=med)) + 
    geom_point(data=y_outlier, aes(x=cell, y=out, colour=id), size=0.5) + 
    scale_color_manual(values=col) + 
    AddStyle(xlab="Cells", ylab="Relative log expression") + 
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
  p
  
  return(p)
}

# Re-define Dotplot for non-scaled dotplots
# Keep this code in here until Seurat fixes this
# Issue: https://github.com/satijalab/seurat/issues/4298
PercentAbove <- function(x, threshold) {
  return(length(x = x[x > threshold]) / length(x = x))
}

DotPlotUpdated = function (object, assay = NULL, features, cols = c("lightgrey", 
                                                                    "blue"), col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6, 
                           idents = NULL, group.by = NULL, split.by = NULL, cluster.idents = FALSE, 
                           scale = TRUE, scale.by = "radius", scale.min = NA, scale.max = NA) 
{
  # assay <- assay %||% DefaultAssay(object = object)
  assay <- DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  split.colors <- !is.null(x = split.by) && !any(cols %in% 
                                                   rownames(x = brewer.pal.info))
  scale.func <- switch(EXPR = scale.by, size = scale_size, 
                       radius = scale_radius, stop("DotPlotUpdated: 'scale.by' must be either 'size' or 'radius'"))
  feature.groups <- NULL
  if (is.list(features) | any(!is.na(names(features)))) {
    feature.groups <- unlist(x = sapply(X = 1:length(features), 
                                        FUN = function(x) {
                                          return(rep(x = names(x = features)[x], each = length(features[[x]])))
                                        }))
    if (any(is.na(x = feature.groups))) {
      warning("Some feature groups are unnamed.", call. = FALSE, 
              immediate. = TRUE)
    }
    features <- unlist(x = features)
    names(x = feature.groups) <- features
  }
  cells <- unlist(x = CellsByIdentities(object = object, idents = idents))
  data.features <- FetchData(object = object, vars = features, 
                             cells = cells)
  data.features$id <- if (is.null(x = group.by)) {
    Idents(object = object)[cells, drop = TRUE]
  }
  else {
    object[[group.by, drop = TRUE]][cells, drop = TRUE]
  }
  if (!is.factor(x = data.features$id)) {
    data.features$id <- factor(x = data.features$id)
  }
  id.levels <- levels(x = data.features$id)
  data.features$id <- as.vector(x = data.features$id)
  if (!is.null(x = split.by)) {
    splits <- object[[split.by, drop = TRUE]][cells, drop = TRUE]
    if (split.colors) {
      if (length(x = unique(x = splits)) > length(x = cols)) {
        stop("DotPlotUpdated: Not enough colors for the number of groups")
      }
      cols <- cols[1:length(x = unique(x = splits))]
      names(x = cols) <- unique(x = splits)
    }
    data.features$id <- paste(data.features$id, splits, sep = "_")
    unique.splits <- unique(x = splits)
    id.levels <- paste0(rep(x = id.levels, each = length(x = unique.splits)), 
                        "_", rep(x = unique(x = splits), times = length(x = id.levels)))
  }
  data.plot <- lapply(X = unique(x = data.features$id), FUN = function(ident) {
    data.use <- data.features[data.features$id == ident, 
                              1:(ncol(x = data.features) - 1), drop = FALSE]
    avg.exp <- apply(X = data.use, MARGIN = 2, FUN = function(x) {
      return(mean(x = expm1(x = x)))
    })
    pct.exp <- apply(X = data.use, MARGIN = 2, FUN = PercentAbove, 
                     threshold = 0)
    return(list(avg.exp = avg.exp, pct.exp = pct.exp))
  })
  names(x = data.plot) <- unique(x = data.features$id)
  if (cluster.idents) {
    mat <- do.call(what = rbind, args = lapply(X = data.plot, 
                                               FUN = unlist))
    mat <- scale(x = mat)
    id.levels <- id.levels[hclust(d = dist(x = mat))$order]
  }
  data.plot <- lapply(X = names(x = data.plot), FUN = function(x) {
    data.use <- as.data.frame(x = data.plot[[x]])
    data.use$features.plot <- rownames(x = data.use)
    data.use$id <- x
    return(data.use)
  })
  data.plot <- do.call(what = "rbind", args = data.plot)
  if (!is.null(x = id.levels)) {
    data.plot$id <- factor(x = data.plot$id, levels = id.levels)
  }
  if (length(x = levels(x = data.plot$id)) == 1) {
    scale <- FALSE
    warning("Only one identity present, the expression values will be not scaled", 
            call. = FALSE, immediate. = TRUE)
  }
  avg.exp.scaled <- sapply(X = unique(x = data.plot$features.plot), 
                           FUN = function(x) {
                             data.use <- data.plot[data.plot$features.plot == 
                                                     x, "avg.exp"]
                             if (scale) {
                               data.use <- scale(x = data.use)
                               data.use <- MinMax(data = data.use, min = col.min, 
                                                  max = col.max)
                             }
                             else {
                               data.use <- log1p(x = data.use)
                             }
                             return(data.use)
                           })
  avg.exp.scaled <- as.vector(x = t(x = avg.exp.scaled))
  if (split.colors) {
    avg.exp.scaled <- as.numeric(x = cut(x = avg.exp.scaled, 
                                         breaks = 20))
  }
  data.plot$avg.exp.scaled <- avg.exp.scaled
  data.plot$features.plot <- factor(x = data.plot$features.plot, 
                                    levels = features)
  data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
  data.plot$pct.exp <- data.plot$pct.exp * 100
  if (split.colors) {
    splits.use <- vapply(X = as.character(x = data.plot$id), 
                         FUN = gsub, FUN.VALUE = character(length = 1L), pattern = paste0("^((", 
                                                                                          paste(sort(x = levels(x = object), decreasing = TRUE), 
                                                                                                collapse = "|"), ")_)"), replacement = "", 
                         USE.NAMES = FALSE)
    data.plot$colors <- mapply(FUN = function(color, value) {
      return(colorRampPalette(colors = c("grey", color))(20)[value])
    }, color = cols[splits.use], value = avg.exp.scaled)
  }
  color.by <- ifelse(test = split.colors, yes = "colors", no = "avg.exp.scaled")
  if (!is.na(x = scale.min)) {
    data.plot[data.plot$pct.exp < scale.min, "pct.exp"] <- scale.min
  }
  if (!is.na(x = scale.max)) {
    data.plot[data.plot$pct.exp > scale.max, "pct.exp"] <- scale.max
  }
  if (!is.null(x = feature.groups)) {
    data.plot$feature.groups <- factor(x = feature.groups[data.plot$features.plot], 
                                       levels = unique(x = feature.groups))
  }
  plot <- ggplot(data = data.plot, mapping = aes_string(x = "features.plot", 
                                                        y = "id")) + geom_point(mapping = aes_string(size = "pct.exp", 
                                                                                                     color = color.by)) + scale.func(range = c(0, dot.scale), 
                                                                                                                                     limits = c(scale.min, scale.max)) + theme(axis.title.x = element_blank(), 
                                                                                                                                                                               axis.title.y = element_blank()) + guides(size = guide_legend(title = "Percent Expressed")) + 
    labs(x = "Features", y = ifelse(test = is.null(x = split.by), 
                                    yes = "Identity", no = "Split Identity")) + cowplot::theme_cowplot()
  if (!is.null(x = feature.groups)) {
    plot <- plot + facet_grid(facets = ~feature.groups, scales = "free_x", 
                              space = "free_x", switch = "y") + theme(panel.spacing = unit(x = 1, 
                                                                                           units = "lines"), strip.background = element_blank())
  }
  if (split.colors) {
    plot <- plot + scale_color_identity()
  }
  else if (length(x = cols) == 1) {
    plot <- plot + scale_color_distiller(palette = cols)
  }
  else {
    plot <- plot + scale_color_gradient(low = cols[1], high = cols[2])
  }
  if (!split.colors) {
    plot <- plot + guides(color = guide_colorbar(title = "Average Expression"))
  }
  return(plot)
}
