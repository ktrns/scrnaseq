# Plotting style
AddStyle = function(title=NULL, col=NULL, fill=NULL, legend_title=NULL, legend_position=NULL, xlab=NULL, ylab=NULL) {
  list(
    theme_light() + theme(panel.border = element_blank()), 
    if (!is.null(title)) ggtitle(title), 
    if (length(col) > 0) scale_colour_manual(values=col),
    if (length(fill) > 0)  scale_fill_manual(values=fill),
    if (!is.null(legend_title)) {
      labs(color=legend_title, fill=legend_title)
    } else {
      theme(legend.title = element_blank()) 
    },
    if (!is.null(legend_position)) theme(legend.position=legend_position),
    if (!is.null(xlab)) xlab(xlab),
    if (!is.null(ylab)) ylab(ylab)
  )
}
# Transform a matrix cells (rows) x htos (cols) into a format that can be understood by 
#   feature_grid: cell class, name hto1, value hto1, name hto2, value hto2
DfAllColumnCombinations = function(x, cell_classification) {
  out = combn(x, 2, simplify=FALSE)
  out = lapply(out, function(o) {
    return(data.frame(cell_classification=unname(cell_classification[rownames(o)]), name1=colnames(o)[1], value1=o[,1], name2=colnames(o)[2], value2=o[,2]))
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