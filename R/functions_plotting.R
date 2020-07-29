# Plotting style 
PlotMystyle = function(p, title=NULL, col=NULL, fill=NULL, legend_title=NULL, legend_position=NULL) {
  p = p + theme_light() + theme(panel.border = element_blank())
  if (!is.null(title)) p = p + ggtitle(title) #+ theme(plot.title = element_text(hjust=0.5))
  if (length(col) > 0) p = p + scale_colour_manual(values=col)
  if (length(fill) > 0) p = p + scale_fill_manual(values=fill)
  if (!is.null(legend_title)) {
    p = p + labs(color=legend_title, fill=legend_title)
  } else {
    p = p + theme(legend.title = element_blank()) 
  }
  if (!is.null(legend_position)) p = p + theme(legend.position=legend_position)
  return(p)
}

# Transform a matrix cells (rows) x htos (cols) into a format that can be understood by 
#   feature_grid: cell class, name hto1, value hto1, name hto2, value hto2
DfAllColumnCombinations = function(x, cell_classification) {
  out = matrix(NA, nrow=0, ncol=4)
  for (i in 1:(ncol(x)-1)) {
    for (j in (i+1):ncol(x)) {
      a = x %>% dplyr::select(i) %>% tidyr::pivot_longer(1)
      b = x %>% dplyr::select(j) %>% tidyr::pivot_longer(1)
      out = rbind(out, cbind(cell_classification, a, b))
    }
  }
  colnames(out) = c("cell_classification", "name1", "value1", "name2", "value2")
  
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
    geom_line(aes(x=cell, y=med))
  p = p + geom_point(data=y_outlier, aes(x=cell, y=out, colour=id), size=0.5) + scale_color_manual(values=col)
  p = PlotMystyle(p) + 
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank()) + 
    xlab("Cells") + ylab("Relative log expression")
  p
  
  return(p)
}