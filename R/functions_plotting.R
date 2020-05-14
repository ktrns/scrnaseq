# Plotting style 
PlotMystyle = function(p, title=NULL, col=NULL, legend_title=NULL, legend_position=NULL) {
  p = p + theme_light() + theme(panel.border = element_blank())
  if (!is.null(title)) p = p + ggtitle(title) #+ theme(plot.title = element_text(hjust=0.5))
  if (length(col) > 0) p = p + scale_fill_manual(values=col)
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
