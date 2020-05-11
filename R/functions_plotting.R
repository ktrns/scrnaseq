# Plotting style 
plot.mystyle = function(p, title=NULL, col=NULL, legend.title=NULL, legend.position=NULL) {
  p = p + theme_light() + theme(panel.border = element_blank())
  if(!is.null(title)) p = p + ggtitle(title) #+ theme(plot.title = element_text(hjust=0.5))
  if(length(col) > 0) p = p + scale_fill_manual(values=col)
  if(!is.null(legend.title)) {
    p = p + labs(color=legend.title, fill=legend.title)
  } else {
    p = p + theme(legend.title = element_blank()) 
  }
  if(!is.null(legend.position)) p = p + theme(legend.position=legend.position)
  return(p)
}

# Transform a matrix cells (rows) x htos (cols) into a format that can be understood by 
#   feature_grid: cell class, name hto1, value hto1, name hto2, value hto2
df.all.col.combinations = function(x, cell.classification) {
  out = matrix(NA, nrow=0, ncol=4)
  for(i in 1:(ncol(x)-1)) {
    for(j in (i+1):ncol(x)){
      a = x %>% select(i) %>% pivot_longer(1)
      b = x %>% select(j) %>% pivot_longer(1)
      out=rbind(out, cbind(cell.classification, a, b))
    }
  }
  colnames(out) = c("cell.classification", "name1", "value1", "name2", "value2")
  
  # Define plot order so that the two levels of interest are always on top, then negatives, doublets, 
  #   and finally all other samples
  out$order = 0
  out[out$cell.classification==out$name1, "order"] = 3
  out[out$cell.classification==out$name2, "order"] = 3
  out[out$cell.classification=="Negative", "order"] = 2
  out[out$cell.classification=="Doublet", "order"] = 1
  out = out %>% group_by(name1, name2) %>% arrange(order)
  out$order = NULL # remove column again -> only needed to order data points
  
  return(out)
}
