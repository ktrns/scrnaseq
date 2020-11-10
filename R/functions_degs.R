#' Introduce signed p-value and sort table of differentially expressed genes per performed test
#' 
#' @param degs Result table of the "Seurat::FindAllMarkers" function (requires column "cluster") 
DegsSort = function(degs) { 
  
  # Introduce signed p-value score and sort table
  degs$p_val_adj_score = -log2(degs$p_val_adj) * sign(degs$avg_log2FC)
  degs = degs %>% 
    dplyr::arrange(cluster, -p_val_adj_score, -avg_log2FC) %>% 
    as.data.frame()
  
  return(degs)
}

#' Filter table of differentially expressed genes
#' 
#' @param degs Result table of the "Seurat::FindAllMarkers" and "Seurat::FindMarkers" functions, further processed with the "DegsConvertLog" function. 
DegsFilter = function(degs, cut_log2FC, cut_padj) { 
  
  # Filter differentially expressed genes based on p-value and fold-change 
  filt = degs %>% 
    dplyr::filter(p_val_adj <= cut_padj) %>% 
    dplyr::filter(abs(avg_log2FC) >= cut_log2FC) %>% 
    as.data.frame()
  
  # Separate up- and down-regulated genes
  down = filt %>% 
    dplyr::filter(avg_log2FC <= -cut_log2FC) %>% 
    as.data.frame()
  up = filt %>% 
    dplyr::filter(avg_log2FC >= cut_log2FC) %>% 
    as.data.frame()  
  
  # Return a list of three tables
  return(list(all=filt, 
              up=up,
              down=down))
}

#' Display top marker genes (=up-regulated genes)
#' 
#' @param degs Result table of the "Seurat::FindAllMarkers" function (requires column "cluster"), further processed with the "DegsConvertLog" function
DegsUpDisplayTop = function(degs, n=5, caption=NULL) { 
  
  # Get top 5 up-regulated markers
  top = degs %>% 
    dplyr::group_by(cluster) %>% 
    dplyr::top_n(n=n, wt=p_val_adj_score) %>% 
    dplyr::ungroup() %>% 
    dplyr::transmute(cluster=cluster,
                     gene=gene,
                     avg_log2FC=round(avg_log2FC, digits=3),
                     p_val=formatC(as.numeric(p_val), format="e", digits=1),
                     p_val_adj=formatC(as.numeric(p_val_adj), format="e", digits=1),
                     pct.1=pct.1,
                     pct.2=pct.2) %>% 
    as.data.frame()
  return(top)
}

#' Add average data per identies
#' 
#' @param sc Seurat object
#' @param genes Gene list for which average data are to be extracted
DegsAvgData = function(sc, genes) { 
  # The standard average log FC is derived from assay="RNA" and slot="data"
  # Add average scaled data per cluster for default assay
  avg_set = list()
  avg_set[["RNA"]] = "counts"
  avg_set[[DefaultAssay(sc)]] = c("data", "scale.data")
  avg_data = matrix(NA+0, nrow=length(genes), ncol=0)
  
  identities = levels(Idents(sc))
  for (as in names(avg_set)) { 
    for (sl in avg_set[[as]]) { 
      if (length(genes) > 0) {
        avg_per_id = mapply(function(id) { 
          id_cells = WhichCells(sc, idents=id)
          id_avg = rowMeans(GetAssayData(sc, assay=as, slot=sl)[genes, id_cells])
          return(id_avg)
        }, identities)
      } else {
        avg_per_id = matrix(NA, nrow=0, ncol=length(identities)) %>% as.data.frame()
      }
      colnames(avg_per_id) = paste0("avg_", as, "_", sl, "_id", identities)
      avg_data = cbind(avg_data, avg_per_id)
    }
  }
  return(avg_data)
}

#' Write differentially expressed genes to file
#' 
#' @param degs Result table of the "Seurat::FindAllMarkers" function (requires column "cluster"). The table is split by the column "cluster", and each split table is written into a separate Excel tab. 
#' @param annot_ensembl Ensembl annotation for all genes
#' @param file Output file name 
DegsWriteToFile = function(degs, annot_ensembl, file) {
  
  # Write results to file if there are degs
  if (nrow(degs) > 0) {
    
    # Convert to list, one table per cluster
    degs_cluster = levels(degs$cluster)
    degs_lst = lapply(degs_cluster, function(x) {degs %>% dplyr::filter(cluster==x)})
    names(degs_lst) = degs_cluster
  
    # Add Ensembl annotation
    for (i in seq(degs_lst)) {
      degs_ensembl = seurat_rowname_to_ensembl[degs_lst[[i]]$gene]
      degs_lst[[i]] = cbind(degs_lst[[i]], annot_ensembl[degs_ensembl,])
    }
    
    # Output in Excel sheet
    openxlsx::write.xlsx(degs_lst, file=file)
  } 
}

#' Plot the number of DEGs per test 
#' 
#' @param degs_up Result table of the "Seurat::FindAllMarkers" function (requires column "cluster"), filtered for up-regulated genes
#' @param degs_down Result table of the "Seurat::FindAllMarkers" function (requires column "cluster"), filtered for down-regulated genes
#' @param title Plot title
DegsPlotNumbers = function(degs_up, degs_down, title=NULL) { 
  
  if ((nrow(degs_up) > 0) | (nrow(degs_down) > 0)) {
    degs_cluster = sort(unique(c(levels(degs_up$cluster), levels(degs_down$cluster))))
    degs_n = cbind(Down=sapply(degs_cluster, function(x) sum(degs_up$cluster == x)), 
                   Up=sapply(degs_cluster, function(x) sum(degs_down$cluster == x))) %>% 
      as.data.frame()
  
    degs_n = cbind(Identity=factor(degs_cluster, levels=degs_cluster), degs_n) %>% 
      tidyr::pivot_longer(cols=c("Down", "Up"), 
                          names_to="Direction", 
                          values_to="n")
    
    p = ggplot(degs_n, aes(x=Identity, y=n, fill=Direction)) + 
      geom_bar(stat="identity") + 
      AddStyle(title=title,
               fill=c("steelblue", "darkgoldenrod1"))
    return(p)
  }
}