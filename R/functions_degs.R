#' Convert natural log of Seurat into log2
#' 
#' @param degs Result table of the "Seurat::FindAllMarkers" and "Seurat::FindMarkers" functions
DegsConvertLog = function(degs) { 
  
  # Convert natural log to log2
  # https://github.com/satijalab/seurat/issues/3346#issue-672668498
  # log_a(b) * log_c(a) = log_c(b)
  #   log_n(fc) = log_2(fc) * log_n(2)
  #   log_2(fc) = log_n(fc) / log_n(2)
  degs$avg_lognfc = degs$avg_logFC
  degs$avg_log2fc = degs$avg_lognfc / log(2)
  
  return(degs)
}

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
    dplyr::filter((avg_log2FC <= -cut_log2FC) | (avg_log2FC >= cut_log2FC)) %>% 
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
    dplyr::select(cluster, gene, avg_log2FC, p_val, p_val_adj, pct.1, pct.2) %>% 
    as.data.frame()
  
  # Format columns of top 5 markers
  top$p_val = formatC(as.numeric(top$p_val), format="e", digits=1)
  top$p_val_adj = formatC(as.numeric(top$p_val_adj), format="e", digits=1)
  top$avg_log2FC = round(top$avg_log2FC, digits=3)
  
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
  avg_set[[DefaultAssay(sc)]] = c("counts", "data", "scale.data")
  avg_data = matrix(NA+0, nrow=length(genes), ncol=0)
  
  identities = levels(Idents(sc))
  for(as in names(avg_set)) { 
    for(sl in avg_set[[as]]) { 
      avg_per_id = mapply(function(id) { 
        id_cells = WhichCells(sc, idents=id)
        id_avg = rowMeans(GetAssayData(sc, assay=as, slot=sl)[genes, id_cells])
        return(id_avg)
      }, identities)
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

#' Plot the number of DEGs per test 
#' 
#' @param degs_up Result table of the "Seurat::FindAllMarkers" function (requires column "cluster"), filtered for up-regulated genes
#' @param degs_down Result table of the "Seurat::FindAllMarkers" function (requires column "cluster"), filtered for down-regulated genes
#' @param title Plot title
DegsPlotNumbers = function(degs_up, degs_down, title=NULL) { 
  
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