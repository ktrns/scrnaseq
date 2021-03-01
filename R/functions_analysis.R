#' Calculate enrichment of cells per sample per cluster.
#' 
#' @param sc Seurat object.
#' @return A table with counts, odd ratios and p-values.
cells_fisher = function(sc) {
  cell_samples = sc[[]] %>% dplyr::pull(orig.ident) %>% unique() %>% sort()
  cell_clusters = sc[[]] %>% dplyr::pull(seurat_clusters) %>% unique() %>% sort()
  out = matrix(0+NA, nrow=0, ncol=length(cell_clusters)) %>% as.data.frame()
  for(s in cell_samples) {
    ft.list = lapply(cell_clusters, function(cl) { 
      a = sc[[]] %>% dplyr::filter(orig.ident==s, seurat_clusters==cl) %>% dplyr::count() %>% as.numeric()
      b = sc[[]] %>% dplyr::filter(orig.ident!=s, seurat_clusters==cl) %>% dplyr::count() %>% as.numeric()
      c = sc[[]] %>% dplyr::filter(orig.ident==s, seurat_clusters!=cl) %>% dplyr::count() %>% as.numeric()
      d = sc[[]] %>% dplyr::filter(orig.ident!=s, seurat_clusters!=cl) %>% dplyr::count() %>% as.numeric()
      tbl.2by2 = matrix(c(a, b, c, d), ncol=2, nrow=2, byrow=TRUE)
      ft = fisher.test(tbl.2by2, alternative="greater")
      return(c(oddsRatio=round(as.numeric(ft$estimate), 2),
               p=formatC(as.numeric(ft$p.value), format="e", digits=1)))
    })
    ft.matrix = purrr::reduce(ft.list, .f=cbind)
    rownames(ft.matrix) = paste0(s, ".", rownames(ft.matrix))
    out = rbind(out, ft.matrix)
  }
  colnames(out) = paste0("Cl_", cell_clusters)
  return(out)
}

#' Calculate cell cycle scores
#' 
#' @param sc Seurat object
#' @param genes_s Vector of gene names characteristic for S-phase
#' @param genes_g2m Vector of gene names characteristic for G2M-phase
#' @param name Name of the dataset, only used to write a sensible warning if necessary (Default "")
#' @return Updated Seurat object
cc_scoring = function(sc, genes_s, genes_g2m, name=""){
  genes_s_exists = genes_s %in% rownames(sc)
  genes_g2m_exists = genes_g2m %in% rownames(sc)
  if (sum(genes_s_exists) >= 10 & sum(genes_s_exists) >= 10){
    sc = Seurat::CellCycleScoring(sc, 
                                  s.features=genes_s[genes_s_exists],
                                  g2m.features=genes_g2m[genes_g2m_exists], 
                                  set.ident=FALSE, verbose=FALSE)
    sc[["CC.Difference"]] = sc[["S.Score", drop=TRUE]] - sc[["G2M.Score", drop=TRUE]]
    sc[["Phase"]] = factor(sc[["Phase", drop=TRUE]], levels=c("G1", "G2M", "S"))
    
  } else {
    sc[["S.Score"]] = NA
    sc[["G2M.Score"]] = NA
    sc[["Phase"]] = factor(NA, levels=c("G1", "G2M", "S"))
    if(name!="") name=paste0(name, " ")
    warning(paste0("There are not enough G2/M and S phase markers in the dataset ", name, "to reliably determine cell cycle scores and phases. Scores and phases will be set to NA and removal of cell cycle effects is skipped."))
  }
  
  return(sc)
}