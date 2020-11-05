prelim.analysis <- function(sc, mt="^MT-", pc_n=10) {
  
  # Remove cells with low counts, or high mitochondrial content
  # I have not gotten it to work to define the cutoffs as variables and get the subset 
  #   function to work - absurd it is
  sc = Seurat::PercentageFeatureSet(sc, pattern=mt, col.name="percent.mt")
  sc = subset(sc, subset=nFeature_RNA>200 & percent.mt<20)
  sc = subset(sc, subset=nFeature_RNA>200 & percent.mt<20)
  
  # Normalize
  sc = Seurat::SCTransform(sc)
  
  # Run PCA
  sc = Seurat::RunPCA(sc, features=VariableFeatures(object=sc))
  
  # Run UMAP 
  sc = Seurat::RunUMAP(sc, dims=1:pc_n)
  
  return(sc)
}

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