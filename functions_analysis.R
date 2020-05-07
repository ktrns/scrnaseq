prelim.analysis <- function(sc, mt="^MT-", pc.n=10) {
  
  # Remove cells with low counts, or high mitochondrial content
  # I have not gotten it to work to define the cutoffs as variables and get the subset 
  #   function to work - absurd it is
  sc = PercentageFeatureSet(sc, pattern=mt, col.name="percent.mt")
  sc = subset(sc, subset=nFeature_RNA>200 & percent.mt<20)
  sc = subset(sc, subset=nFeature_RNA>200 & percent.mt<20)
  
  # Normalize
  sc = SCTransform(sc)
  
  # Run PCA
  sc = RunPCA(sc, features=VariableFeatures(object=sc))
  
  # Run UMAP 
  sc = RunUMAP(sc, dims=1:pc.n)
  
  return(sc)
}

