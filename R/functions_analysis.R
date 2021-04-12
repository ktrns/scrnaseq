#' Calculate enrichment of cells per sample per cluster.
#' 
#' @param sc Seurat object.
#' @return A table with counts, odd ratios and p-values.
CellsFisher = function(sc) {
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
CCScoring = function(sc, genes_s, genes_g2m, name=""){
  # Subset to CC genes that are in the object
  genes_s_exists = genes_s %in% rownames(sc)
  genes_g2m_exists = genes_g2m %in% rownames(sc)
  if (sum(genes_s_exists) >= 10 & sum(genes_s_exists) >= 10){
    
    # Calculate CC scores
    sc = Seurat::CellCycleScoring(sc, 
                                  s.features=genes_s[genes_s_exists],
                                  g2m.features=genes_g2m[genes_g2m_exists], 
                                  set.ident=FALSE, verbose=FALSE)
    sc[["CC.Difference"]] = sc[["S.Score", drop=TRUE]] - sc[["G2M.Score", drop=TRUE]]
    sc[["Phase"]] = factor(sc[["Phase", drop=TRUE]], levels=c("G1", "G2M", "S"))
    
    # Add to 'gene_lists' slot in the misc slot of the Seurat object
    sc = ScAddLists(sc, lists=list(CC_S_phase=genes_s[genes_s_exists], CC_G2M_phase=genes_g2m[genes_g2m_exists]), lists_slot="gene_lists")
    
  } else {
    sc[["S.Score"]] = sc[["G2M.Score"]] = sc[["CC.Difference"]] = NA
    sc[["Phase"]] = factor(NA, levels=c("G1", "G2M", "S"))
    if(name!="") name=paste0(name, " ")
    warning(paste0("There are not enough G2/M and S phase markers in the dataset ", name, "to reliably determine cell cycle scores and phases. Scores and phases will be set to NA and removal of cell cycle effects is skipped."))
  }
  
  return(sc)
}

#' Integrate multiple samples using Seurat's integration strategy. In short, a set of features (e.g. genes) is used to anchor cells that are in a matched biological state. Within these anchored set of cells, technical effects are then removed. There are three possible ways to run the integration process:
#' - Default: Anchors are computed for all pairs of datasets. This will give all datasets the same weight during dataset integration but can be computationally intensive.
#' - Provide a reference: One dataset is used as reference and anchors are computed for all other datasets. This approach is computational faster but less accurate.
#' - Use reciprocal PCA: Anchors are not computed based on features (e.g. genes) but in PCA space which reduces the complexity. This approach is much faster and recommended for large datasets (>50000 cells). However, it is also less accurate.
#' 
#' @param sc List of seurat objects with variable features
#' @param ndims Number of dimensions used for integration (Default: min(30, minimum number of cells in a sample - 1))
#' @param reference Use one or more datasets as reference. Separate multiple datasets by comma (Default: NULL)
#' @param use_reciprocal_pca Use reciprocal PCA for cell anchoring (Default: FALSE)
#' @param verbose Be verbose (Default: FALSE)
#' @param assay Can be 'RNA' or 'SCT' (Default: RNA)
#' @param k_filter How many neighbors to use when filtering anchors. If NULL, automatically set to min(200, minimum number of cells in a sample)).
#' @param k_weight Number of neighbors to consider when weighting anchors. If NULL, automatically set to min(100, minimum number of cells in a sample)).
#' @param k_anchor How many neighbors to use when picking anchors. If NULL, automatically set to min(5, minimum number of cells in a sample)).
#' @param vars_to_regress For reciprocal PCA: when doing the scaling, which variables should be regressed out when doing the scaling
#' @param min_cells For reciprocal PCA: when doing the scaling for SCT, the minimum number of cells a gene should be expressed
#' @return A Seurat object with an integrated assay (RNAintegrated or SCTintegrated) and a merged assay (RNA or SCT)
RunIntegration = function(sc, ndims=30, reference=NULL, use_reciprocal_pca=FALSE, verbose=FALSE, assay="RNA", k_filter=NULL, k_weight=NULL, k_anchor=NULL, vars_to_regress=NULL, min_cells=1) {
  # THIS FUNCTION NEEDS TO BE REVIEWED
  
  # Note: Assay names should only have numbers and letters
  # Warning: Keys should be one or more alphanumeric characters followed by an underscore (seurat/R/object.R)
  
  # Set defaults
  # Param k.filter: How many neighbors to use when filtering anchors
  if (is.null(k_filter)) k_filter = min(200, purrr::map_int(sc, ncol))
  # Param k.weight: Number of neighbors to consider when weighting anchors
  if (is.null(k_weight)) k_weight = min(100, purrr::map_int(sc, ncol))
  # Param k.anchor: How many neighbors to use when picking anchors
  if (is.null(k_anchor)) k_anchor = min(5, purrr::map_int(sc, ncol))
  
  # Param ndims: Number of dimensions cannot be larger than number of cells; note: Seurat gives an error that ndims must smaller than the number of cells
  ndims = min(c(ndims, purrr::map_int(sc, ncol)-1))
  
  # Param reference: split by comma and test if all entries are indices
  if (!is.null(reference)) {
    reference = trimws(unlist(strsplit(reference, ",")))
    if (any(is.na(suppressWarnings(as.numeric(reference))))) {
      # At least one entry is not an index, assume that the entries are names and match to get indices
      reference = match(reference, names(datasets))
      if (any(is.na(reference))) stop("RunIntegration: At least one reference dataset name could not be found in the datasets!")
      
    } else {
      # All entries are indices, just check that they are valid
      reference = as.numeric(reference)
      if (max(reference) > length(sc) | min(reference) < 1) stop("RunIntegration: Numeric index for the reference is either <1 or >length(samples)!")
    }
  }
  
  # Depending on assay RNA or SCT, integration is different
  if (assay == "RNA") {
    # Find RNA features that are consistently variable
    integrate_RNA_features = SelectIntegrationFeatures(object.list=sc,
                                                       nfeatures=3000,
                                                       verbose=verbose)
    
    # If reciprocal PCA is used for integration, scale integration features and use them to run PCA for each sample
    if (use_reciprocal_pca) {
      sc = purrr::map(sc, function(s) {
        s = Seurat::ScaleData(s,
                              assay="RNA",
                              features=integrate_RNA_features,
                              vars.to.regress=vars_to_regress,
                              verbose=verbose)
        s = Seurat::RunPCA(x,
                           assay="RNA",
                           npcs=min(c(50, ncol(s))),
                           features=integrate_RNA_features,
                           verbose=verbose)
        return(s)
      })
    }
    
    # Find integration anchors for assay RNA
    integrate_RNA_anchors = Seurat::FindIntegrationAnchors(object.list=sc, 
                                                           dims=1:ndims,
                                                           normalization.method="LogNormalize",
                                                           anchor.features=integrate_RNA_features, 
                                                           k.filter=k_filter,
                                                           k.anchor=k_anchor, 
                                                           reference=reference,
                                                           verbose=verbose,
                                                           reduction=ifelse(use_reciprocal_pca, "rpca", "cca"))
    
    # Integrate RNA data 
    sc = Seurat::IntegrateData(integrate_RNA_anchors, 
                               new.assay.name="RNAintegrated",
                               normalization.method="LogNormalize",
                               dims=1:ndims,
                               k.weight=k_weight,
                               verbose=verbose)
    
    rm(integrate_RNA_anchors, integrate_RNA_features)
  } else if (assay == "SCT") {
    # Find integration anchors for assay SCT
    integrate_SCT_features = Seurat::SelectIntegrationFeatures(object.list=sc, 
                                                       nfeatures=3000,
                                                       verbose=verbose)
    
    # If reciprocal PCA is used for integration, run SCTransform on the integration features and use them to run PCA for each sample
    # NOTE: This might be that we run SCTransform twice
    if (use_reciprocal_pca) {
      min_npcs = 
      sc = purrr::map(sc, function(s) {
        s = SCTransform(s, 
                    vars.to.regress=vars_to_regress, 
                    min_cells=min_cells, 
                    verbose=verbose, 
                    return.only.var.genes=FALSE,
                    method=ifelse(packages_installed("glmGamPoi"), "glmGamPoi", "poisson")) 
        s = Seurat::RunPCA(x,
                           assay="SCT",
                           npcs=min(c(50, ncol(s))),
                           features=integrate_SCT_features,
                           verbose=verbose)
        return(s)
      })
    }
    
    # Ensure that sctransform residuals have been computed for each feature and subset scale.data to only contain residuals for the integration features
    sc = Seurat::PrepSCTIntegration(object.list=sc, 
                                 anchor.features=integrate_SCT_features, 
                                 assay=rep("SCT",length(sc)),
                                 verbose=verbose)
    
    # Find integration anchors for assay SCT
    integrate_SCT_anchors = Seurat::FindIntegrationAnchors(object.list=sc,
                                                   dims=1:ndims, 
                                                   normalization.method="SCT", 
                                                   anchor.features=integrate_SCT_features, 
                                                   k.filter=k_filter,
                                                   k.anchor=k_anchor,
                                                   reference=reference,
                                                   verbose=verbose,
                                                   reduction=ifelse(use_reciprocal_pca, "rpca", "cca"))
    
    # Integrate SCT data
    sc = Seurat::IntegrateData(integrate_SCT_anchors, 
                                    new.assay.name="SCTintegrated",
                                    normalization.method="SCT",
                                    dims=1:ndims, 
                                    k.weight=k_weight,
                                    verbose=verbose)

    rm(integrate_SCT_features, integrate_SCT_anchors)
  }
  
  # Call garbage collector to free memory (hope it helps)
  gc(verbose=verbose)
  return(sc)
}

#' R implementation of the Seurat LogNormalize strategy but with additional options.
#' Normalisation: Counts for each cell are divided by the total counts for that cell and multiplied by a scale factor followed by natural-log transformation.
#' Note: It is likely slower than the C++ implementation of Seurat.
#' 
#' @param sc Seurat object.
#' @param assay The assay to normalize. Default is RNA.
#' @param slot The assay slot to use for normalization. Default is counts.
#' @param exclude_highly_expressed Whether to exclude highly expressed genes from size factor computation. Default is FALSE.
#' @param max_fraction Maximum fraction of counts for a gene not to be considered highly expressed. Default is 0.05.
#' @param exclude_genes A list of genes to exclude by default. Default is empty.
#' @param scale_factor The scale factor. Default is 10000.
#' @return A seurat object with normalized counts in the data slot of the assay.
LogNormalizeCustom = function(sc, assay="RNA", slot="counts", exclude_highly_expressed=FALSE, max_fraction=0.05, exclude_genes=NULL, scale_factor=10000) {
  # Get assay data
  counts_assay = Seurat::GetAssayData(sc, slot="counts", assay=assay)
  
  # For each cell: determine highly expressed genes and return the counts sum but excluding highly expressed genes and genes that should be excluded by default
  total_counts_per_cell_excl_high = apply(counts_assay, 2, function(cell_counts) {
    is_highly_expressed = exclude_highly_expressed & cell_counts/sum(cell_counts)>max_fraction
    is_excluded_by_default = names(cell_counts) %in% exclude_genes
    return(sum(cell_counts[!highly_expressed_genes & !is_excluded_by_default]))
  })
  
  # Size factors for scaling the counts
  size_factors = 1/(total_counts_per_cell_excl_high*scale_factor)
  
  # Multiply counts by size factors
  # Use R matrix multiplication for speed: matrix %*% diag(v)
  # To divide, multiply by the inverse
  # Finally log1p
  sc = Seurat::SetAssayData(sc, slot="data", 
                            assay=assay, 
                            new.data=log1p(Matrix::t(Matrix::t(counts_assay) * size_factors)))
  return(sc)
}
