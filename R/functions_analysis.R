#' Calculates the median of rows or columns of a sparse (dgCMatrix) or iterable (IterableMatrix) matrix. 
#'
#' @param matrix Sparse (dgCMatrix) or iterable (IterableMatrix) matrix
#' @param margin Margin for which to calculate median. Can be: 1 - rows, 2 - columns. Default is 1.
#' @param chunk_size Iterable matrices will be converted into sparse matrics. To avoid storing the entire matrix in memory, only process this number of rows/columns at once. Default is no chunks.
#' @return A named vector with medians.
CalculateMedians = function(matrix, margin=1, chunk_size=NULL){
  # Checks
  assertthat::assert_that(margin %in% c("1", "2"),
                          msg=FormatMessage("Margin can only be 1 - rows or 2 - columns."))
  
  # Define chunks
  if (margin == 1) {
    indices = 1:nrow(matrix)
  } else {
    indices = 1:ncol(matrix)
  }
  
  # Define chunks
  chunks = NULL
  if (!is.null(chunk_size)) {
    if (chunk_size < length(indices)) {
      chunks = split(indices, ceiling(seq_along(indices)/chunk_size))
    }
  }
  
  if (!is.null(chunks)) {
    # Calculate medians per chunk
    medians = purrr::map(chunks, function(c) {
      if (margin == 1) {
        mt = as(matrix[c, ], "dgCMatrix")
        mds = sparseMatrixStats::rowMedians(mt)
      } else {
        mt = as(matrix[, c], "dgCMatrix")
        mds = sparseMatrixStats::colMedians(mt)
      }
      return(mds)
    }) %>% purrr::flatten_dbl()
  } else {
    # Convert to sparse matrix
    if (!is(matrix, "dgCMatrix")) {
      matrix = as(matrix, "dgCMatrix")
    }
    
    # Calculate medians
    if (margin == 1) {
      medians = sparseMatrixStats::rowMedians(matrix)
    } else {
      medians = sparseMatrixStats::colMedians(matrix)
    }
  }
  
  return(medians)
}

#' Calculates the boxplot statistics for rows or columns of a sparse (dgCMatrix) or iterable (IterableMatrix) matrix. 
#' These are: min, q25, q50, q75, max, IQR (q75-q25), uppper whisker (q50 + 1.5 IQR), lower whisker (q50 - 1.5 IQR)
#'
#' @param matrix Sparse (dgCMatrix) or iterable (IterableMatrix) matrix
#' @param margin Margin for which to calculate median. Can be: 1 - rows, 2 - columns. Default is 1.
#' @param chunk_size Iterable matrices will be converted into sparse matrics. To avoid storing the entire matrix in memory, only process this number of rows/columns at once. Default is no chunks.
#' @return A named vector with medians.
CalculateBoxplotStats = function(matrix, margin=1, chunk_size=NULL){
  # Checks
  assertthat::assert_that(margin %in% c("1", "2"),
                          msg=FormatMessage("Margin can only be 1 - rows or 2 - columns."))
  
  # Define chunks
  if (margin == 1) {
    indices = 1:nrow(matrix)
  } else {
    indices = 1:ncol(matrix)
  }
  
  # Define chunks
  chunks = NULL
  if (!is.null(chunk_size)) {
    if (chunk_size < length(indices)) {
      chunks = split(indices, ceiling(seq_along(indices)/chunk_size))
    }
  }
  
  if (!is.null(chunks)) {
    # Calculate medians per chunk
    boxplot_stats = purrr::map_dfr(chunks, function(c) {
      if (margin == 1) {
        mt = as(matrix[c, ], "dgCMatrix")
        mds = sparseMatrixStats::rowQuantiles(mt)
      } else {
        mt = as(matrix[, c], "dgCMatrix")
        mds = sparseMatrixStats::colQuantiles(mt)
      }
      return(as.data.frame(mds))
    })
  } else {
    # Convert to sparse matrix
    if (!is(matrix, "dgCMatrix")) {
      matrix = as(matrix, "dgCMatrix")
    }
    
    # Calculate medians
    if (margin == 1) {
      boxplot_stats = sparseMatrixStats::rowQuantiles(matrix)
    } else {
      boxplot_stats = sparseMatrixStats::colQuantiles(matrix)
    }
    boxplot_stats = as.data.frame(boxplot_stats)
  }
  
  colnames(boxplot_stats) = c("min", "q25", "q50", "q75", "max")
  boxplot_stats$IQR = boxplot_stats$q75 - boxplot_stats$q25
  
  boxplot_stats$lower_whisker = purrr::pmap_dbl(boxplot_stats, function(min, q50, IQR, ...) {
    if (IQR>0) {
      lower_whisker = q50 - 1.5*IQR
      if (lower_whisker < min) lower_whisker = min
    } else {
      lower_whisker = min
    }
    return(lower_whisker)
  })
  
  boxplot_stats$upper_whisker = purrr::pmap_dbl(boxplot_stats, function(q50, max, IQR, ...) {
    if (IQR > 0){
      upper_whisker = q50 + 1.5*IQR
      if (upper_whisker > max) upper_whisker = max
    } else {
      upper_whisker = max
    }
    return(upper_whisker)
  })
  
  return(boxplot_stats)
}



#' Calculate cell cycle scores.
#' 
#' @param sc Seurat v5 object
#' @param genes_s Vector of gene names characteristic for S-phase
#' @param genes_g2m Vector of gene names characteristic for G2M-phase
#' @param assay Assay to use
#' @return Updated Seurat object with cell cycle phase as well as scores.
CCScoring = function(sc, genes_s, genes_g2m, assay=NULL){
  if (is.null(assay)) assay = Seurat::DefaultAssay(sc)
  
  # For each layer (dataset)
  # In this case, we need to split the Seurat object 
  # (since CellCycleScoring and AddModuleScore still cannot work with layers)
  sc_split = suppressMessages(Seurat::SplitObject(sc, split.by="orig.ident"))
    
  cell_cycle_scores = purrr::map_dfr(names(sc_split), function(n) {
    # Check that the genes exist
    genes_s_exists = genes_s %in% rownames(sc_split[[n]][[assay]])
    genes_g2m_exists = genes_g2m %in% rownames(sc_split[[n]][[assay]])
    
    
    if (sum(genes_s_exists) >= 20 & sum(genes_g2m_exists) >= 20) {
      sc_split[[n]] = Seurat::CellCycleScoring(sc_split[[n]],
                                               s.features=genes_s[genes_s_exists],
                                               g2m.features=genes_g2m[genes_g2m_exists],
                                               assay=assay,
                                               verbose=FALSE)
      cc_scores = sc_split[[n]][[c("Phase", "S.Score", "G2M.Score")]]
      cc_scores[["CC.Difference"]] = cc_scores[["S.Score"]] - cc_scores[["G2M.Score"]]
    } else {
      cc_scores = data.frame(Phase=character(), S.Score=numeric(), G2M.Score=numeric(), CC.Difference=numeric())
    }
    cc_scores[["Phase"]] = factor(cc_scores[["Phase"]], levels=c("G1", "G2M", "S"))
    return(cc_scores)
  })
    
  # Add to barcode metadata
  sc = Seurat::AddMetaData(sc, cell_cycle_scores)
  
  # Add to 'gene_lists' slot in the misc slot of the Seurat object
  sc = ScAddLists(sc, lists=list(CC_S_phase=genes_s, CC_G2M_phase=genes_g2m),
                  lists_slot="gene_lists")
  
  # Log command
  sc = Seurat::LogSeuratCommand(sc)
  
  return(sc)
}

#' Apply scran normalization (using pooled size factors) to counts data.
#' 
#' @param sc Seurat v5 object.
#' @param assay Assay to normalize. If NULL, will be default assay of the Seurat object.
#' @param layer Layer to normalize. Default is "counts" meaning all raw counts layers.
#' @param save Name of the normalized layers. Default is "data" meaning new layers will be named data.X, data.Y, ...
#' @param chunk_size Maximum number of barcodes for which to compute size factors at once. Large counts matrices will be split into chunks to save memory.
#' @return Seurat v5 object with new layers data.X, data.Y, ...
NormalizeDataScran = function(sc, assay=NULL, layer="counts", save="data", chunk_size=50000) {
  if (is.null(assay)) assay = Seurat::DefaultAssay(sc)
  
  #suppressWarnings({Seurat::Misc(sc[[assay]], slot="scran_size_factors") = c()})
  
  # Iterate over layers (datasets) and get size factors
  olayers = layers = unique(layer)
  layers = SeuratObject::Layers(sc[[assay]], layer)
  if (length(save) != length(layers)) {
    save = make.unique(names=gsub(pattern=olayers, replacement=save, x=layers))
  }
  
  for (i in seq_along(layers)) {
    l = layers[i]
    
    barcodes = SeuratObject::Cells(sc[[assay]], layer=l)
    chunks = split(barcodes, ceiling(seq_along(barcodes)/chunk_size))
    
    # Calculate size factors per chunk
    size_factors = purrr::map(chunks, function(c) {
      counts = SeuratObject::LayerData(sc[[assay]], layer=l, fast=NA, cells=c)
      
      # Convert to dgCMatrix and create
      if (!is(counts, "dgCMatrix")) {
        counts = as(counts, "dgCMatrix")
      }
      
      # Convert to SingleCellExperiment, cluster and compute size factors per cluster
      # Note: These are already centered
      sce = SingleCellExperiment::SingleCellExperiment(list(counts=counts))
      clusters = scran::quickCluster(sce, min.size=100)
      sce = scuttle::computePooledFactors(sce, clusters=clusters)
      sf = scuttle::pooledSizeFactors(sce, clusters=clusters)
      sf = setNames(sf, c)
      return(sf)
    })
    size_factors = purrr::flatten(size_factors) %>% unlist()
    
    # Store size factors in misc slot of assay
    #misc_lst = Seurat::Misc(sc[[assay]], slot="scran_size_factors")
    #misc_lst = c(misc_lst, size_factors)
    #suppressWarnings({Seurat::Misc(sc[[assay]], slot="scran_size_factors") = misc_lst})
    
    # Now apply size factors and normalise
    counts = SeuratObject::LayerData(sc[[assay]], layer=l, fast=NA)
    counts = log1p(Matrix::t(Matrix::t(counts) / size_factors))
    
    LayerData(sc[[assay]], 
              layer=save[i], 
              features=SeuratObject::Features(sc[[assay]], layer=l),
              cells=SeuratObject::Cells(sc[[assay]], layer=l)) = counts
    
  }
  
  # Log command
  sc = Seurat::LogSeuratCommand(sc)
  
  return(sc)
}

#' Identify highly variable features with the scran method (mean - var analysis).
#' 
#' @param sc Seurat v5 object.
#' @param assay Assay to analyze. If NULL, will be default assay of the Seurat object.
#' @param nfeatures Number of features to identify. Default is 2000.
#' @param combined If TRUE, analyze all data together. Recommended but requires more memory. If FALSE, features will be identified by dataset and sets will then be combined.
#' @return Seurat v5 object with highly variable features for assay.
FindVariableFeaturesScran = function(sc, assay=NULL, nfeatures=2000, combined=TRUE) {
  if (is.null(assay)) assay = Seurat::DefaultAssay(sc)
  
  # Checks
  layers = SeuratObject::Layers(sc[[assay]], "data")
  assertthat::assert_that(length(layers) > 0,
                          msg=FormatMessage("Could not find normalized data for assay {assay}."))
  
  # Get normalized data
  if (combined) {
    # Find features using all layers (dataset) together
    data = purrr::map(layers, function(l) {
      dt = SeuratObject::LayerData(sc, assay=assay, layer=l)
      
      # Convert to dgCMatrix matrix
      if (!is(dt, "dgCMatrix")) {
        dt = as(dt, "dgCMatrix")
      }
      return(dt)
    })
    
    # Record to which dataset a cell belongs
    block = purrr::map(seq_along(layers), function(i) {
      b = rep(layers[i], ncol(data[[i]]))
      return(b)
    }) %>% purrr::flatten_chr()
    block = factor(block, levels=unique(block))
    
    # Find variable features using the entire dataset
    data = purrr::reduce(data, cbind)
    hvf_info = scran::modelGeneVar(data, block=block)
    hvf = scran::getTopHVGs(hvf_info, n=nfeatures)
    
    # Combine hvf_info tables
    hvf_info = purrr::map(layers, function(l) {
      hvfi = hvf_info[["per.block"]][[l]]
      h = scran::getTopHVGs(hvfi, n=nfeatures)
      
      hvfi = as.data.frame(hvfi)
      hvfi$variable = rownames(hvfi) %in% h
      hvfi$rank = match(rownames(hvfi), h)
      colnames(hvfi) = paste0("vf_scran_", l, "_", colnames(hvfi))
      return(hvfi)
    })
    hvf_info = dplyr::bind_cols(hvf_info_lst)
  } else {
    # Find for each layer (dataset) separately
    
    # Find highly variable features per layer (dataset)
    hvf_info_lst = purrr::map(layers, function(l) {
      data = SeuratObject::LayerData(sc, assay=assay, layer=l)
      
      # Convert to dgCMatrix matrix
      if (!is(data, "dgCMatrix")) {
        data = as(data, "dgCMatrix")
      }
      
      # Find highly variable genes with scran
      hvf_info = scran::modelGeneVar(data)
      hvf = scran::getTopHVGs(hvf_info, n=nfeatures)
      
      hvf_info = as.data.frame(hvf_info)
      hvf_info$variable = rownames(hvf_info) %in% hvf
      hvf_info$rank = match(rownames(hvf_info), hvf)
      
      colnames(hvf_info) = paste0("vf_scran_", l, "_", colnames(hvf_info))
      return(hvf_info)
    })
    names(hvf_info_lst) = layers
    
    # Get variable features per layer (dataset) and rank
    hvf_lst =purrr::map(layers, function(l) {
      var_col = paste0("vf_scran_", l, "_variable")
      rank_col = paste0("vf_scran_", l, "_rank")
      hvf = hvf_info_lst[[l]] %>% 
        dplyr::filter(!!sym(var_col)) %>%
        dplyr::arrange(!!sym(rank_col)) %>%
        rownames()
      return(hvf)
    })
    
    # Find a good overall set (from SeuratObject:::VariableFeatures.StdAssay)
    hvf = SeuratObject:::.SelectFeatures(hvf_lst, 
                                         all.features=intersect(x = slot(sc[[assay]], name = "features")[, layers]), 
                                         nfeatures = nfeatures)
    
    # Combine hvf_info tables
    hvf_info = dplyr::bind_cols(hvf_info_lst)
  }
  
  # Add to feature metadata
  hvf_info$var.features = rownames(hvf_info) %in% hvf
  hvf_info$var.features.rank = match(rownames(hvf_info), hvf)
  sc[[assay]] = SeuratObject::AddMetaData(sc[[assay]], metadata=hvf_info)
  
  # Set variable features
  SeuratObject::VariableFeatures(sc[[assay]]) = hvf
  
  # Log command
  sc = Seurat::LogSeuratCommand(sc)
  
  return(sc)
}







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
#' @param k_score How many neighbors for calculating scores. If NULL, automatically set to min(30, minimum number of cells in a sample)).
#' @param vars_to_regress For reciprocal PCA: when doing the scaling, which variables should be regressed out when doing the scaling
#' @param min_cells For reciprocal PCA: when doing the scaling for SCT, the minimum number of cells a gene should be expressed
#' @return A Seurat object with an integrated assay (RNAintegrated or SCTintegrated) and a merged assay (RNA or SCT)
RunIntegration = function(sc, ndims=30, reference=NULL, use_reciprocal_pca=FALSE, verbose=FALSE, assay="RNA", k_filter=NULL, k_weight=NULL, k_anchor=NULL, k_score=NULL, vars_to_regress=NULL, min_cells=1) {
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
  # Param k_score: How many neighbors to use fo calculating scores
  if (is.null(k_score)) k_score = min(30, purrr::map_int(sc, ncol))
  
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
                                                           k.score=k_score,
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
                                                   k.score=k_score,
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
LogNormalizeCustom = function(sc, assay="RNA", exclude_highly_expressed=FALSE, max_fraction=0.05, exclude_genes=NULL, scale_factor=10000) {
  # Get assay data
  counts = Seurat::GetAssayData(sc, slot="counts", assay=assay)
  
  # Determine highly expressed genes
  if (exclude_highly_expressed) { 
    
    # Calculate fractions for genes per cell, as a TRUE/FALSE table
    is_highly_expressed = Matrix::t(Matrix::t(counts) * (1/Matrix::colSums(counts))) > max_fraction
    
    # Which genes are highly expressed in at least one cell, as a TRUE/FALSE vector
    is_highly_expressed = Matrix::rowSums(is_highly_expressed) > 0
    
  } else { 
    is_highly_expressed = FALSE
  }
  
  # Convert excluded genes into TRUE/FALSE vector
  is_excluded_by_default = rownames(counts) %in% exclude_genes
  
  # For each cell return counts sum excluding genes that are highly expressed or should be excluded anyhow
  total_counts_per_cell = Matrix::colSums(counts[!is_highly_expressed & !is_excluded_by_default, ])
  
  # Size factors for scaling the counts
  size_factors = 1 / total_counts_per_cell * scale_factor
  
  # Multiply counts by size factors
  # Use R matrix multiplication for speed: matrix %*% diag(v)
  # To divide, multiply by the inverse
  # Finally log1p
  sc = Seurat::SetAssayData(sc, 
                            slot="data", 
                            assay=assay, 
                            new.data=log1p(Matrix::t(Matrix::t(counts) * size_factors)))
  return(sc)
}