#' Our plotting style.
#'
#' @param title The plot title.
#' @param col A vector of colours to use.
#' @param fill A vector of fill colours to use.
#' @param legend_title The legend title.
#' @param legend_position The legend position.
#' @param xlab The title of the x-axis.
#' @param ylab The title of the y-axis.
#' @return None, add as theme.
AddPlotStyle = function(title=NULL, col=NULL, fill=NULL, legend_title=NULL, legend_position=NULL, xlab=NULL, ylab=NULL) {
  style = list(
    # Basic theme
    ggplot2::theme_light(11),
    ggplot2::theme(panel.border = element_blank())
  )
  
  # Title
  if (!is.null(title)) {
    if (nchar(title) > 0) {
      style = c(style, list(ggplot2::ggtitle(title)))
    } else {
      style = c(style, list(ggplot2::theme(title=element_blank())))
    }
  }
  
  # Colour
  if (!is.null(col)) style = c(style, list(ggplot2::scale_colour_manual(values=col)))
  
  # Fill
  if (!is.null(fill)) style = c(style, list(ggplot2::scale_fill_manual(values=fill)))
    
  # Legend title
  if (!is.null(legend_title)) {
    if (nchar(legend_title) > 0) {
      style = c(style, list(labs(color=legend_title, fill=legend_title)))
    } else {
      style = c(style, list(ggplot2::theme(legend.title=element_blank())))
    }
  }
  
  # Legend position
  if (!is.null(legend_position)) style = c(style, list(ggplot2::theme(legend.position=legend_position)))
  
  # Axis labels
  if (!is.null(xlab)) {
    if (nchar(xlab) > 0) {
      style = c(style, list(ggplot2::xlab(xlab)))
    } else {
      style = c(style, list(ggplot2::theme(axis.title.x=element_blank())))
    }
  }
  
  if (!is.null(ylab)) {
    if (nchar(ylab) > 0) {
      style = c(style, list(ggplot2::ylab(ylab)))
    } else {
      style = c(style, list(ggplot2::theme(axis.title.y=element_blank())))
    }
  }
  
  return(style)
}

#' Helper function to generate plot captions.
#'
#' @param plot_names A list of plot names.
#' @param remove Part to remove before matching. Can be a string or regex.
#' @param split If not NULL, split plot names to generate a X vs Y caption.
#' @param capitalize If TRUE, capitalize the first letter of the caption.
#' @return A list of captions.
GeneratePlotCaptions = function(plot_names, remove=NULL, split=NULL, capitalize=TRUE) {
  # Split names into two parts if requested
  if (!is.null(split)) {
    plot_names = strsplit(plot_names, split)
  } else {
    plot_names = as.list(plot_names)
  }
  
  # Remove parts
  if (!is.null(remove)) {
    plot_names = purrr::map(plot_names, function(x) {
      return(gsub(pattern=remove, replacement="", x))
    })
  }
  
  # Match and generate captions
  captions = purrr::map(plot_names, function(x) {
    # Add captions here
    capts = dplyr::case_match(x,
                                 "nCount" ~ "number of counts",
                                 "nFeature" ~ "number of features",
                                 "pCountsTop50" ~ "percent counts in top50 features",
                                 "pMito" ~ "percent counts in mitochondrial genes",
                                 "pRibosomal" ~ "percent counts in ribosomal genes",
                                 "pGlobin" ~ "percent counts in globin genes",
                                 "pERCC" ~ "percent counts in ERCC controls",
                                 "pXIST" ~ "percent counts in XIST gene",
                                 "pChrY" ~ "percent counts in chrY genes",
                                 "S.Score" ~ "S phase score",
                                 "G2M.Score" ~ "G2M phase score",
                                 "Phase" ~ "cell cycle phase",
                                 "varFeatures" ~ "variable features",
                                 "Expression" ~ "Expression of",
                                 .default = NULL
    )
    idx = which(is.na(capts))
    capts[idx] = paste0('Metric "', x[idx], '"')
    return(capts)
  })
  
  # If there are two parts, collapse with " Vs "
  captions = purrr::map(captions, function(x) {
    for(i in seq_along(x)) {
      # Skip first one
      if (i==1) next
      
      # Make first character lower-case unless the first two characters are upper-case (meaning an abbreviation)
      # if (!grepl("^[[:upper:]]{2,}", x[i])) {
      #   s = unlist(strsplit(x[i], ""))
      #   s[1] = tolower(s[1])
      #   x[i] = paste(s, collapse="")
      # }
    }
    
    # Collapse
    x = paste(x, collapse=" Vs ")
    
    # If requested, make first letter upper-case
    if (capitalize) {
      s = unlist(strsplit(x, ""))
      s[1] = toupper(s[1])
      x = paste(s, collapse="")
    }
    
    return(x)
  }) %>% unlist()
  
  return(captions)
}

#' Plots barcode metadata for QC. Numeric columns will plotted as violin plots and non-numeric
#' columns will be plotted as bar plots.
#'
#' @param sc Seurat v5 object.
#' @param qc Barcode metadata columns to plot.
#' @param filter A nested list where the first level is the barcode metadata column and the second levels 
#' contains filters per dataset. Filters for numeric columns must numeric vectors with min and max. Filter
#' for character/factor columns must be character vectors with the values that should be kept. 
#' level contains the filter values 
#' @return A list of ggplot2 objects.
PlotBarcodeQC = function(sc, qc, filter=NULL) {
  barcode_metadata = sc[[]]
  
  # Determine QC type (numeric or non-numeric)
  is_numeric = purrr::map_lgl(qc, function(q) return(is.numeric(barcode_metadata[, q, drop=TRUE])))
  
  # Get filter thresholds per QC metrics (numeric)
  qc_thresholds_numeric = purrr::map(qc[is_numeric], function(f) {
    tresh = purrr::map_dfr(names(filter), function(n) {
      tr = data.frame(qc_feature=character(), ident=character(), 
                      threshold=character(), value=numeric())
      
      if (f %in% names(filter[[n]])) {
        tr = data.frame(qc_feature=f, 
                        ident=n, 
                        min=filter[[n]][[f]][1], 
                        max=filter[[n]][[f]][2])  %>% 
          tidyr::pivot_longer(c(min, max), names_to="threshold", values_to="value")
      }
      
      tr$ident = factor(tr$ident, levels=orig_idents)
      return(tr)
    })
    return(tresh)
  })
  names(qc_thresholds_numeric) = qc[is_numeric]
  
  qc_thresholds_other = purrr::map(qc[!is_numeric], function(f) {
    tresh = purrr::map(names(filter), function(n) {
      return(filter[[n]][[f]])
    })
    names(tresh) = names(filter)
    return(tresh)
  })
  names(qc_thresholds_other) = qc[!is_numeric]
  
  # Plot numeric QC
  plist_numeric = Seurat::VlnPlot(sc, features=qc[is_numeric], combine=FALSE, pt.size=0, raster=FALSE, layer="counts")
  names(plist_numeric) = qc[is_numeric]
  
  plist_numeric = purrr::map(names(plist_numeric), function(n) {
    # Add style
    p = plist_numeric[[n]] + 
      AddPlotStyle(legend_position="none", xlab="", fill=ScColours(sc, "orig.ident")) +
      theme(axis.text.x=element_text(angle=45, hjust=1))
    
    # Add filter thresholds
    qc_threshold_segments = purrr::pmap(qc_thresholds_numeric[[n]], function(qc_feature, ident, threshold, value) {
      return(annotate(geom="segment", x=as.integer(ident)-0.5, xend=as.integer(ident)+0.5, y=value, yend=value, linetype=ifelse(threshold=="max", 6, 2), colour="firebrick"))
    })
    p = p + qc_threshold_segments
    return(p)
  })
  names(plist_numeric) = qc[is_numeric]
  
  # Plot non-numeric QC
  plist_other = purrr::map(qc[!is_numeric], function(n) {
    # Collect plot_data
    plot_data = barcode_metadata %>% 
      dplyr::select(x=orig.ident, y=!!sym(n)) %>%
      dplyr::count(x, y) %>%
      dplyr::group_by(x) %>%
      dplyr::mutate(perc=n/sum(n)*100)
    plot_data$status = "filter"
    for(i in 1:nrow(plot_data)){
      allowed_values = qc_thresholds_other[[n]]
      allowed_values = allowed_values[[plot_data$x[i]]]
      plot_data$status[i] = ifelse(plot_data$y[i] %in% allowed_values, "keep", "filter")
    }
    plot_data$status = factor(plot_data$status, levels=c("keep", "filter"))
    
    # Make plot
    p = ggplot(plot_data, aes(x=x, y=perc, fill=y, alpha=status, col=status)) +
      geom_bar(stat="identity", position="stack") +
      scale_x_discrete("Identity") +
      scale_color_manual(values=c("keep"="black", "filter"="grey")) +
      scale_alpha_manual(values=c("keep"=1, "filter"=0.2)) +
      AddPlotStyle(fill=ScColours(sc, "orig.ident")) +
      theme(axis.title.y=element_blank())
    
    return(p)
  })
  names(plist_other) = qc[!is_numeric]
  
  # Adjust order
  plist = c(plist_numeric, plist_other)
  plist = plist[qc]
  
  return(plist)
}

#' Plots two barcode metadata columns for QC. Supports only numeric columns.
#'
#' @param sc Seurat v5 object.
#' @param qc Pairs of barcode metadata columns to plot.
#' @param filter A nested list where the first level is the barcode metadata column and the second levels 
#' contains filters per dataset. Filters for numeric columns must numeric vectors with min and max. Filter
#' for character/factor columns must be character vectors with the values that should be kept.
#' @param raster Whether to raster the points.
#' @return A list of ggplot2 objects.
PlotBarcodeQCCor = function(sc, qc, filter=NULL, raster=FALSE) {
  barcode_metadata = sc[[]]
  
  # Check QC type (only numeric allowed)
  qc_cols = purrr::flatten(qc) %>%
    unlist() %>%
    unique()
  f = purrr::map_lgl(qc_cols, function(c) return(is.numeric(barcode_metadata[, c, drop=TRUE])))
  assertthat::assert_that(all(f),
                          msg=FormatMessage("Barcode metadata columns {qc_cols[!f]*} are not numeric! Function PlotBarcodeQCCor can only plot numeric data."))
  
  # Get filter thresholds per QC metrics (numeric)
  qc_thresholds = purrr::map(qc_cols, function(f) {
    tresh = purrr::map_dfr(names(filter), function(n) {
      tr = data.frame(qc_feature=character(), ident=character(), 
                      threshold=character(), value=numeric())
      
      if (f %in% names(filter[[n]])) {
        tr = data.frame(qc_feature=f, 
                        ident=n, 
                        min=filter[[n]][[f]][1], 
                        max=filter[[n]][[f]][2])  %>% 
          tidyr::pivot_longer(c(min, max), names_to="threshold", values_to="value")
      }
      
      tr$ident = factor(tr$ident, levels=orig_idents)
      return(tr)
    })
    return(tresh)
  })
  names(qc_thresholds) = qc_cols
  
  plist = purrr::map(qc, function(c) {
    f1 = c[1]
    f2 = c[2]
    
    # Plot QC feature f1 vs f2
    p = Seurat::FeatureScatter(sc, feature1=f1, feature2=f2, shuffle=TRUE, seed=getOption("random_seed"), raster=raster)
    p = p + AddPlotStyle(col=ScColours(sc, "orig.ident"))
    
    # Add filter thresholds for f1
    qc_threshold_segments = purrr::pmap(qc_thresholds[[f1]], function(qc_feature, ident, threshold, value) {
      return(geom_vline(xintercept=value, linetype=ifelse(threshold=="max", 6, 2), colour="firebrick"))
    })
    p = p + qc_threshold_segments
    
    # Add filter thresholds for f2
    qc_threshold_segments = purrr::pmap(qc_thresholds[[f2]], function(qc_feature, ident, threshold, value) {
      return(geom_hline(yintercept=value, linetype=ifelse(threshold=="max", 6, 2), colour="firebrick"))
    })
    p = p + qc_threshold_segments
    
    return(p)
  })
  names(plist) = purrr::map(qc, paste, collapse="__") %>% unlist()
  
  return(plist)
}

#' Plots the variable features for each layer (dataset).
#'
#' @param sc Seurat v5 object.
#' @param method Method used to find variable features. Can be: 'vst' (Seurat standard), 'sct' (SCTransform) or 'scran' (Scran).
#' @param assay Assay. If NULL, defaults to default assay of the Seurat object.
#' @param top The top genes that should be labeled.
#' @return A list of ggplot2 objects.
PlotVariableFeatures = function(sc, method, assay=NULL, top=10) {
  if (is.null(assay)) assay = Seurat::DefaultAssay(sc)
  
  # Checks
  valid_methods = c("vst", "sct", "scran")
  assertthat::assert_that(method %in% valid_methods,
                          msg=FormatMessage("Method is {method} but must be one of: {valid_methods*}."))
  
  layers = SeuratObject::Layers(sc[[assay]], "^data")
  assertthat::assert_that(length(layers) > 0,
                          msg=FormatMessage("Could not find normalized data for assay {assay}."))
  
  # Make plots per layer (dataset)
  orig_idents =  levels(sc$orig.ident)
  
  plist = purrr::map(orig_idents, function(n) {
    # Collect information about highly variable genes
    if (method == "sct") {
      assertthat::assert_that(.hasSlot(sc[[assay]], "SCTModel.list"),
                              msg=FormatMessage("No variable feature information in slot SCTModel.list available for {assay}."))
      hvf_info = Seurat::SCTResults(sc[[assay]], slot="feature.attributes", model=n)
      hvf_info = hvf_info[, c("gmean", "variance", "residual_variance")]
      hvf_info$variable = rownames(hvf_info) %in% Seurat::VariableFeatures(sc[[assay]])
      hvf_info$rank = match(rownames(hvf_info), Seurat::VariableFeatures(sc[[assay]]))
    } else if (method == "vst") {
      hvf_info = SeuratObject::HVFInfo(sc[[assay]],
                                       method="vst",
                                       layer=paste("data", n, sep="."),
                                       status=TRUE)
      assertthat::assert_that(!is.null(hvf_info),
                              msg=FormatMessage("No variable feature information available for {assay}."))
    } else if (method == "scran") {
      hvf_info = SeuratObject::HVFInfo(sc[[assay]],
                                       method="scran",
                                       layer=paste("data", n, sep="."),
                                       status=TRUE)
      assertthat::assert_that(!is.null(hvf_info),
                              msg=FormatMessage("No variable feature information available for {assay}."))
      
    }
    
    # Get top genes
    top_genes = which(hvf_info$rank<=top)
    top_genes = rownames(hvf_info)[top_genes]
    
    # Define columns to plot
    if (method == "scran") {
      hvf_info = hvf_info[, c(1, 2, ncol(hvf_info)-1, ncol(hvf_info))]
      xlab = "Average Expression"
      ylab = "Dispersion"
    } else if (method == "vst") {
      hvf_info = hvf_info[, c(1, 4, ncol(hvf_info)-1, ncol(hvf_info))]
      xlab = "Average Expression"
      ylab = "Standardized Variance"
    } else if (method == "sct") {
      hvf_info = hvf_info[, c(1, 3, ncol(hvf_info)-1, ncol(hvf_info))]
      xlab = "Geometric Mean of Expression"
      ylab = "Residual Variance"
    }
    
    xvar = colnames(hvf_info)[1]
    yvar = colnames(hvf_info)[2]
    
    # Make plot
    p = ggplot(hvf_info, aes(x=!!sym(xvar), y=!!sym(yvar), col=variable)) +
      geom_point() +
      scale_x_log10(xlab) +
      scale_y_continuous(ylab) +
      scale_color_manual("Variable", values=c("black", "red")) +
      AddPlotStyle() + 
      theme(legend.position="none", legend.background=element_rect(fill=alpha("white", 0.0)))
    p = Seurat::LabelPoints(plot=p, points=top_genes, repel=TRUE, xnudge=0, ynudge=0)
    return(p)
  })
  names(plist) = paste("varFeatures", orig_idents, sep="_")
  
  return(plist)
}

#' Plots the relative log expression features for each layer (dataset).
#'
#' @param sc Seurat v5 object.
#' @param assay Assay. If NULL, defaults to default assay of the Seurat object.
#' @param layer Type of layer. Can be counts or data but also a specific layer.
#' @param nbarcodes Number of barcodes to plot.
#' @param is_log Is the data already in log? If not, will be logged.
#' @return A ggplot2 object.
PlotRLE = function(sc, assay=NULL, layer="counts", nbarcodes=500, is_log=FALSE) { 
  if (is.null(assay)) assay = Seurat::DefaultAssay(sc)
  
  # Checks
  layers = SeuratObject::Layers(sc[[assay]], layer)
  assertthat::assert_that(length(layers) > 0,
                          msg=FormatMessage("Could not find data for layer {{layer}} of assay {assay}."))
  
  if (!is(sc[[assay]], "SCTAssay")) {
    # Standard assays
    
    # Sample
    data = purrr::map(layers, function(l) {
      dt = SeuratObject::LayerData(sc, assay=assay, layer=l)
      set.seed(getOption("random_seed"))
      idx = 1:ncol(dt)
      smp = sample(idx, min(length(idx), nbarcodes))
      return(dt[, smp])
    })
    
    # Get identity
    orig_idents = purrr::map(seq_along(layers), function(i) {
      return(rep(layers[i], ncol(data[[i]])))
    }) %>% unlist()
    orig_idents = gsub(pattern=paste0(layer,"."), replacement="", x=orig_idents)
    orig_idents = factor(orig_idents, levels=levels(sc$orig.ident))
    
    # Merge
    data = purrr::reduce(data, cbind)
  } else {
    # SCT assay
    
    # Get data
    data = SeuratObject::LayerData(sc, assay=assay, layer=layer)
    
    # Sample
    barcode_metadata = sc[[]]
    barcode_metadata = barcode_metadata[, ]
    orig_idents = split(colnames(data), barcode_metadata$orig.ident)
    orig_idents = purrr::map(orig_idents, function(bcs) {
      set.seed(getOption("random_seed"))
      bcs = sample(bcs, min(length(bcs), nbarcodes))
      return(bcs)
    })
    
    # Subset data
    data = data[, orig_idents %>% purrr::flatten_chr()]
    data = as(data, "dgCMatrix")
    
    # Get identity
    orig_idents = purrr::map(names(orig_idents), function(n) {
      return(rep(n, length(orig_idents[[n]])))
    }) %>% purrr::flatten_chr()
    orig_idents = factor(orig_idents, levels=levels(sc$orig.ident))
  }
    
  if (!is_log) {
    data = log1p(data + 1)
  }
  
  # Calculate feature medians across all barcodes
  row_medians = CalculateMedians(matrix=data, margin=1, chunk_size=NULL)
  
  # Subtract gene median from gene count
  data = data - row_medians
  
  # Calculate RLE stats per cell
  rle_stats = CalculateBoxplotStats(data, margin=2, chunk_size=NULL)
  rle_stats$orig.ident = orig_idents
  rle_stats = rle_stats %>% 
    dplyr::arrange(orig.ident)
  rle_stats$x = 1:nrow(rle_stats)
  
  
  # Actual plotting
  p = ggplot(rle_stats) +
    geom_segment(aes(x=x, xend=x, y=lower_whisker , yend=upper_whisker, col=orig.ident)) +
    geom_segment(aes(x=x, xend=x, y=q25-0.01 , yend=q75+0.01), colour="grey20") +
    geom_point(aes(x=x, y=q50), shape=1) +
    AddPlotStyle(xlab="Cells", ylab="Relative log expression", col=ScColours(sc, "orig.ident")) + 
    theme(axis.text.x=element_blank(),
          axis.ticks.x=element_blank())
  
  return(p)
}

#' Transform a matrix cells (rows) x htos (cols) into a format that can be understood by feature_grid: cell class, name hto1, value hto1, name hto2, value hto2
#' 
#' @param x: A matrix cells (rows) x htos (cols).
#' @param cell_classification A vector of cell classifications.
DfAllColumnCombinations = function(x, cell_classification) {
  out = combn(x, 2, simplify=FALSE)
  out = lapply(out, function(o) {
    return(data.frame(cell_classification=unname(cell_classification[rownames(o)]), name1=colnames(o)[1], value1=o[, 1], name2=colnames(o)[2], value2=o[, 2]))
  }) %>% dplyr::bind_rows()
  
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

# Re-define Dotplot for non-scaled dotplots
# Keep this code in here until Seurat fixes this
# Issue: https://github.com/satijalab/seurat/issues/4298
PercentAbove <- function(x, threshold) {
  return(length(x = x[x > threshold]) / length(x = x))
}

DotPlotUpdated = function (object, assay = NULL, features, cols = c("lightgrey", 
                                                                    "blue"), col.min = -2.5, col.max = 2.5, dot.min = 0, dot.scale = 6, 
                           idents = NULL, group.by = NULL, split.by = NULL, cluster.idents = FALSE, 
                           scale = TRUE, scale.by = "radius", scale.min = NA, scale.max = NA) 
{
  # assay <- assay %||% DefaultAssay(object = object)
  assay <- DefaultAssay(object = object)
  DefaultAssay(object = object) <- assay
  split.colors <- !is.null(x = split.by) && !any(cols %in% 
                                                   rownames(x = brewer.pal.info))
  scale.func <- switch(EXPR = scale.by, size = scale_size, 
                       radius = scale_radius, stop("DotPlotUpdated: 'scale.by' must be either 'size' or 'radius'"))
  feature.groups <- NULL
  if (is.list(features) | any(!is.na(names(features)))) {
    feature.groups <- unlist(x = sapply(X = 1:length(features), 
                                        FUN = function(x) {
                                          return(rep(x = names(x = features)[x], each = length(features[[x]])))
                                        }))
    if (any(is.na(x = feature.groups))) {
      warning("Some feature groups are unnamed.", call. = FALSE, 
              immediate. = TRUE)
    }
    features <- unlist(x = features)
    names(x = feature.groups) <- features
  }
  cells <- unlist(x = CellsByIdentities(object = object, idents = idents))
  data.features <- FetchData(object = object, vars = features, 
                             cells = cells)
  data.features$id <- if (is.null(x = group.by)) {
    Idents(object = object)[cells, drop = TRUE]
  }
  else {
    object[[group.by, drop = TRUE]][cells, drop = TRUE]
  }
  if (!is.factor(x = data.features$id)) {
    data.features$id <- factor(x = data.features$id)
  }
  id.levels <- levels(x = data.features$id)
  data.features$id <- as.vector(x = data.features$id)
  if (!is.null(x = split.by)) {
    splits <- object[[split.by, drop = TRUE]][cells, drop = TRUE]
    if (split.colors) {
      if (length(x = unique(x = splits)) > length(x = cols)) {
        stop("DotPlotUpdated: Not enough colors for the number of groups")
      }
      cols <- cols[1:length(x = unique(x = splits))]
      names(x = cols) <- unique(x = splits)
    }
    data.features$id <- paste(data.features$id, splits, sep = "_")
    unique.splits <- unique(x = splits)
    id.levels <- paste0(rep(x = id.levels, each = length(x = unique.splits)), 
                        "_", rep(x = unique(x = splits), times = length(x = id.levels)))
  }
  data.plot <- lapply(X = unique(x = data.features$id), FUN = function(ident) {
    data.use <- data.features[data.features$id == ident, 
                              1:(ncol(x = data.features) - 1), drop = FALSE]
    avg.exp <- apply(X = data.use, MARGIN = 2, FUN = function(x) {
      return(mean(x = expm1(x = x)))
    })
    pct.exp <- apply(X = data.use, MARGIN = 2, FUN = PercentAbove, 
                     threshold = 0)
    return(list(avg.exp = avg.exp, pct.exp = pct.exp))
  })
  names(x = data.plot) <- unique(x = data.features$id)
  if (cluster.idents) {
    mat <- do.call(what = rbind, args = lapply(X = data.plot, 
                                               FUN = unlist))
    mat <- scale(x = mat)
    id.levels <- id.levels[hclust(d = dist(x = mat))$order]
  }
  data.plot <- lapply(X = names(x = data.plot), FUN = function(x) {
    data.use <- as.data.frame(x = data.plot[[x]])
    data.use$features.plot <- rownames(x = data.use)
    data.use$id <- x
    return(data.use)
  })
  data.plot <- do.call(what = "rbind", args = data.plot)
  if (!is.null(x = id.levels)) {
    data.plot$id <- factor(x = data.plot$id, levels = id.levels)
  }
  if (length(x = levels(x = data.plot$id)) == 1) {
    scale <- FALSE
    warning("Only one identity present, the expression values will be not scaled", 
            call. = FALSE, immediate. = TRUE)
  }
  avg.exp.scaled <- sapply(X = unique(x = data.plot$features.plot), 
                           FUN = function(x) {
                             data.use <- data.plot[data.plot$features.plot == 
                                                     x, "avg.exp"]
                             if (scale) {
                               data.use <- scale(x = data.use)
                               data.use <- MinMax(data = data.use, min = col.min, 
                                                  max = col.max)
                             }
                             else {
                               data.use <- log1p(x = data.use)
                             }
                             return(data.use)
                           })
  avg.exp.scaled <- as.vector(x = t(x = avg.exp.scaled))
  if (split.colors) {
    avg.exp.scaled <- as.numeric(x = cut(x = avg.exp.scaled, 
                                         breaks = 20))
  }
  data.plot$avg.exp.scaled <- avg.exp.scaled
  data.plot$features.plot <- factor(x = data.plot$features.plot, 
                                    levels = features)
  data.plot$pct.exp[data.plot$pct.exp < dot.min] <- NA
  data.plot$pct.exp <- data.plot$pct.exp * 100
  if (split.colors) {
    splits.use <- vapply(X = as.character(x = data.plot$id), 
                         FUN = gsub, FUN.VALUE = character(length = 1L), pattern = paste0("^((", 
                                                                                          paste(sort(x = levels(x = object), decreasing = TRUE), 
                                                                                                collapse = "|"), ")_)"), replacement = "", 
                         USE.NAMES = FALSE)
    data.plot$colors <- mapply(FUN = function(color, value) {
      return(colorRampPalette(colors = c("grey", color))(20)[value])
    }, color = cols[splits.use], value = avg.exp.scaled)
  }
  color.by <- ifelse(test = split.colors, yes = "colors", no = "avg.exp.scaled")
  if (!is.na(x = scale.min)) {
    data.plot[data.plot$pct.exp < scale.min, "pct.exp"] <- scale.min
  }
  if (!is.na(x = scale.max)) {
    data.plot[data.plot$pct.exp > scale.max, "pct.exp"] <- scale.max
  }
  if (!is.null(x = feature.groups)) {
    data.plot$feature.groups <- factor(x = feature.groups[data.plot$features.plot], 
                                       levels = unique(x = feature.groups))
  }
  plot <- ggplot(data = data.plot, mapping = aes_string(x = "features.plot", 
                                                        y = "id")) + geom_point(mapping = aes_string(size = "pct.exp", 
                                                                                                     color = color.by)) + scale.func(range = c(0, dot.scale), 
                                                                                                                                     limits = c(scale.min, scale.max)) + theme(axis.title.x = element_blank(), 
                                                                                                                                                                               axis.title.y = element_blank()) + guides(size = guide_legend(title = "Percent Expressed")) + 
    labs(x = "Features", y = ifelse(test = is.null(x = split.by), 
                                    yes = "Identity", no = "Split Identity")) + cowplot::theme_cowplot()
  if (!is.null(x = feature.groups)) {
    plot <- plot + facet_grid(facets = ~feature.groups, scales = "free_x", 
                              space = "free_x", switch = "y") + theme(panel.spacing = unit(x = 1, 
                                                                                           units = "lines"), strip.background = element_blank())
  }
  if (split.colors) {
    plot <- plot + scale_color_identity()
  }
  else if (length(x = cols) == 1) {
    plot <- plot + scale_color_distiller(palette = cols)
  }
  else {
    plot <- plot + scale_color_gradient(low = cols[1], high = cols[2])
  }
  if (!split.colors) {
    plot <- plot + guides(color = guide_colorbar(title = "Average Expression"))
  }
  return(plot)
}
