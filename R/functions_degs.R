#' Introduce signed p-value and sort table of differentially expressed genes per performed test
#' 
#' @param degs Result table of the "Seurat::FindAllMarkers" function (requires column "cluster")
#' @param group Group results first by column(s) before sorting
#' @return Sorted table with differentially expressed genes. 
DegsSort = function(degs, group=NULL) { 
  # Introduce signed p-value score and group (if requested)
  degs$p_val_adj_score = -log10(degs$p_val_adj) * sign(degs$avg_log2FC)
  if (!is.null(group)) degs = degs %>% dplyr::group_by(dplyr::across(dplyr::all_of(group)))
  
  # Now sort table
  degs = degs %>% dplyr::arrange(-p_val_adj_score, -avg_log2FC, .by_group=TRUE) %>% as.data.frame()
  
  return(degs)
}

#' Filter table of differentially expressed genes (degs).
#' 
#' @param degs Result table of the "Seurat::FindAllMarkers" and "Seurat::FindMarkers" functions, further processed with the "DegsConvertLog" function.
#' @param cut_log2FC Log2 fold change threshold
#' @param cut_padj Adjusted p-value threshold
#' @param split_by_dir Split filtered table into a table with all degs, a table with up-regulated degs and a table down-regulated degs.
#' @return If split_by_dir is set to FALSE filtered table else list of filtered tables with all, up-regulated and down-regulated degs.
DegsFilter = function(degs, cut_log2FC, cut_padj, split_by_dir=TRUE) { 
  
  # Filter differentially expressed genes based on p-value and fold-change 
  filt = degs %>% 
    dplyr::filter(p_val_adj <= cut_padj) %>% 
    dplyr::filter(abs(avg_log2FC) >= cut_log2FC) %>% 
    as.data.frame()
  
  # Separate up- and down-regulated genes (if split_by_dir is TRUE)
  if (split_by_dir) {
    down = filt %>% 
      dplyr::filter(avg_log2FC <= -cut_log2FC) %>% 
      as.data.frame()
    up = filt %>% 
      dplyr::filter(avg_log2FC >= cut_log2FC) %>% 
      as.data.frame()
    filt = list(all=filt, up=up, down=down)
  }
  
  return(filt)
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
  avg_set[[DefaultAssay(sc)]] = c(avg_set[[DefaultAssay(sc)]],"data")
  avg_data = matrix(NA+0, nrow=length(genes), ncol=0)
  
  identities = levels(Idents(sc))
  for (as in names(avg_set)) { 
    for (sl in avg_set[[as]]) {
      if (length(genes) > 0) {
        avg_per_id = mapply(function(id) { 
          id_cells = WhichCells(sc, idents=id)
          if (sl=="data") {
            id_avg = log(Matrix::rowMeans(exp(Seurat::GetAssayData(sc[, id_cells], assay=as, slot=sl)[genes, ])))
          } else if (sl=="counts") {
            id_avg = Matrix::rowMeans(Seurat::GetAssayData(sc[, id_cells], assay=as, slot=sl)[genes,])
          }
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
    
    # Add README
    readme_table = data.frame(Column=c("p_val"), Description=c("Uncorrected p-value"))
    readme_table = rbind(readme_table, c("avg_log2FC", "Mean log2 fold change group 1 vs group2"))
    readme_table = rbind(readme_table, c("pct.1", "Fraction cells expressing gene in group 1"))
    readme_table = rbind(readme_table, c("pct.2", "Fraction cells expressing gene in group 2"))
    readme_table = rbind(readme_table, c("p_val_adj", "Adjusted p-value"))
    readme_table = rbind(readme_table, c("group", "Cluster or group"))
    readme_table = rbind(readme_table, c("gene", "Gene"))
    readme_table = rbind(readme_table, c("p_val_adj_score", "Score calculated as follows: -log10(p_val_adj)*sign(log2FC)"))
    readme_table = rbind(readme_table, c("avg_<assay>_<slot>_id<group>", "Average expression value for group; <assay> can be RNA or SCT; <slot> can be (raw) counts or (normalised) data"))
    readme_table = rbind(readme_table, c("ensembl_gene_id","Ensembl gene id (if available as annotation)"))
    readme_table = rbind(readme_table, c("external_gene_name","Gene symbol (if available as annotation)"))
    readme_table = rbind(readme_table, c("chromosome_name","Chromosome name of gene (if available as annotation)"))    
    readme_table = rbind(readme_table, c("start_position","Start position of gene (if available as annotation)"))
    readme_table = rbind(readme_table, c("end_position","End position of gene (if available as annotation)"))
    readme_table = rbind(readme_table, c("percentage_gene_gc_content","GC content of gene (if available as annotation)"))
    readme_table = rbind(readme_table, c("gene_biotype","Biotype of gene (if available as annotation)"))    
    readme_table = rbind(readme_table, c("strand","Strand of gene (if available as annotation)"))
    readme_table = rbind(readme_table, c("description","Description of gene (if available as annotation)"))
    
    degs_lst = c(list("README"=readme_table), degs_lst)
    
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

#' Returns an empty deg test table with the columns 'p_val', 'avg_log2FC', 'pct.1', 'pct.2', 'p_val_adj' and 'gene'.
#' 
#' @return A R data.frame with the columns 'p_val', 'avg_log2FC', 'pct.1', 'pct.2', 'p_val_adj' and 'gene'.
EmptyDegResultsTable = function() {
  empty_table = data.frame(p_val=as.numeric(), avg_log2FC=as.numeric(), pct.1=as.numeric(), pct.2=as.numeric(), p_val_adj=as.numeric(), gene=as.character())
  return(empty_table)
}

#' Splits a specification string. Multiple levels are specified by semicolons. Combining levels is done with the plus sign. Trims leading and trailing whitespaces.
#' 
#' @param spec_string: A specification string.
#' @first_level_separator: A fixed character string for the first string split. 
#' @second_level_separator: A fixed character string for the second string split. applied to the results of the first string split. Can be NULL. 
#' @return: A list or a list of lists.
SplitSpecificationString = function(spec_string, first_level_separator=";", second_level_separator=NULL) {
  # first split
  spec_list = trimws(unlist(strsplit(spec_string, split=first_level_separator, fixed=TRUE)))
  
  # second split (optional)
  if (!is.null(second_level_separator)) spec_list = purrr::map(spec_list, function(s) {
    unlist(purrr::map(strsplit(s, split=second_level_separator, fixed=TRUE), trimws))})
  return(spec_list)
}

#' Joins a specification list or a list of lists into a string. Sublists are joined with the plus sign and the resulting strings are joined with semicolons.
#' 
#' @param spec_list: A specification list of lists.
#' @first_level_separator: A fixed character string for the first level string join. 
#' @second_level_separator: A fixed character string for the second level string join split. Can be NULL in which case the first_level_separator will be used. 
#' @return: A specification string.
JoinSpecificationList = function(spec_list, first_level_separator=";", second_level_separator=NULL) {
  is_list_of_lists = any(purrr::map_lgl(spec_list, is.list))
  
  if (is_list_of_lists) {
    if (is.null(second_level_separator)) second_level_separator = first_level_separator
    return(paste(purrr::map(spec_list, paste, collapse=second_level_separator), collapse=first_level_separator))
  } else {
    return(paste(spec_list, collapse=first_level_separator))
  }
}

#' Tests two sets of cells for differential expression using Seurat::FindMarkers.
#' 
#' @param object A Seurat assay object or a Seurat DimReduc object.
#' @param slot If object is a Seurat assay object, which slot to use.
#' @param cells_1 The cell names in set1.
#' @param cells_2 The cell names in set2.
#' @param is_reduction Object is a Seurat DimReduc object. 
#' @param ... Additional parameters passed on to FindMarkers.
#' @return A table with the columns p_val, avg_log2FC, pct.1, pct.2, p_val_adj and gene.
DegTestCellSets = function(object, slot="data", cells_1=NULL, cells_2=NULL, is_reduction=FALSE, ...){
  # Additional arguments for FindMarkers in the three-dots construct
  additional_arguments = list(...)
  
  # Make sure that object, assay, slot, reduction, cells.1 and cells.2 are not in the additional_arguments list
  additional_arguments[["object"]] = NULL
  additional_arguments[["slot"]] = NULL
  additional_arguments[["cells.1"]] = NULL
  additional_arguments[["cells.2"]] = NULL

  # Create empty table to return when there are no results
  no_degs_results = EmptyDegResultsTable()
  
  # Check that there are at least 3 cell names and that all cell names are part of the Seurat object
  if (is.null(cells_1) || length(cells_1)<3 || any(!cells_1 %in% colnames(object))) return(no_degs_results)
  if (is.null(cells_2) || length(cells_2)<3 || any(!cells_2 %in% colnames(object))) return(no_degs_results)
  
  # Run Seurat::FindMarkers
  if (!is_reduction) {
    arguments = c(list(object=object, slot=slot, cells.1=cells_1, cells.2=cells_2), additional_arguments)
  } else {
    arguments = c(list(object=object, cells.1=cells_1, cells.2=cells_2), additional_arguments)
  }
  
  deg_results = suppressMessages((do.call(Seurat::FindMarkers, arguments)))
  if (nrow(deg_results)==0) return(no_degs_results)
  
  # Fix base 2 for log fold change
  if (!"avg_log2FC" %in% colnames(deg_results)) {
    lfc_idx = grep("avg_log\\S*FC", colnames(deg_results))
    deg_results[,lfc_idx] = deg_results[,lfc_idx] / log(2)
    col_nms = colnames(deg_results)
    col_nms[2] = "avg_log2FC"
    colnames(deg_results) = col_nms
  }
  
  # Add column gene
  deg_results$gene = rownames(deg_results)
  return(deg_results)
}

#' Given the DEG contrasts table as an R data.frame table or as an Excel file, prepares a list with DEG contrasts to do. Parses test parameter, checks and establishes defaults.
#' 
#' @param sc A Seurat single cell object.
#' @param contrasts_table A table as an R data.frame or as an Excel file. The table must at least contain the columns 'condition_column', 'condition_group1' and 'condition_group2'.
#' @param latent_vars Global variables to account for when testing. Can be overwritten by table. Can be NULL or empty.
#' @return A list with contrasts to be analysed with the DegTestCondition function. Will return an empty list if the table/file is empty/not existing. If test parameters fail, then it will an empty result table and an error message to be shown.
SetupDegContrastsList = function(sc, contrasts_table, latent_vars=NULL) {
  contrasts_list = list()
  cell_metadata = sc[[]]
  
  # If parameter is a file, read the table with the contrasts
  if (is.character(contrasts_table) && file.exists(contrasts_table)) {
    contrasts_table = openxlsx::read.xlsx(contrasts_table)
  }
  
  # If invalid or null or empty, return empty list
  if (is.null(contrasts_table) || !is.data.frame(contrasts_table) || nrow(contrasts_table)==0) return(list())
  
  # Convert into list, do checks and set up defaults
  contrasts_list = split(contrasts_table, seq(nrow(contrasts_table)))
  contrasts_list = purrr::map(seq(contrasts_list), function(i) {
    contrast = unlist(contrasts_list[[i]])
    
    error_messages = c()
    
    # deal with leading/trailing whitespaces
    contrast = purrr::map(contrast, trimws)
    
    # condition_column
    if (!contrast[["condition_column"]] %in% colnames(cell_metadata)) error_messages = c(error_messages, paste("The 'condition_column' column value '",contrast[["condition_column"]], "'  in row", i, "of the deg contrasts table is not part of the cell metadata."))
    
    # condition_group1; multiple levels can be combined with the plus sign; can be empty string to use all levels not in the condition group2 combined
    if (nchar(contrast[["condition_group1"]])==0) {
      contrast[["condition_group1"]] = list(NULL)
    } else {
      condition_group1 = SplitSpecificationString(contrast[["condition_group1"]], first_level_separator="+")
      if (any(!condition_group1 %in% unique(cell_metadata[, contrast[["condition_column"]], drop=TRUE]))) error_messages = c(error_messages, paste("At least one value of column 'condition_group1' in row",i,"of the deg contrasts table cannot be found in the column '", contrast[["condition_column"]] , "' of the cell metadata."))
      contrast[["condition_group1"]] = JoinSpecificationList(condition_group1, first_level_separator="+")
    }

    # condition_group2; multiple levels can be combined with the plus sign; can be empty string to use all levels not in the condition group1 combined (see the actual DEG functions)
    if (nchar(contrast[["condition_group2"]])==0) {
      contrast[["condition_group2"]] = list(NULL)
    } else {
      condition_group2 = SplitSpecificationString(contrast[["condition_group2"]], first_level_separator="+")
      if (any(!condition_group2 %in% unique(cell_metadata[, contrast[["condition_column"]], drop=TRUE]))) error_messages = c(error_messages, paste("At least one value of column 'condition_group2' in row",i,"of the deg contrasts table cannot be found in the column '", contrast[["condition_column"]] , "' of the cell metadata."))
      contrast[["condition_group2"]] = paste(condition_group2, collapse="+")
    }
    
    # subset_column and subset_group (if available); subsets cells based on this column prior to testing 
    if ("subset_column" %in% names(contrast)) {
      valid = TRUE
      if (!"subset_group" %in% names(contrast)) {
        error_messages = c(error_messages, paste("The 'subset_column' column must be used together with the 'subset_group' column (in row", i ,"of the deg contrasts table)."))
        valid = FALSE
      }
      
      if (valid && !contrast[["subset_column"]] %in% colnames(cell_metadata)) {
        error_messages = c(error_messages, paste("The 'subset_column' column value '",contrast[["subset_column"]], "'  in row", i, "of the deg contrasts table is not part of the cell metadata."))
        valid = FALSE
      }
      
      if (valid) {
        # subset_group; can have multiple levels separated by semicolons; can have multiple levels combined with the plus sign; can be empty string to use analyse all levels
        if (nchar(contrast[["subset_group"]])==0) {
          contrast[["subset_group"]] = JoinSpecificationList(unique(cell_metadata[, contrast[["subset_column"]], drop=TRUE]), first_level_separator=";")
        } else {
          subset_group = SplitSpecificationString(contrast[["subset_group"]], first_level_separator=";", second_level_separator="+")
          if (any(!unlist(subset_group) %in% unique(cell_metadata[, contrast[["subset_column"]], drop=TRUE]))) error_messages = c(error_messages, paste("At least one value of column 'subset_group' in row",i,"of the deg contrasts table cannot be found in the column '", contrast[["subset_column"]] , "' of the cell metadata."))
          contrast[["subset_group"]] = JoinSpecificationList(subset_group, first_level_separator=";", second_level_separator="+")
        }
      }
    }
    
    # test
    if (!"test" %in% names(contrast) || is.na(contrast[["test"]])) contrast[["test"]] = "wilcox"
    if (!contrast[["test"]] %in% c("wilcox", "bimod", "roc", "t", "negbinom", "poisson", "LR", "MAST", "DESeq2")) error_messages = c(error_messages, paste("The 'test' column value '",contrast[["test"]], "'in row", i, "of the deg contrasts table must be one of: 'wilcox', 'bimod', 'roc', 't', 'negbinom', 'poisson', 'LR', 'MAST', 'DESeq2'."))
    
    # padj
    if (!"padj" %in% names(contrast) || is.na(contrast[["padj"]])) contrast[["padj"]] = 0.05
    
    # log2FC
    if (!"log2FC" %in% names(contrast) || is.na(contrast[["log2FC"]])) contrast[["log2FC"]] = 0
    
    # min_pct
    if (!"min_pct" %in% names(contrast) || is.na(contrast[["min_pct"]])) contrast[["min_pct"]] = 0.1
    
    # assay
    if (!"assay" %in% names(contrast) || is.na(contrast[["assay"]])) contrast[["assay"]] = "RNA"
    
    if (contrast[["assay"]] %in% Assays(sc)) {
      contrast[["use_reduction"]] = FALSE
    } else if (contrast[["assay"]] %in% Reductions(sc)) {
      contrast[["use_reduction"]] = TRUE
    } else {
      c(error_messages, paste("The 'assay' column value '",contrast[["assay"]], "'in row", i, "of the deg contrasts table is neither a Seurat assay nor a Seurat reduction."))
    }
    
    # slot
    if (!"slot" %in% names(contrast) || is.na(contrast[["slot"]])) contrast[["slot"]] = "data"
    if (!contrast[["slot"]] %in% c("counts", "data", "scale.data")) c(error_messages, paste("The 'slot' column value '",contrast[["assay"]], "'in row", i, "of the deg contrasts table must be 'counts', 'data' or 'scale.data'."))
    
    # latent_vars
    if (!"latent_vars" %in% names(contrast) || is.na(contrast[["latent_vars"]])) {
      
      if (length(latent_vars)>0) contrast[["latent_vars"]] = latent_vars
      else contrast[["latent_vars"]] = NULL
    } else {
      contrast[["latent_vars"]] = SplitSpecificationString(contrast[["subset_group"]], first_level_separator=";")
    }
    
    if (!is.null(contrast[["latent_vars"]]) && length(contrast[["latent_vars"]])>0) {
      if (any(!contrast[["latent_vars"]] %in% colnames(cell_metadata))) {
        error_messages = c(error_messages, paste("At least one value of column 'latent_vars' in row",i,"of the deg contrasts table or a global latent var cannot be found cell metadata column."))
      }
    } 
    
    # add error messages
    contrast = c(contrast, list(error_messages = error_messages))
    
    # add contrast rows
    contrast[["contrast_row"]] = i
    
    return(contrast)
  })
  
  # Expand subsets so that there is now one row per subset
  contrasts_list = purrr::map(contrasts_list, function(contrast) {
    if ("subset_column" %in% names(contrast) & "subset_group" %in% names(contrast)) {
      subset_group = SplitSpecificationString(contrast[["subset_group"]], first_level_separator=";")
      contrasts_expanded = purrr::map(subset_group, function(g) {
        con = contrast
        con[["subset_group"]] = g
        return(con)
      })
    } else {
      contrasts_expanded = contrast
    }
    return(contrasts_expanded)
  })
  names(contrasts_list) = seq(contrasts_list)
  
  # Now get cell indices per contrast; this will later be used for more efficient parallelisation 
  contrasts_list = purrr::map(contrasts_list, function(contrast) {
    error_messages = c()
    
    # subset
    if ("subset_column" %in% names(contrast) && "subset_group" %in% names(contrast)) {
      subset_group_values = SplitSpecificationString(contrast[["subset_group"]], first_level_separator="+")
      is_in_subset_group = cell_metadata[,subset_column, drop=TRUE] %in% subset_group_values
    } else {
      is_in_subset_group = rep(TRUE, nrow(cell_metadata))
    }
    
    # condition_group1
    if (!is.null(contrast[["condition_group1"]])) {
      condition_group1_values = SplitSpecificationString(contrast[["condition_group1"]], first_level_separator=";")
      is_in_condition_group1 = cell_metadata[,contrast[["condition_column"]], drop=TRUE] %in% condition_group1_values
      is_in_condition_group1 = is_in_condition_group1 & is_in_subset_group
    } else {
      cells_condition_group1 = NA
    }
    
    # condition_group2
    if (!is.null(contrast[["condition_group2"]])) {
      condition_group2_values = SplitSpecificationString(contrast[["condition_group2"]], first_level_separator=";")
      is_in_condition_group2 = cell_metadata[,contrast[["condition_column"]], drop=TRUE] %in% condition_group2_values
      is_in_condition_group2 = is_in_condition_group2 & is_in_subset_group
    } else {
      cells_condition_group2 = NA
    }    
    
    # If one of the condition groups is NULL, set the cell names to the complement of the cells in the other condition group (and potential subset)
    if (is.null(contrast[["condition_group1"]]) & !is.null(contrast[["condition_group2"]])) {
      is_in_condition_group2 = !is_in_condition_group1 & is_in_subset_group
    } else if (!is.null(contrast[["condition_group1"]]) & is.null(contrast[["condition_group2"]])) {
      is_in_condition_group1 = !is_in_condition_group2 & is_in_subset_group
    }
    
    # Get cell indices and also associated latent vars
    if (sum(is_in_condition_group1)>=3) {
      contrast[["cells_group1_idx"]] = as.integer(which(is_in_condition_group1))
    } else {
      error_messages = c(error_messages,  paste("Condition group 1 in one of the contrasts specified in row", contrast[["contrast_row"]], "of the deg contrasts table has less than three cells"))
    }
    
    if (sum(is_in_condition_group2)>=3) {
      contrast[["cells_group2_idx"]] = as.integer(which(is_in_condition_group2))
    } else {
      error_messages = c(error_messages,  paste("Condition group 2 in one of the contrasts specified in row", contrast[["contrast_row"]], "of the deg contrasts table has less than three cells"))
    }
      
    contrast = c(contrast, list(error_messages = c(contrast[["error_messages"]], error_messages)))
    return(contrast)
  })

  return(contrasts_list)
}