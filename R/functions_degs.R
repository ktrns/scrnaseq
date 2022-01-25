#' Sorts table of differentially expressed genes per performed test. Introduces a signed p-value score calculated as follows:
#' p_val_adj_score = -log10(p_val_adj) * sign(avg_log2FC).
#' 
#' @param degs Result table of the "Seurat::FindAllMarkers" function or the "Seurat::FindMarkers" function.
#' @param group Group results first by column(s) before sorting.
#' @return Sorted table with differentially expressed genes. 
DegsSort = function(degs, group=NULL) { 
  # Introduce signed p-value score and group (if requested)
  degs$p_val_adj_score = -log10(degs$p_val_adj) * sign(degs$avg_log2FC)
  if (!is.null(group)) degs = degs %>% dplyr::group_by(dplyr::across(dplyr::all_of(group)))
  
  # Now sort table
  degs = degs %>% dplyr::arrange(-p_val_adj_score, -avg_log2FC, .by_group=TRUE) %>% as.data.frame()
  
  return(degs)
}

#' Filter table of differentially expressed genes.
#' 
#' @param degs Result table of the "Seurat::FindAllMarkers" or the "Seurat::FindMarkers" functions.
#' @param cut_log2FC Log2 fold change threshold.
#' @param cut_padj Adjusted p-value threshold.
#' @param split_by_dir Split filtered table into a table with all degs, a table with up-regulated degs and a table down-regulated degs.
#' @return If split_by_dir is set to FALSE filtered table else list of filtered tables with all, up-regulated and down-regulated degs.
DegsFilter = function(degs, cut_log2FC, cut_padj, split_by_dir=TRUE) { 
  
  # Filter differentially expressed genes based on p-value and fold change 
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

#' Display top marker genes (=up-regulated genes).
#' 
#' @param degs Result table of the "Seurat::FindAllMarkers" function (requires column "cluster")
#' @param n Number of top genes to show
#' @param column_1 First column to sort genes on
#' @param column_2 Second column to sort genes on
#' @return Data.frame of top genes 
DegsUpDisplayTop = function(degs, n=5, column_1="p_val_adj_score", column_2="pct.diff") { 
  
  # Calculate difference in percentage of cells that express the gene
  degs = degs %>% dplyr::mutate(pct.diff=abs(pct.1-pct.2))
  
  # Get top 5 up-regulated markers
  top = degs %>% 
    dplyr::group_by(cluster) %>% 
    dplyr::top_n(n=n, wt=get(column_1)) %>% 
    dplyr::top_n(n=n, wt=get(column_2)) %>% 
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

#' Compute average gene expression data per identity class for a set of genes.
#' 
#' @param sc Seurat object.
#' @param genes Gene list for which average data are to be extracted.
#' @return A table with average RNA counts and data per identity class for each gene.
DegsAvgDataPerIdentity = function(sc, genes, assay="RNA") { 
  # The standard average log FC is derived from assay and slot="data"
  # Add average scaled data per cluster for default assay
  avg_set = list()
  avg_set[[assay]] = "counts"
  avg_set[[DefaultAssay(sc)]] = c(avg_set[[DefaultAssay(sc)]], "data")
  avg_data = matrix(NA+0, nrow=length(genes), ncol=0)
  
  identities = levels(Idents(sc))
  for (as in names(avg_set)) { 
    for (sl in avg_set[[as]]) {
      if (length(genes) > 0) {
        avg_per_id = mapply(function(id) { 
          id_cells = WhichCells(sc, idents=id)
          if (sl=="data") {
            # This calculation is in accordance with what Seurat is doing
            id_avg = log2(Matrix::rowMeans(exp(Seurat::GetAssayData(sc[, id_cells], assay=as, slot=sl)[genes, ])) + 1)
          } else if (sl=="counts") {
            id_avg = Matrix::rowMeans(Seurat::GetAssayData(sc[, id_cells], assay=as, slot=sl)[genes, ])
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


#' Compute average gene expression data for a set of cells and a set of genes.
#' 
#' @param object Seurat assay object.
#' @param cells Cells to be used. NULL if all cells should be used.
#' @param genes Gene list for which average data are to be extracted. Can be NULL where all genes will be calculated.
#' @param slot Slot to be used (data). Can be 'counts' or 'data' or both (vector).
#' @return A table with average data for each gene.
DegsAvgData = function(object, cells=NULL, genes=NULL, slot="data") {
  if ("data" %in% slot) avg_data = as.numeric() else avg_data = NULL
  if ("counts" %in% slot) avg_counts = as.numeric() else avg_counts = NULL
  
  # If object is NULL or cells==0 or genes==0 return empty data.frame
  if (is.null(object) | (!is.null(cells) && length(cells)==0) | (!is.null(genes) && length(genes)==0)) {
    return(cbind(avg_counts, avg_data))
  }
  
  # If genes is NULL set defaults
  if (is.null(genes)) genes = rownames(object)
  
  # Subset by cells first (if requested)
  if (!is.null(cells)) object = subset(object, cells=cells)
  
  if ("data" %in% slot) {
    if (length(genes)>0) {
      # This calculation is in accordance with what Seurat is doing
      avg_data = log2(Matrix::rowMeans(exp(Seurat::GetAssayData(object, slot="data")[genes, , drop=FALSE])) + 1)
    }
  }
  
  if ("counts" %in% slot) {
    if (length(genes)>0) {
      avg_counts = Matrix::rowMeans(Seurat::GetAssayData(object, slot="data")[genes, , drop=FALSE])
    }
  }
  
  avg = as.data.frame(cbind(avg_counts, avg_data))
  if (length(genes)>0) rownames(avg) = genes
  
  return(avg)
}


#' Write differentially expressed genes or markers to an Excel file.
#' 
#' @param degs Result table of the "Seurat::FindMarkers" or "Seurat::FindAllMarkers" functions. Can also be a list of tables so that each table is written into an extra Excel tab. 
#' @param annot_ensembl Ensembl annotation for all genes with Ensembl IDs as rownames.
#' @param gene_to_ensembl Named vector for translating the Seurat gene names of the result table(s) to Ensembl IDs.
#' @param file Output file name.
#' @param additional_readme A data.frame to describe additional columns. Should contain columns 'Column' and 'Description'. Can be NULL.
#' @return Output file name.
DegsWriteToFile = function(degs, annot_ensembl, gene_to_ensembl, file, additional_readme=NULL) {
  # Always use list and add annotation
  if (is.data.frame(degs)) degs_lst = list(All=degs) else degs_lst = degs
  
  # Add Ensembl annotation
  for (i in seq(degs_lst)) {
    degs_ensembl = gene_to_ensembl[as.character(degs_lst[[i]]$gene)]
    degs_lst[[i]] = cbind(degs_lst[[i]], annot_ensembl[degs_ensembl, ])
  }
  
  # Add README
  readme_table = data.frame(Column=c("p_val"), Description=c("Uncorrected p-value"))
  readme_table = rbind(readme_table, c("avg_log2FC", "Mean log2 fold change group 1 vs group 2"))
  readme_table = rbind(readme_table, c("pct.1", "Fraction cells expressing gene in group 1"))
  readme_table = rbind(readme_table, c("pct.2", "Fraction cells expressing gene in group 2"))
  readme_table = rbind(readme_table, c("p_val_adj", "Adjusted p-value"))
  readme_table = rbind(readme_table, c("gene", "Gene"))
  readme_table = rbind(readme_table, c("p_val_adj_score", "Score calculated as follows: -log10(p_val_adj)*sign(avg_log2FC)"))
  readme_table = rbind(readme_table, c("avg.1", "Mean expression in group 1"))
  readme_table = rbind(readme_table, c("avg.2", "Mean expression in group 2"))
  readme_table = rbind(readme_table, c("ensembl_gene_id","Ensembl gene id (if available)"))
  readme_table = rbind(readme_table, c("external_gene_name","Gene symbol (if available)"))
  readme_table = rbind(readme_table, c("chromosome_name","Chromosome name of gene (if available)"))    
  readme_table = rbind(readme_table, c("start_position","Start position of gene (if available)"))
  readme_table = rbind(readme_table, c("end_position","End position of gene (if availablen)"))
  readme_table = rbind(readme_table, c("percentage_gene_gc_content","GC content of gene (if available)"))
  readme_table = rbind(readme_table, c("gene_biotype","Biotype of gene (if available)"))    
  readme_table = rbind(readme_table, c("strand","Strand of gene (if available)"))
  readme_table = rbind(readme_table, c("description","Description of gene (if available)"))
  
  if (!is.null(additional_readme)) readme_table = dplyr::bind_rows(readme_table, additional_readme)
  
  degs_lst = c(list("README"=readme_table), degs_lst)
  
  # Fix names that are more 31bp (which is too long for Excel)
  names(degs_lst) = strtrim(names(degs_lst), 31)
  
  # Output in Excel sheet
  openxlsx::write.xlsx(degs_lst, file=file)
  return(file)
}

#' Plot the number of DEGs per test.
#' 
#' @param degs Result table of the "Seurat::FindAllMarkers" function.
#' @param group Group results by column for plotting.
#' @param title Plot title.
#' @return A ggplot object.
DegsPlotNumbers = function(degs, group=NULL, title=NULL) {
  
  degs_up = degs %>% dplyr::filter(avg_log2FC > 0)
  degs_down = degs %>% dplyr::filter(avg_log2FC < 0)
  
  if ((nrow(degs_up) > 0) | (nrow(degs_down) > 0)) {
    degs_n = rbind(data.frame(degs_up, Direction=rep("Up", nrow(degs_up))), data.frame(degs_down, Direction=rep("Down", nrow(degs_down))))
    
    if (!is.null(group)) {
      degs_group_levels = unique(levels(degs_up[, group, drop=TRUE]), levels(degs_down[, group, drop=TRUE]))
      degs_n$Identity = factor(degs_n[, group, drop=TRUE], levels=degs_group_levels)
    } else {
      degs_n$Identity = as.factor("Genes")
    }
    
    degs_n = degs_n %>% dplyr::group_by(Identity, Direction) %>% dplyr::summarise(n=length(Direction))
    p = ggplot(degs_n, aes(x=Identity, y=n, fill=Direction)) + 
      geom_bar(stat="identity") +
      xlab(group) +
      AddStyle(title=title,
               fill=setNames(c("steelblue", "darkgoldenrod1"), c("Down", "Up")))
    return(p)
  }
}

#' Returns an empty deg test table with the columns 'p_val', 'avg_log2FC', 'pct.1', 'pct.2', 'p_val_adj' and 'gene'.
#' 
#' @param col_def Additional columns.
#' @return An R data.frame with the columns 'p_val', 'avg_log2FC', 'pct.1', 'pct.2', 'p_val_adj' and 'gene'.
DegsEmptyResultsTable = function() {
  empty_table = data.frame(p_val=as.numeric(), avg_log2FC=as.numeric(), pct.1=as.numeric(), pct.2=as.numeric(), p_val_adj=as.numeric(), gene=as.character())
  return(empty_table)
}

#' Returns an empty deg marker test table with the columns 'p_val', 'avg_log2FC', 'pct.1', 'pct.2', 'p_val_adj', 'cluster' and 'gene'.
#' 
#' @param clusters: Cluster names to set factor levels of empty 'cluster' column.
#' @return An R data.frame with the columns 'p_val', 'avg_log2FC', 'pct.1', 'pct.2', 'p_val_adj', 'cluster' and 'gene'.
DegsEmptyMarkerResultsTable = function(clusters) {
  empty_table = DegsEmptyResultsTable()
  empty_table$cluster = factor(as.character(), levels=clusters)
  return(empty_table[ c('p_val','avg_log2FC','pct.1','pct.2','p_val_adj','cluster','gene')])
}

#' Splits a specification string. Multiple levels are specified by semicolons. Combining levels is done with the plus sign. Trims leading and trailing whitespaces.
#' 
#' @param spec_string: A specification string.
#' @first_level_separator: A fixed character string for the first string split. 
#' @second_level_separator: A fixed character string for the second string split. Applied to the results of the first string split. Can be NULL. 
#' @return: A list or a list of lists.
SplitSpecificationString = function(spec_string, first_level_separator=";", second_level_separator=NULL) {
  # First split
  spec_list = trimws(unlist(strsplit(spec_string, split=first_level_separator, fixed=TRUE)))
  
  # Second split (optional)
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
  is_list_of_lists = any(purrr::map_lgl(spec_list, function(l) {return(is.vector(l) | is.list(l))} ))
  
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
DegsTestCellSets = function(object, slot="data", cells_1=NULL, cells_2=NULL, is_reduction=FALSE, ...){
  # Additional arguments for FindMarkers in the three-dots construct
  additional_arguments = list(...)
  
  # Make sure that object, assay, slot, reduction, cells.1 and cells.2 are not in the additional_arguments list
  additional_arguments[["object"]] = NULL
  additional_arguments[["slot"]] = NULL
  additional_arguments[["cells.1"]] = NULL
  additional_arguments[["cells.2"]] = NULL

  # Create empty table to return when there are no results
  no_degs_results = DegsEmptyResultsTable()
  
  # Check that there are at least 3 cell names and that all cell names are part of the Seurat object
  if (is.null(cells_1) || length(cells_1) < 3 || any(!cells_1 %in% colnames(object))) return(no_degs_results)
  if (is.null(cells_2) || length(cells_2) < 3 || any(!cells_2 %in% colnames(object))) return(no_degs_results)
  
  # Run Seurat::FindMarkers
  if (!is_reduction) {
    arguments = c(list(object=object, slot=slot, cells.1=cells_1, cells.2=cells_2), additional_arguments)
  } else {
    arguments = c(list(object=object, cells.1=cells_1, cells.2=cells_2), additional_arguments)
  }
  
  deg_results = suppressMessages(do.call(Seurat::FindMarkers, arguments))
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
DegsSetupContrastsList = function(sc, contrasts_table, latent_vars=NULL) {
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
    # Unlist does not keep the classes of the variables
    # Hence we do a manual unlisting
    contrast = list()
    for (j in colnames(contrasts_list[[i]])) {
      contrast[[j]] = contrasts_list[[i]][, j, drop=TRUE]
      if (class(contrast[[j]]) == "character") contrast[[j]] = trimws(contrast[[j]])
    }

    error_messages = c()
    
    # Get condition_column
    if (contrast[["condition_column"]] %in% colnames(cell_metadata)) {
      
      # Get condition_group1; multiple levels can be combined with the plus sign; can be empty string to use all levels not in the condition group2 combined
      if (nchar(contrast[["condition_group1"]])==0) {
        contrast[["condition_group1"]] = as.character(NA)
      } else {
        
        condition_group1 = SplitSpecificationString(contrast[["condition_group1"]], first_level_separator="+")
        if (any(!condition_group1 %in% unique(cell_metadata[, contrast[["condition_column"]], drop=TRUE]))) {
          error_messages = c(error_messages, paste("At least one value of column 'condition_group1' in row",i,"of the deg contrasts table cannot be found in the column '", contrast[["condition_column"]] , "' of the cell metadata."))
        }
        contrast[["condition_group1"]] = JoinSpecificationList(condition_group1, first_level_separator="+")
        
      }
      
      
      # Get condition_group2; multiple levels can be combined with the plus sign; can be empty string to use all levels not in the condition group1 combined (see the actual DEG functions)
      if (nchar(contrast[["condition_group2"]])==0) {
        contrast[["condition_group2"]] = as.character(NA)
      } else {
        
        condition_group2 = SplitSpecificationString(contrast[["condition_group2"]], first_level_separator="+")
        if (any(!condition_group2 %in% unique(cell_metadata[, contrast[["condition_column"]], drop=TRUE]))) {
          error_messages = c(error_messages, paste("At least one value of column 'condition_group2' in row",i,"of the deg contrasts table cannot be found in the column '", contrast[["condition_column"]] , "' of the cell metadata."))
        }
        contrast[["condition_group2"]] = paste(condition_group2, collapse="+")
        
      }
    } else {
      
      # Error if condition_column not column of cell metadata
      error_messages = c(error_messages, paste0("The 'condition_column' column value '",contrast[["condition_column"]], "'  in row ", i, " of the deg contrasts table is not part of the cell metadata."))
    }

    # Get subset_column and subset_group (if available); subsets cells based on this column prior to testing
    if (!"subset_column" %in% names(contrast)) contrast[["subset_column"]] = as.character(NA)
    if (!"subset_group" %in% names(contrast)) contrast[["subset_group"]] = as.character(NA)
    
    if (!is.na(contrast[["subset_column"]])) {
      valid = TRUE
      
      # Error when subset_column but not subset_group specified
      if (is.na(contrast[["subset_column"]])) {
        error_messages = c(error_messages, paste("The 'subset_column' column must be used together with the 'subset_group' column (in row", i ,"of the deg contrasts table)."))
        valid = FALSE
      }
      
      # Error when subset_column is not a column of cell metadata
      if (valid && !contrast[["subset_column"]] %in% colnames(cell_metadata)) {
        error_messages = c(error_messages, paste0("The 'subset_column' column value '",contrast[["subset_column"]], "' in row ", i, " of the deg contrasts table is not part of the cell metadata."))
        valid = FALSE
      }
      
      if (valid) {
        if (nchar(contrast[["subset_group"]])==0) {
          # If subset_group is empty: get all levels or unique values and join them by semicolon 
          subset_group_values = levels(cell_metadata[, contrast[["subset_column"]], drop=TRUE])
          if (length(subset_group_values)==0) subset_group_values = unique(sort(cell_metadata[, contrast[["subset_column"]], drop=TRUE])) 
          contrast[["subset_group"]] = JoinSpecificationList(subset_group_values, first_level_separator=";")
        } else {
          # If subset_group is not empty: split by semicolon, then by plus; then validate values 
          subset_group = SplitSpecificationString(contrast[["subset_group"]], first_level_separator=";", second_level_separator="+")
          if (any(!unlist(subset_group) %in% unique(cell_metadata[, contrast[["subset_column"]], drop=TRUE]))) {
            error_messages = c(error_messages, paste("At least one value of column 'subset_group' in row",i,"of the deg contrasts table cannot be found in the column '", contrast[["subset_column"]] , "' of the cell metadata."))
          }
          contrast[["subset_group"]] = JoinSpecificationList(subset_group, first_level_separator=";", second_level_separator="+")
        }
      }
    }
    
    # Get test method
    if (!"test" %in% names(contrast) || is.na(contrast[["test"]])) contrast[["test"]] = "wilcox"
    if (!contrast[["test"]] %in% c("wilcox", "bimod", "roc", "t", "negbinom", "poisson", "LR", "MAST", "DESeq2")) error_messages = c(error_messages, paste0("The 'test' column value '",contrast[["test"]], "' in row ", i, " of the deg contrasts table must be one of: 'wilcox', 'bimod', 'roc', 't', 'negbinom', 'poisson', 'LR', 'MAST', 'DESeq2'."))
    
    # Get padj
    if (!"padj" %in% names(contrast) || is.na(contrast[["padj"]])) contrast[["padj"]] = 0.05
    contrast[["padj"]] = as.numeric(contrast[["padj"]])
    
    # Get log2FC
    if (!"log2FC" %in% names(contrast) || is.na(contrast[["log2FC"]])) contrast[["log2FC"]] = 0
    contrast[["log2FC"]] = as.numeric(contrast[["log2FC"]])
    
    # Get min_pct
    if (!"min_pct" %in% names(contrast) || is.na(contrast[["min_pct"]])) contrast[["min_pct"]] = 0.1
    contrast[["min_pct"]] = as.numeric(contrast[["min_pct"]])
    
    # Get assay
    if (!"assay" %in% names(contrast) || is.na(contrast[["assay"]])) contrast[["assay"]] = "RNA"
    
    if (contrast[["assay"]] %in% Assays(sc)) {
      contrast[["use_reduction"]] = FALSE
    } else if (contrast[["assay"]] %in% Reductions(sc)) {
      contrast[["use_reduction"]] = TRUE
    } else {
      c(error_messages, paste0("The 'assay' column value '",contrast[["assay"]], "' in row ", i, " of the deg contrasts table is neither a Seurat assay nor a Seurat reduction."))
    }
    
    # Get slot
    if (!"slot" %in% names(contrast) || is.na(contrast[["slot"]])) contrast[["slot"]] = "data"
    if (!contrast[["slot"]] %in% c("counts", "data", "scale.data")) c(error_messages, paste0("The 'slot' column value '",contrast[["assay"]], "' in row ", i, " of the deg contrasts table must be 'counts', 'data' or 'scale.data'."))
    
    # Get latent_vars
    if (!"latent_vars" %in% names(contrast) || is.na(contrast[["latent_vars"]])) {
      contrast[["latent_vars"]] = NULL
    } else {
      contrast[["latent_vars"]] = SplitSpecificationString(contrast[["latent_vars"]], first_level_separator=";")
    }
    
    if (!is.null(contrast[["latent_vars"]]) && length(contrast[["latent_vars"]]) > 0) {
      if (any(!contrast[["latent_vars"]] %in% colnames(cell_metadata))) {
        error_messages = c(error_messages, paste("At least one value of column 'latent_vars' in row",i,"of the deg contrasts table or a global latent var cannot be found cell metadata column."))
      }
    } 
    
    # Get number of cells to downsample
    if (!"downsample_cells_n" %in% names(contrast) || is.na(contrast[["downsample_cells_n"]])) {
      contrast[["downsample_cells_n"]] = NULL
    } else {
      contrast[["downsample_cells_n"]] = as.integer(contrast[["downsample_cells_n"]])
    }
    
    # Add error messages
    contrast = c(contrast, list(error_messages = error_messages))
    
    # Add contrast row number
    contrast[["contrast_row"]] = i
    
    return(contrast)
  })
  
  
  # Expand subsets so that there is now one row per subset
  contrasts_list = purrr::map(contrasts_list, function(contrast) {
    if (!is.na(contrast[["subset_column"]]) && !is.na(contrast[["subset_group"]])) {
      subset_group = SplitSpecificationString(contrast[["subset_group"]], first_level_separator=";")
      contrasts_expanded = purrr::map(subset_group, function(g) {
        con = contrast
        con[["subset_group"]] = g
        return(con)
      })
    } else {
      contrasts_expanded = list(contrast)
    }
    return(contrasts_expanded)
  })
  contrasts_list = purrr::flatten(contrasts_list)
  names(contrasts_list) = seq(contrasts_list)
  
  
  # Now get cell indices per contrast; this will later be used for more efficient parallelisation 
  contrasts_list = purrr::map(contrasts_list, function(contrast) {
    error_messages = c()
    
    # If there were already errors, just return
    if (length(contrast[["error_messages"]]) > 0) return(c(contrast, list(cells_group1_idx=as.integer(), cells_group2_idx=as.integer())))
    
    # Do subset
    if (!is.na(contrast[["subset_column"]]) && !is.na(contrast[["subset_group"]])) {
      subset_group_values = SplitSpecificationString(contrast[["subset_group"]], first_level_separator="+")
      is_in_subset_group = cell_metadata[, contrast[["subset_column"]], drop=TRUE] %in% subset_group_values
    } else {
      is_in_subset_group = rep(TRUE, nrow(cell_metadata))
    }
    
    # Do condition_group1
    if (!is.na(contrast[["condition_group1"]])) {
      condition_group1_values = SplitSpecificationString(contrast[["condition_group1"]], first_level_separator="+")
      is_in_condition_group1 = cell_metadata[, contrast[["condition_column"]], drop=TRUE] %in% condition_group1_values
      is_in_condition_group1 = is_in_condition_group1 & is_in_subset_group
    } else {
      cells_condition_group1 = as.character(NA)
    }
    
    # Do condition_group2
    if (!is.na(contrast[["condition_group2"]])) {
      condition_group2_values = SplitSpecificationString(contrast[["condition_group2"]], first_level_separator="+")
      is_in_condition_group2 = cell_metadata[, contrast[["condition_column"]], drop=TRUE] %in% condition_group2_values
      is_in_condition_group2 = is_in_condition_group2 & is_in_subset_group
    } else {
      cells_condition_group2 = as.character(NA)
    }    
    
    # If one of the condition groups is NA, set the cell names to the complement of the cells in the other condition group (and potential subset)
    if (is.na(contrast[["condition_group1"]]) && !is.na(contrast[["condition_group2"]])) {
      is_in_condition_group2 = !is_in_condition_group1 & is_in_subset_group
    } else if (!is.na(contrast[["condition_group1"]]) && is.na(contrast[["condition_group2"]])) {
      is_in_condition_group1 = !is_in_condition_group2 & is_in_subset_group
    }
    
    # If both condition groups are NA, return with an error
    if (is.na(contrast[["condition_group1"]]) && is.na(contrast[["condition_group2"]])) {
      error_messages = c(error_messages,  paste("Only one condition group can be empty (to specify the complement of the other group) in row", contrast[["contrast_row"]], "of the deg contrasts table."))
    }
    
    # Get cell indices and also associated latent vars
    contrast[["cells_group1_idx"]] = as.integer(which(is_in_condition_group1))
    if (length(contrast[["cells_group1_idx"]]) < 3) {
      error_messages = c(error_messages,  paste("Condition group 1 in one of the contrasts specified in row", contrast[["contrast_row"]], "of the deg contrasts table has less than three cells."))
    }
    
    contrast[["cells_group2_idx"]] = as.integer(which(is_in_condition_group2))
    if (length(contrast[["cells_group2_idx"]]) < 3) {
      error_messages = c(error_messages,  paste("Condition group 2 in one of the contrasts specified in row", contrast[["contrast_row"]], "of the deg contrasts table has less than three cells."))
    }
    
    # Downsample groups if requested
    if (!is.null(contrast[["downsample_cells_n"]])) {
      contrast[["cells_group1_idx"]] = sample(x=contrast[["cells_group1_idx"]], 
                                              size=min(contrast[["downsample_cells_n"]], length(contrast[["cells_group1_idx"]])))
      contrast[["cells_group2_idx"]] = sample(x=contrast[["cells_group2_idx"]], 
                                              size=min(contrast[["downsample_cells_n"]], length(contrast[["cells_group2_idx"]])))
    }
      
    contrast = c(contrast, list(error_messages = c(contrast[["error_messages"]], error_messages)))
    return(contrast)
  })

  return(contrasts_list)
}

#' Returns an empty Enrichr results table.
# '
# '
#' @param overlap_split: If TRUE, then table will contain the two columns 'In.List' and 'In.Annotation' (which result from splitting 'Overlap') instead of the column 'Overlap'.
#' @return An empty Enrichr results dataframe.
EmptyEnrichrDf = function(overlap_split=FALSE) {
  if (overlap_split) {
    return(data.frame(Term=as.character(), In.List=as.numeric(), In.Annotation=as.numeric(), P.value=as.numeric(), Adjusted.P.value=as.numeric(), Odds.Ratio=as.numeric(), Combined.Score=as.numeric(), Genes=as.character())) 
  } else {
    return(data.frame(Term=as.character(), Overlap=as.character(), P.value=as.numeric(), Adjusted.P.value=as.numeric(), Odds.Ratio=as.numeric(), Combined.Score=as.numeric(), Genes=as.character())) 
  }
}

#' Tests a list of entrez gene symbols for functional enrichment via Enrichr.
# '
#' @param genes: A vector of entrez gene symbols.
#' @param databases: A vector of Enrichr databases with functional annotation.
#' @param padj: Maximum adjusted p-value (0.05).
#' @return The path to the data directory.
EnrichrTest = function(genes, databases, padj=0.05) {
  empty_enrichr_df = EmptyEnrichrDf()
  
  # Run enrichr
  if (length(genes) >= 3 & length(databases) > 0) {
    enrichr_results = suppressMessages(enrichR::enrichr(unique(unname(genes)), databases=databases))
  } else {
    enrichr_results = purrr::map(databases, function(d) return(empty_enrichr_df)) # no good gene lists
    names(enrichr_results) = databases
  }
  
  databases = setNames(databases, databases)
  enrichr_results = purrr::map(databases, function(d) {
    if (nrow(enrichr_results[[d]]) == 0) { # table is empty
      return(empty_enrichr_df)
    } else if (ncol(enrichr_results[[d]]) == 1) { # there was an actual error
      warning(paste("Enrichr returns with error: ", enrichr_results[[d]][1, 1]))
      return(empty_enrichr_df)
    } else {
      return(enrichr_results[[d]])
    }
  })
  
  # Drop Old.P.value and Old.Adjusted.P.value columns if they exist
  enrichr_results = purrr::map(enrichr_results, function(df) {
    df[["Old.P.value"]] = NULL
    df[["Old.Adjusted.P.value"]] = NULL
    return(df)
  })
  
  # Some enrichr server do not return all columns - add missing columns and set them to NA
  enrichr_results = purrr::map(enrichr_results, function(df) {
    character_cols = c("Term", "Overlap", "Genes")
    missing = !character_cols %in% colnames(df)
    if (any(missing)) {
      for (c in character_cols[missing]) df[[c]] = as.character(NA)
    }
    
    numeric_cols = c("P.value", "Adjusted.P.value", "Odds.Ratio", "Combined.Score")
    missing = !numeric_cols %in% colnames(df)
    if (any(missing)) {
      for (c in numeric_cols[missing]) df[[c]] = as.numeric(NA)
    }
    
    df = df %>% dplyr::select(Term, Overlap, P.value, Adjusted.P.value, Odds.Ratio, Combined.Score, Genes)
    return(df)
  })
  
  # Filter by adjusted p-value
  enrichr_results = purrr::map(enrichr_results, dplyr::filter, Adjusted.P.value < padj)
  
  # Split column "Overlap" into numbers
  enrichr_results = purrr::map(enrichr_results, tidyr::separate, col=Overlap, into=c("In.List", "In.Annotation"), sep="/", convert=TRUE)
  
  # Remap entrez genes to seurat rownames (if seurat rownames are the names of the genes vector)
  if (!is.null(names(genes))) {
    value2names = split(names(genes), genes)
    
    enrichr_results = purrr::map(enrichr_results, function(df) {
      values_mapped_to_names = purrr::map(strsplit(df$Genes, split=";", fixed=TRUE), function(g) {
        return(paste(unlist(value2names[g]), collapse=";")) 
      })
      df$Genes.Species = values_mapped_to_names
      return(df)
    })
  }
  
  return(enrichr_results)
}

#' Writes the enrichr results to an Excel file(s). Includes a README.
# '
#' @param enrichr_results: A list with enrichr results.
#' @param file: A file path.
#' @return The path to the file.
EnrichrWriteResults = function(enrichr_results, file) {
  
  # README
  readme_table = data.frame(Column=c("Term"), Description=c("Functional annotation in database"))
  readme_table = rbind(readme_table, c("In.List", "Number of genes in list of interest with this functional annotation"))
  readme_table = rbind(readme_table, c("In.Annotation", "Total number of genes with this functional annotation"))
  readme_table = rbind(readme_table, c("P.value", "P-value (uncorrected)"))
  readme_table = rbind(readme_table, c("Adjusted.P.value", "P-value (adjusted for multiple testing), use this one"))
  readme_table = rbind(readme_table, c("Odds.Ratio", "How to interpret whether or not the gene list is just random (<1 less than by chance, =1 equals chance, >1 more than by chance)"))
  readme_table = rbind(readme_table, c("Combined.Score", "Combined enrichr score"))
  readme_table = rbind(readme_table, c("Genes", "Enrichr genes in functional annotation"))
  readme_table = rbind(readme_table, c("Genes.Species", "Species genes in functional annotation (which are translated to Enrichr genes)"))
  enrichr_results = c(list("README"=readme_table), enrichr_results)
  
  # Fix names that are more 31bp (which is too long for Excel)
  names(enrichr_results) = strtrim(names(enrichr_results), 31)
  
  # Output in Excel sheet
  openxlsx::write.xlsx(enrichr_results, file=file)
  
  return(file)
}


#' Flattens the list of Enrichr results into one data.frame.
# '
#' @param enrichr_results: A list with enrichr results.
#' @return A data.frame with the columns reported by Enrichr as well as the db
FlattenEnrichr = function(enrichr_results) {
  enrichr_results_flat = purrr::map(names(enrichr_results), function(n) {
    return(data.frame(Database=rep(n, nrow(enrichr_results[[n]])), enrichr_results[[n]]))
  })
  enrichr_results_flat = purrr::invoke(dplyr::bind_rows, enrichr_results_flat)
  return(enrichr_results_flat)
}
