#' Reads one or more counts tables (e.g. provided by SmartSeq-2) and converts them into Seurat objects.
#' 
#' @param path Path to a counts table. Cell metadata can be passed by a file metadata.tsv.gz which must be in the same directory and where the first column is the cell name.
#' @param project A project name for the dataset ("SeuratProject").
#' @param project A project name for the dataset.
#' @param row_name_column Name or index of the column which should be used as row names in the Seurat object (2).
#' @param convert_row_names Named vector for converting the row names (e.g. from Ensembl ids to gene symbols). Does not need to contain all row names. Can be NULL in which case row names are not converted.
#' @param feature_type_column Name or index of the column which should be used for feature types (NULL). Can be NULL in which case feature types are "Gene expression" and "ERCC".
#' @param columns_to_drop Names or indices of the columns which should be dropped (excluding feature_name_column or feature_type_column).
#' @param drop_non_numeric_columns After dropping columns manually with columns_to_drop: Should any remaining non-numeric columns (excluding feature_name_column or feature_type_column) be droppped automatically? If TRUE, they will be dropped silently. If FALSE, then an error message will be printed. 
#' @param feature_type_to_assay_name How should the assays for the different feature type be named? Default is: "Gene expression" = "RNA","Antibody Capture" = "ADT","CRISPR Guide Capture" = "Crispr", "ERCC" = "ERCC" and "Custom" = "Custom". Also sets the order in which the assays are loaded.
#' @param col_sep Separator used in the counts table (tab).
#' @param parse_plate_information If the cell names contain plate information: parse information and add to metadata.
#' @param plate_information_regex If the cell names contain plate information: a regular expression pattern with capture groups for sample name, plate number, row or column which must match the entire cell name ("^(\\S+)_(\\d+)_([A-Z])(\\d+)$").
#' @param sample_name_group If the cell names contain plate information: index of the capture group which contains the sample name of the cells ("1"). Can be NULL in which case it will not be evaluated and the sample name will be NA.
#' @param plate_number_group If the cell names contain plate information: index of the capture group which contains the plate number name ("2"). Can be NULL in which case it will not be evaluated and the plate number will be NA.
#' @param row_name_group If the cell names contain plate information: index of the capture group which contains the plate row name ("3"). Can be NULL in which case it will not be evaluated and the plate row name will be NA.
#' @param col_name_group If the cell names contain plate information: index of the capture group which contains the plate col name ("4"). Can be NULL in which case it will not be evaluated and the plate col name will be NA.
#' @param return_samples_as_datasets If the cell names contain plate information: return each of the parsed samples as separate dataset ("TRUE").
#' @param cellnames_suffix If not NULL, add this string as suffix to the cell names.
#' @return A list of Seurat objects.
ReadCountsTable = function(counts_table, project="SeuratProject", row_name_column=1, convert_row_names=NULL, feature_type_column=NULL, columns_to_drop=NULL, drop_non_numeric_columns=TRUE, feature_type_to_assay_name=NULL, sep="\t", parse_plate_information=TRUE, plate_information_regex='^(\\S+)_(\\d+)_([A-Z])(\\d+)$', sample_name_group=1, plate_number_group=2, row_name_group=3, col_name_group=4, return_samples_as_datasets=TRUE, cellnames_suffix=NULL) {
  library(magrittr)
  
  # Defaults
  if (is.null(feature_type_to_assay_name)) feature_type_to_assay_name = c("Gene Expression"="RNA", 
                                                                          "Antibody Capture"="ADT",
                                                                          "Multiplexing Capture"="HTO",
                                                                          "CRISPR Guide Capture"="Crispr", 
                                                                          "Custom"="Custom", 
                                                                          "ERCC"="ERCC")
  
  # Setting this is neccessary for big tables
  Sys.setenv("VROOM_CONNECTION_SIZE" = 131072*10)
  
  # Read counts data
  if (!file.exists(counts_table)) stop(sprintf("ReadCountsTable: The counts file %s does not exist!",counts_table))
  feature_data = readr::read_delim(counts_table, delim=sep, col_names=TRUE, comment="#", progress=FALSE, col_types = readr::cols())
  
  # Variables row_name_column, feature_type_column and columns_to_drop: if (numeric) index, convert to column name (easier)
  if (!is.null(row_name_column) & is.numeric(row_name_column)) row_name_column = colnames(feature_data)[row_name_column]
  if (!is.null(feature_type_column) & is.numeric(feature_type_column)) feature_type_column = colnames(feature_data)[feature_type_column]
  if (!is.null(columns_to_drop) & length(columns_to_drop) > 0 & is.numeric(columns_to_drop)) columns_to_drop = colnames(feature_data)[columns_to_drop]
  
  # Drop unwanted columns
  if (!is.null(columns_to_drop) & length(columns_to_drop) > 0) {
    feature_data = feature_data[, !colnames(feature_data) %in% columns_to_drop]
  }
  
  # Check if there are other non-numeric columns
  invalid = unlist(lapply(colnames(feature_data), function(n) {
    !(n %in% c(row_name_column, feature_type_column) | is.numeric(feature_data[, n, drop=TRUE]))
  }))
  if (sum(invalid) > 0) {
    if (drop_non_numeric_columns) {
      # If yes and drop_non_numeric_columns: just remove them
      feature_data = feature_data[, !invalid]
    } else {
      # Otherwise print an error
      invalid = invalid[which(invalid)]
      stop(sprintf("ReadCountsTable: Some columns in the counts table do not contain numeric data: %s!", 
                                first_n_elements_to_string(names(invalid))))
    }
  }
  
  # Create a data frame for feature id, feature name and feature type (similar to features.tsv in 10x datasets)
  # Variable feature_type_column: indicates the type of feature; if null, assume, there is only "Gene Expression" features; 
  #   Only check for ERCC controls which will be separate
  features_ids_types = data.frame(feature_id = feature_data[, row_name_column, drop=TRUE], 
                                  feature_name = feature_data[, row_name_column, drop=TRUE],
                                  stringsAsFactors=FALSE)
  feature_data[, row_name_column] = NULL
  if (is.null(feature_type_column)) {
    features_ids_types$feature_type = ifelse(grepl(pattern="^ERCC-", x=features_ids_types$feature_name), "ERCC", "Gene Expression")
  } else{
    features_ids_types$feature_type = feature_data[, feature_type_column]
    feature_data[, feature_type_column] = NULL
  }
  features_ids_types[, "feature_id"] = make.unique(features_ids_types[, "feature_id"])
  
  # Define a named vector for the row names: its names are the feature names in the dataset and its values are the final names in the Seurat object
  seurat_row_names = setNames(features_ids_types[, "feature_id"],features_ids_types[, "feature_id"])
  if (!is.null(convert_row_names) & length(convert_row_names) > 0) {
    convert_ids = which(seurat_row_names %in% names(convert_row_names))
    seurat_row_names[convert_ids] = convert_row_names[seurat_row_names[convert_ids]]
  }
  seurat_row_names = setNames(make.unique(gsub(pattern="_", replacement="-", x=seurat_row_names, fixed=TRUE)),names(seurat_row_names))
  rownames(features_ids_types) = unname(seurat_row_names)
  
  # Split by feature type and add Seurat-compatible row names
  feature_data = as.data.frame(feature_data)
  rownames(feature_data) = rownames(features_ids_types)
  feature_data = split(feature_data, as.factor(features_ids_types$feature_type))
  
  # Feature type to assay name
  feature_types = names(feature_data)
  missed = feature_types[!feature_types %in% names(feature_type_to_assay_name)]
  if (length(missed) > 0) stop(sprintf("ReadCountsTable: The 'feature_type_to_assay_name' argument misses some feature types: %s!",
                                     first_n_elements_to_string(missed)))
  
  # Sort feature types by their order in feature_type_to_assay_name
  feature_types = feature_types[order(match(feature_types, names(feature_type_to_assay_name)))]
  
  # Read metadata file (if available)
  metadata_file = file.path(dirname(path), "metadata.tsv.gz")
  if (file.exists(metadata_file)) {
    metadata_table = readr::read_delim(metadata_file, delim="\t", col_names=TRUE, comment="#", progress=FALSE, col_types = readr::cols())
    metadata_table = as.data.frame(metadata_table)
    
    # Check that it contains all cell names
    barcodes = colnames(feature_data[[1]])
    missed = barcodes[!barcodes %in% metadata_table[, 1, drop=TRUE]]
    if (length(missed) > 0) stop(sprintf("ReadCountsTable: The 'metadata.tsv.gz' file misses some cell names: %s!", 
                                                    first_n_elements_to_string(missed)))
    rownames(metadata_table) = metadata_table[, 1, drop=TRUE]
    
    # Some column names are not allowed since they would be overwritten later
    colnames_metadata_table = colnames(metadata_table)
    invalid = colnames_metadata_table[colnames_metadata_table %in% c("orig.ident", "nCount_RNA", "nFeature_RNA", "nCount_SCT", "nFeature_SCT", "percent_mt", "percent_ercc", "S.Score", "G2M.Score", "Phase", "CC.Difference", "SampleName", "PlateNumber", "PlateRow", "PlateCol", "seurat_clusters")]
    if (length(invalid)>0) stop(sprintf("ReadCountsTable: Some column names in 'metadata.tsv.gz' are not allowed since they would be overwritten: %s!", 
                                        first_n_elements_to_string(invalid)))
    
  } else {
    metadata_table = data.frame(Cells=colnames(feature_data[[1]]) ,stringsAsFactors=FALSE)
    rownames(metadata_table) = metadata_table[, 1, drop=TRUE]
  }
  metadata_table$orig.dataset = factor(project)
  
  # If a regex for parsing the plate information from the names has been provided, parse: sample name, plate number, plate row and plate column
  if (parse_plate_information) {
    plate_information = ParsePlateInformation(colnames(feature_data[[1]]), pattern=plate_information_regex, sample_name_group=sample_name_group, plate_number_group=plate_number_group, row_name_group=row_name_group, col_name_group=col_name_group)
    
    # If the pattern did not match for names, warn; then set plate_information$SampleName to the unparsed name
    invalid = rownames(plate_information[is.na(plate_information$SampleName), ])
    if (any(is.na(plate_information$SampleName))) {
      stop(sprintf("ReadCountsTable: Could not parse plate information from the following cell names: %s. Check the regular expression ('plate_information_regex')!", first_n_elements_to_string(invalid)))
    }
    
    # Then add to metadata
    metadata_table = merge(metadata_table,plate_information,by="row.names", all.x=TRUE)
    rownames(metadata_table) = metadata_table$Row.names
    metadata_table$Row.names = NULL
  }
  
  # Group cells by samples for datasets
  samples_to_process = list()
  if (return_samples_as_datasets) {
    if (parse_plate_information) {
      samples_to_process = split(rownames(metadata_table), metadata_table$SampleName)
    } else {
      warning("ReadCountsTable: The option 'return_samples_as_datasets' can only work when the sample name is parsed from the cell name.")
    }
  } else {
    samples_to_process[[project]] = rownames(metadata_table)
  }
  
  # Now create Seurat objects
  sc = list()
  for (s in names(samples_to_process)) {
    
    # New name
    #n = paste(project, s, sep=".")
    n = s

    # Create Seurat object with first assay
    # Include cell and feature metadata
    c = samples_to_process[[s]] # cells to include for sample
    f = feature_types[1] # feature type
    a = unname(feature_type_to_assay_name[f]) # assay name
    d = as(as.matrix(feature_data[[f]][,c, drop=FALSE]), "dgCMatrix") # cell data
    m = metadata_table[c, -1, drop=FALSE] # metadata for cells, also drop first column which contains cell names
    if (ncol(m)==0 | nrow(m)==0) m=NULL
    sc[[n]] = Seurat::CreateSeuratObject(counts=d, 
                                  project=n, 
                                  assay=a, 
                                  min.cells=0, 
                                  min.features=0, 
                                  names.delim=NULL, 
                                  names.field=NULL, 
                                  meta.data=m)
    
    # Check that the feature symbols are the same and add feature meta information
    nms = rownames(sc[[n]][[a]])
    missed = nms[!nms %in% rownames(features_ids_types)]
    if (length(missed) > 0) stop(sprintf("ReadCountsTable: The 'CreateSeuratObject' method modifies feature symbols for assay %s not as expected: %s!", 
                                                  a, first_n_elements_to_string(missed)))
    sc[[n]][[a]] = Seurat::AddMetaData(sc[[n]][[a]], features_ids_types[rownames(sc[[n]][[a]]), ])
    
    # Now add remaining assays
    for (f in feature_types[-1]) {
      a = unname(feature_type_to_assay_name[f])
      d = as(as.matrix(feature_data[[f]][,c, drop=FALSE]), "dgCMatrix")
      sc[[n]][[a]] = Seurat::CreateAssayObject(counts=d,
                                          min.cells=0,
                                          min.features=0)
      nms = rownames(sc[[n]][[a]])
      missed = nms[!nms %in% rownames(features_ids_types)]
      if (length(missed) > 0) stop(sprintf("ReadCountsTable: The 'CreateAssayObject' method modifies feature symbols for assay %s not as expected: %s!", 
                                                    a, first_n_elements_to_string(missed)))
      sc[[n]][[a]] = Seurat::AddMetaData(sc[[n]][[a]],features_ids_types[rownames(sc[[n]][[a]]), ])
    }
    
    # Add suffixes if requested
    if (!is.null(cellnames_suffix)) {
      sc[[n]] =  Seurat::RenameCells(sc[[n]], new.names=paste0(colnames(sc[[n]]), cellnames_suffix))
    }
  }
  
  return(sc)
}

#' Reads a sparse matrix dataset (e.g. provided by 10X) and converts it into a Seurat object.
#' 
#' @param path Path to the directory with sparse matrix dataset. Must contain matrix.mtx.gz, features.tsv.gz and barcodes.tsv.gz. Can contain additional metadata in a file metadata.tsv.gz where the first column is the cell name.
#' @param project A project name for the dataset.
#' @param row_name_column Name or index of the column in features.tsv.gz which should be used as row names in the Seurat object (2).
#' @param convert_row_names Named vector for converting the row names obtained from features.tsv.gz (e.g. from Ensembl ids to gene symbols). Does not need to contain all row names. Can be NULL in which case row names are not converted.
#' @param feature_type_column Name or index of the column which should be used for feature types (3).
#' @param feature_type_to_assay_name How should the assays for the different feature type be named? Default is: "Gene expression" = "RNA","Antibody Capture" = "ADT","CRISPR Guide Capture" = "Crispr", "ERCC" = "ERCC" and "Custom" = "Custom". Also sets the order in which the assays are loaded.
#' @param hto_names If Antibody Capture data are used as hashtags for multiplexing, a vector with feature names to be used as hashtags. If a named vector is provided, the hashtags will be renamed accordingly.
#' @param hto_regex If Antibody Capture data are used as hashtags for multiplexing, a regular expression to identify the hashtag (e.g. HashTag).
#' @param return_samples_as_datasets If there are multiple samples in the dataset: return each as separate dataset ("TRUE").
#' @param cellnames_suffix If not NULL, add this string as suffix to the cell names. Note that 10x-style suffixes (e.g. ACGTACGCTCAGT-1) will be removed.
#' @return A Seurat object.
ReadSparseMatrix = function(path, project="SeuratProject", row_name_column=2, convert_row_names=NULL, feature_type_column=3, feature_type_to_assay_name=NULL, hto_names=NULL, hto_regex=NULL, cellnames_suffix=NULL) {
  library(magrittr)
  
  # Defaults
  if (is.null(feature_type_to_assay_name)) feature_type_to_assay_name = c("Gene Expression"="RNA",
                                                                          "Antibody Capture"="ADT",
                                                                          "Multiplexing Capture"="HTO",
                                                                          "CRISPR Guide Capture"="Crispr",
                                                                          "Custom"="Custom",
                                                                          "ERCC"="ERCC")
  
  # Check that files exist and read barcodes.tsv.gz and features.tsv.gz separately
  if (!file.exists(file.path(path, "barcodes.tsv.gz"))) stop(sprintf("ReadSparseMatrix: Could not find file 'barcodes.tsv.gz' at directory '%s'!", path))
  barcodes = read.delim(file.path(path, "barcodes.tsv.gz"), header=FALSE, stringsAsFactors=FALSE)[, 1]
  
  if (!file.exists(file.path(path, "features.tsv.gz"))) stop(sprintf("ReadSparseMatrix: Could not find file 'features.tsv.gz' at directory '%s'!", path))
  features_ids_types = read.delim(file.path(path, "features.tsv.gz"), header=FALSE, stringsAsFactors=FALSE)
  colnames(features_ids_types) = c("feature_id", "feature_name", "feature_type")
  features_ids_types[, row_name_column] = make.unique(features_ids_types[, row_name_column])
  
  if (!file.exists(file.path(path, "matrix.mtx.gz"))) stop(sprintf("ReadSparseMatrix: Could not find file 'matrix.mtx.gz' at directory '%s'!", path))
  
  # Define a named vector for the row names: its names are the feature names in the dataset and its values are the final names in the Seurat object
  seurat_row_names = setNames(features_ids_types[, row_name_column],features_ids_types[, row_name_column])
  if (!is.null(convert_row_names) & length(convert_row_names) > 0) {
    convert_ids = which(seurat_row_names %in% names(convert_row_names))
    seurat_row_names[convert_ids] = convert_row_names[seurat_row_names[convert_ids]]
  }
  seurat_row_names = setNames(make.unique(gsub(pattern="_", replacement="-", x=seurat_row_names, fixed=TRUE)), names(seurat_row_names))
  rownames(features_ids_types) = unname(seurat_row_names)

  # Read feature data in a list and change row names of the data so that they are Seurat compatible
  feature_data = Seurat::Read10X(path, gene.column=row_name_column)
  if (!is.list(feature_data)) {
    # only one feature type: need to know which and convert into list
    feature_types = unique(features_ids_types[, feature_type_column]) 
    feature_data = list(feature_data)
    names(feature_data) = feature_types[1]
  }
  for (n in names(feature_data)) {
    rownames(feature_data[[n]]) = unname(seurat_row_names[rownames(feature_data[[n]])])
  }

  # special case: if hashtags are specified and if there is a Antibody Capture data, then split and add as separate assay
  # if a regular expression is defined for HTO, overwrite hto_names 
  if ("Antibody Capture" %in% names(feature_data) && !is.null(hto_regex) && nchar(hto_regex)>0) {
    hto_names = grep(pattern=hto_regex, x=rownames(feature_data[["Antibody Capture"]]), v=TRUE, ignore.case=TRUE)
    if (length(hto_names)==0) stop(sprintf("ReadSparseMatrix: Could not find HTO names with the 'hto_regex' argument '%s'!", hto_regex))
  }
  
  if ("Antibody Capture" %in% names(feature_data) & length(hto_names)>0) {
    # check to avoid special chars
    invalid = grep(pattern="[-_]", x=hto_names, v=TRUE)
    if (length(invalid)>0) stop(sprintf("ReadSparseMatrix: The 'hto_names' argument contains invalid (not allowed: -,_) names: %s!", 
                                                     first_n_elements_to_string(invalid)))
    
    # If vector without names, just add values as names
    if (is.null(names(hto_names))) hto_names = setNames(hto_names, hto_names)
    
    # Test if the provided names are in the assay
    missed = names(hto_names[!names(hto_names) %in% rownames(feature_data[["Antibody Capture"]])])
    if (length(missed)>0) stop(sprintf("ReadSparseMatrix: Some of names in the 'hto_names' argument are not present in the 'Antibody Capture' assay: %s!", 
                                                    first_n_elements_to_string(missed)))
    
    # Split assay
    is_hashtag = rownames(feature_data[["Antibody Capture"]]) %in% names(hto_names)
    feature_data[["_HashTags_"]] = feature_data[["Antibody Capture"]][is_hashtag,, drop=FALSE]
    feature_data[["Antibody Capture"]] = feature_data[["Antibody Capture"]][!is_hashtag,, drop=FALSE]
    
    # rename if neccessary
    rownames(feature_data[["_HashTags_"]]) = unname(hto_names[rownames(feature_data[["_HashTags_"]])])
    nms = rownames(features_ids_types)
    nms[nms %in% names(hto_names)] = hto_names[nms[nms %in% names(hto_names)]]
    rownames(features_ids_types) = unname(nms)
    
    # This will be the HashTag assay
    feature_type_to_assay_name = c(feature_type_to_assay_name, c("_HashTags_" = "HTO"))
    
    # Drop Antibody Capture if empty
    if (nrow(feature_data[["Antibody Capture"]])==0) feature_data[["Antibody Capture"]] = NULL
  }
  
  # Read metadata file (if available)
  metadata_table = NULL
  metadata_file = file.path(path, "metadata.tsv.gz")
  if (file.exists(metadata_file)) {
    metadata_table = readr::read_delim(metadata_file, delim="\t", col_names=TRUE, comment="#", progress=FALSE, col_types = readr::cols())
    metadata_table = as.data.frame(metadata_table)
    
    # Check that it contains all cell names
    barcodes = colnames(feature_data[[1]])
    missed = barcodes[!barcodes %in% metadata_table[, 1, drop=TRUE]]
    if (length(missed) > 0) stop(sprintf("ReadSparseMatrix: The 'metadata.tsv.gz' file misses some cell names: %s!", 
                                       first_n_elements_to_string(missed)))
    rownames(metadata_table) = metadata_table[, 1, drop=TRUE]
    
    # Some column names are not allowed since they would be overwritten later
    colnames_metadata_table = colnames(metadata_table)
    invalid = colnames_metadata_table[colnames_metadata_table %in% c("orig.ident", "nCount_RNA", "nFeature_RNA", "nCount_SCT", "nFeature_SCT", "percent_mt", "percent_ercc", "S.Score", "G2M.Score", "Phase", "CC.Difference", "SampleName", "PlateNumber", "PlateRow", "PlateCol", "seurat_clusters")]
    if (length(invalid) > 0) stop(sprintf("ReadSparseMatrix: Some column names in 'metadata.tsv.gz' are not allowed since they would be overwritten: %s!", 
                                        first_n_elements_to_string(invalid)))
  } else {
    metadata_table = data.frame(Cells=colnames(feature_data[[1]]), stringsAsFactors=FALSE)
    rownames(metadata_table) = metadata_table[, 1, drop=TRUE]
  }
  metadata_table$orig.dataset = factor(project)
  
  # Do feature type to assay name
  feature_types = names(feature_data)
  missed = feature_types[!feature_types %in% names(feature_type_to_assay_name)]
  if (length(missed) > 0) stop(sprintf("ReadSparseMatrix: The 'feature_type_to_assay_name' argument misses some feature types: %s!", 
                                                  first_n_elements_to_string(missed)))
  
  # Sort feature types by their order in feature_type_to_assay_name
  feature_types = feature_types[order(match(feature_types, names(feature_type_to_assay_name)))]
  
  # Create Seurat object with first assay
  # Include cell and feature metadata
  sc = list()
  f = feature_types[1] # feature type
  a = feature_type_to_assay_name[f] # cell data
  n = project # project name
  m = metadata_table[, -1, drop=FALSE] # metadata for cells, also drop first column which contains cell names

  sc[[n]] = Seurat::CreateSeuratObject(counts=feature_data[[f]], 
                                  project=project, 
                                  assay=a, 
                                  min.cells=0, 
                                  min.features=0, 
                                  names.delim=NULL, 
                                  names.field=NULL, 
                                  meta.data=m)
  
  # Check that the feature symbols are the same and add feature meta information
  nms = rownames(sc[[n]][[a]])
  missed = nms[!nms %in% rownames(features_ids_types)]
  if (length(missed)>0) stop(sprintf("ReadSparseMatrix: The 'CreateSeuratObject' method modifies feature symbols for assay %s not as expected: %s!", 
                                                  a, first_n_elements_to_string(missed)))
  
  sc[[n]][[a]] = Seurat::AddMetaData(sc[[n]][[a]], features_ids_types[rownames(sc[[n]][[a]]), ])
  
  # Now add remaining assays
  for (f in feature_types[-1]) {

    a = unname(feature_type_to_assay_name[f])
    sc[[n]][[a]] = Seurat::CreateAssayObject(counts=feature_data[[f]], min.cells=0, min.features=0)
    
    nms = rownames(sc[[n]][[a]])
    missed = nms[!nms %in% rownames(features_ids_types)]
    if (length(missed) > 0) stop(sprintf("ReadSparseMatrix: The 'CreateSeuratObject' method modifies feature symbols for assay %s not as expected: %s!", 
                                                    a, first_n_elements_to_string(missed)))
    
    sc[[n]][[a]] = Seurat::AddMetaData(sc[[n]][[a]], features_ids_types[rownames(sc[[n]][[a]]), ])
  }
  
  # Add suffixes if requested
  if (!is.null(cellnames_suffix)) {
    sc[[n]] =  Seurat::RenameCells(sc[[n]], new.names=paste0(gsub("-\\d+$", "", colnames(sc[[n]])), cellnames_suffix))
  }
  
  return(sc)
}

#' Saves the assay data of a Seurat object as a sparse matrix dataset in a format similar to that of 10X.
#' 
#' @param object A seurat object.
#' @param dir Name of the directory.
#' @param assays Which assays to export. Default is all but can be specified.
#' @param slot Which slot to use. Default is counts.
#' @param feature_type_to_assay_name How should the feature types for the different assays be named? Default is: "RNA" = "Gene expression","ADT" = "Antibody Capture" ,"HTO" = "Antibody Capture","Crispr" = "CRISPR Guide Capture", "ERCC" = "ERCC" and "Custom" = "Custom". Also sets the order in which the feature types are written.
#' @param include_cell_metadata_cols Vector with names of cell metadata columns to be written to the metadata.tsv.gz. Can be NULL in which case no file is created.
#' @param include_feature_metadata_cols Vector with names of feature metadata columns to be written. Can be NULL in which case only the first three columns (feature id, feature symbol and feature type) are written.
#' @param metadata_prefix Prefix for cell metadata column names. If NULL, no prefix is prepended. The prefix should end with a separator character that separates the prefix from the original column name.
#' @return The path to the data directory.
ExportSeuratAssayData = function(sc, dir="data", assays=NULL, slot="counts", assay_name_to_feature_type=NULL, include_cell_metadata_cols=NULL, include_feature_metadata_cols=NULL, metadata_prefix=NULL) {
  # Defaults
  if (is.null(assay_name_to_feature_type)) assay_name_to_feature_type = c("RNA"="Gene Expression",
                                                                          "ADT"="Antibody Capture",
                                                                          "HTO"="Multiplexing Capture",
                                                                          "Crispr"="CRISPR Guide Capture",
                                                                          "ERCC"="ERCC",
                                                                          "Custom"="Custom")
  
  # Make the directory and write data
  d = file.path(dir)
  dir.create(d, showWarnings=FALSE)
  
  # Collect assay data and write
  if (is.null(assays)) assays = Seurat::Assays(sc)
  missed = assays[!assays %in% names(assay_name_to_feature_type)]
  if (length(missed) > 0) stop(sprintf("ExportSeuratAssayData: The 'assay_name_to_feature_type' argument misses some assays: %s!", 
                                                  first_n_elements_to_string(missed)))
  
  assays = assays[order(match(assays, names(assay_name_to_feature_type)))]
  feature_data = do.call(rbind, lapply(assays, function(a) { Seurat::GetAssayData(sc, assay=a, slot=slot) }))
  mh = file.path(d, "matrix.mtx")
  Matrix::writeMM(feature_data, file=mh)
  R.utils::gzip(mh, overwrite=TRUE)
  
  # Collect cell barcodes and write
  bh = gzfile(file.path(d, "barcodes.tsv.gz"), open="wb")
  write(colnames(feature_data), file=bh)
  close(bh)
  
  # Collect feature names (plus additional feature meta data) and write
  # If not set, include columns 1:3
  if (is.null(include_feature_metadata_cols)) {
    include_feature_metadata_cols = c(1:3)
  }
  assay_feature_meta_data_df = dplyr::bind_rows(lapply(assays, function(a) { sc[[a]][[include_feature_metadata_cols]] }))
  fh = gzfile(file.path(d, "features.tsv.gz"), open="wb")
  write.table(assay_feature_meta_data_df, file=fh, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
  close(fh)
  
  # Collect cell metadata and write
  if (!is.null(include_cell_metadata_cols) & length(include_cell_metadata_cols) > 0) {
    missed = include_cell_metadata_cols[!include_cell_metadata_cols %in% colnames(sc[[]])]
    if (length(missed) > 0) stop(sprintf("ExportSeuratAssayData: The 'include_cell_metadata_cols' argument contains columns which are not in the metadata ob the Seurat object: %s!", 
                                                    first_n_elements_to_string(missed)))
    
    metadata_table = sc[[include_cell_metadata_cols]]
    if (!is.null(metadata_prefix)) {
      colnames(metadata_table) = paste0(metadata_prefix, colnames(metadata_table))
    }
    metadata_table_nms = colnames(metadata_table)

    metadata_table$CELL_ID_METADATA = rownames(metadata_table)
    metadata_table = metadata_table[, c("CELL_ID_METADATA", metadata_table_nms)]
    
    metah = gzfile(file.path(d, "metadata.tsv.gz"), open="wb")
    write.table(metadata_table, file = metah, row.names=FALSE, col.names=TRUE, quote=TRUE, sep="\t")
    close(metah)
  }  
  
  return(d)
}

#' Parses plate information from the cell names. Mainly used for Smartseq2 datasets where this information is often included in the cell name.
#' 
#' @param cell_names A vector with cell names.
#' @param pattern A regular expression pattern with capture groups for sample name, plate number, row or column which must match the entire cell name. Default is '^(\\S+)_(\\d+)_([A-Z])(\\d+)$'. If the pattern does not match, all information will be set to NA.
#' @param sample_name_group Index of the capture group which contains the sample name. Can be NULL in which case it will not be evaluated and the sample name will be the cell name. Default is 1.
#' @param plate_number_group Index of the capture group which contains the plate number name. Can be NULL in which case it will not be evaluated and the plate number will be NA. Default is 2.
#' @param row_name_group Index of the capture group which contains the plate row name. Can be NULL in which case it will not be evaluated and the plate row name will be NA. Default is 3.
#' @param col_name_group Index of the capture group which contains the plate col name. Can be NULL in which case it will not be evaluated and the plate col name will be NA. Default is 4.
#' @return A data frame with plate information.
ParsePlateInformation = function(cell_names, pattern='^(\\S+)_(\\d+)_([A-Z])(\\d+)$', sample_name_group=1, plate_number_group=2, row_name_group=3, col_name_group=4) {
  # Split cell names
  cell_names_matched_and_split = as.data.frame(stringr::str_match(string=cell_names, pattern=pattern), stringsAsFactors=FALSE)
  cell_names_matched_and_split = merge(x=data.frame(V1=cell_names), y=cell_names_matched_and_split, by=1, all.x=TRUE)
  
  # Now prepare plate information
  plate_information = data.frame(rowname=cell_names_matched_and_split$V1)
  rownames(plate_information) = plate_information$rowname
  
  # Get sample name
  if (!is.null(sample_name_group)) {
    plate_information$SampleName = cell_names_matched_and_split[, 1+sample_name_group]
  } else {
    plate_information$SampleName = cell_names
  }
  plate_information$SampleName = as.factor(plate_information$SampleName)
  
  # Get plate number
  if (!is.null(plate_number_group)) {
    plate_information$PlateNumber = cell_names_matched_and_split[, 1+plate_number_group]
  } else {
    plate_information$PlateNumber = NA
  }
  plate_information$PlateNumber = as.integer(plate_information$PlateNumber)
  
  # Get row
  if (!is.null(row_name_group)) {
    plate_information$PlateRow = cell_names_matched_and_split[, 1+row_name_group]
  } else {
    plate_information$PlateRow = NA
  }
  plate_information$PlateRow = as.character(plate_information$PlateRow)
  
  # Get col
  if (!is.null(col_name_group)) {
    plate_information$PlateCol = cell_names_matched_and_split[, 1+col_name_group]
  } else {
    plate_information$PlateCol = NA
  }
  plate_information$PlateCol = as.integer(plate_information$PlateCol)
  
  # Decide on plate layout
  if ("Q" %in% plate_information$PlateRow | max(c(-Inf,plate_information$PlateCol), na.rm=T) > 24) {
    # super plate?
    plate_information$PlateRow = factor(plate_information$PlateRow, 
                                        levels=c("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z"), 
                                        ordered=TRUE)
    plate_information$PlateCol = factor(plate_information$PlateCol, levels=1:max(c(-Inf,plate_information$PlateCol), na.rm=T), ordered=TRUE)
  } else if ("I" %in% plate_information$PlateRow | max(c(-Inf,plate_information$PlateCol), na.rm=T) > 12) {
    # 384 plate
    plate_information$PlateRow = factor(plate_information$PlateRow, 
                                        levels=c("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P"), 
                                        ordered=TRUE)
    plate_information$PlateCol = factor(plate_information$PlateCol, levels=1:24, ordered=TRUE)
  } else {
    plate_information$PlateRow = factor(plate_information$PlateRow, ordered=TRUE)
    plate_information$PlateCol = factor(plate_information$PlateCol, ordered=TRUE)
  }
  plate_information$rowname = NULL
  
  return(plate_information)
}

#' Exports a Seurat object for visualisation with the cerebroApp (https://github.com/romanhaa/Cerebro).
#' 
#' @param sc Seurat object.
#' @param path File path for export.
#' @param assay Assay to export; default "RNA".
#' @param assay_raw Assay containing the raw data; default "RNA";
#' @param delayed_array Use delayed array for large datasets; default FALSE.
#' @return TRUE if export was successful otherwise FALSE.
ExportToCerebro = function(sc, path, assay="RNA", assay_raw="RNA", delayed_array=FALSE) {
  # Here is an example how the cerebro object is organised 
  #examp_file= system.file('extdata', 'v1.3', 'example.crb', package = 'cerebroApp')
  #examp = readRDS(examp_file)
  
  # Cerebro requires that if certain analyses are done based metadata columns , that these columns are in the "groups" argument
  # For now, this is only neccessary for the trees of the different clustering resolutions
  tree_names = names(Misc(sc, "trees"))
  cerebro_groups = c("orig.ident", "seurat_clusters", setdiff(tree_names, c("orig.ident", "seurat_clusters")))
  
  # Export to cerebro format
  cerebroApp::exportFromSeurat(sc,
                               assay=assay,
                               slot="data",
                               file=path,
                               experiment_name="na",
                               organism="na",
                               groups=cerebro_groups,
                               cell_cycle=c("Phase"),
                               nUMI=paste0("nCount_", assay_raw), 
                               nGene=paste0("nFeature_", assay_raw), 
                               add_all_meta_data=TRUE,
                               use_delayed_array=delayed_array,
                               verbose=FALSE)
  
  # Now read in and modify
  cerebro = readRDS(path)
  
  # Gene lists
  gene_lists = Misc(sc, "gene_lists")
  for(n in names(gene_lists)) cerebro$addGeneList(n, gene_lists[[n]])
  
  # Experiment
  cerebro$addExperiment("experiment_name", Misc(sc, "experiment")$project_id)
  species = Misc(sc, "experiment")$species
  cerebro$addExperiment("organism", ifelse(grepl("sapiens", species), "Hg", ifelse(grepl("musculus", species), "Mm", species)))
  cerebro$addExperiment("date_of_analysis", Misc(sc, "experiment")$date)
  cerebro$addExperiment("date_of_export", Misc(sc, "experiment")$date)
  
  # Marker genes
  if (!is.null(Seurat::Misc(sc, "markers"))) {
    marker_test = Seurat::Misc(sc, "markers")$test
    marker_table = Seurat::Misc(sc, "markers")$results
    
    if (nrow(marker_table) > 0) {
      marker_table = marker_table[, c("cluster", "gene", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")]
      marker_table$on_cell_surface = NA
      colnames(marker_table) = c("seurat_clusters", "gene", "p_val", "avg_logFC", "pct.1", "pct.2", "p_val_adj", "on_cell_surface")
      cerebro$addMarkerGenes(method=marker_test, name="seurat_clusters", table=marker_table)
    }
    
    # Pathways for marker genes
    if ("enrichr" %in% names(Seurat::Misc(sc, "markers")) & nrow(Seurat::Misc(sc, "markers")$enrichr) > 0) {
      marker_enrichr_table = Seurat::Misc(sc, "markers")$enrichr
      marker_enrichr_table$ClusterDirection = paste(marker_enrichr_table$Cluster, marker_enrichr_table$Direction)
      marker_enrichr_table$Overlap = paste(marker_enrichr_table$In.List, marker_enrichr_table$In.Annotation, sep="/")
      marker_enrichr_table$Old.P.value = marker_enrichr_table$P.value
      marker_enrichr_table$Old.Adjusted.P.value= marker_enrichr_table$Adjusted.P.value
      marker_enrichr_table = marker_enrichr_table[, c("ClusterDirection", "Database", "Term", "Overlap", "P.value", "Adjusted.P.value", "Odds.Ratio", "Combined.Score", "Genes")]
      nms = colnames(marker_enrichr_table)
      nms[1] = "seurat_clusters"
      nms[2] = "db"
      colnames(marker_enrichr_table) = nms
      cerebro$addEnrichedPathways(method="enrichr", name="seurat_clusters", table=marker_enrichr_table)
    }
  }
  
  # Cluster tree(s)
  trees = Seurat::Misc(sc, "trees")
  if (!is.null(trees)) {
    if ("orig.ident" %in% names(trees)) cerebro$addTree(group_name="orig.ident", tree=trees$orig.ident)
    if ("seurat_clusters" %in% names(trees)) cerebro$addTree(group_name="seurat_clusters", tree=trees$seurat_clusters)
  }
  
  # Add PCA as projection
  cerebro$addProjection("pca", as.data.frame(Seurat::Embeddings(sc, reduction="pca")[, 1:2]))
  
  # Add DEGs and their pathways as tables
  degs = Misc(sc, "degs")
  if (!is.null(degs)) {
    degs = purrr::flatten(degs)
    for(i in seq(degs)) {
      name = paste0(degs[[i]]$condition_column, ": ", degs[[i]]$condition_group1, " vs ", degs[[i]]$condition_group2)
      if ("subset_column" %in% names(degs[[i]]) && !is.na(degs[[i]]$subset_column)) {
        name = paste0(name, " (", degs[[i]]$subset_column, ": ", degs[[i]]$subset_group, ")")
      }
      
      # DEGs
      if (nrow(degs[[i]]$results) > 0) {
        cerebro$addExtraTable(name, degs[[i]]$results)
      }
      
      # Pathways
      if ("enrichr" %in% names(degs[[i]]) && nrow(degs[[i]]$enrichr) > 0) {
        name = paste(name, "Pathways")
        cerebro$addExtraTable(name, degs[[i]]$enrichr)
      }
    }
  }
  saveRDS(cerebro, path)
  
  return(path)
}

#' Reads a metrics_summary.csv file generated by cellranger count or cellranger multi.
#' 
#' @param file Path to the metrics_summary.csv file.
#' @return A table with a column 'metric' and a column 'value'.
Read10XMetrics = function(file) {
  # Read content and check whether metrics was generated by cellranger count or multi
  content = read.delim(file, sep=",", header=FALSE)
  
  if (nrow(content) == 2) {
    # cellranger count file - use all
    content = content %>% t() %>% as.data.frame()
  } else {
    # cellranger multi file - extract only sample-related stuff
    sample_rows = which(content[, 1] == "Sample")
    content = content[sample_rows, ]
    content = content[, c(5, 6)]
  }
  colnames(content) = c("metric", "value")
  
  idx_num = which(suppressWarnings(!is.na(as.numeric(na.omit(content$value)))))
  content$value[idx_num] = round(as.numeric(content$value[idx_num]), 2)
  
  return(content)
}
