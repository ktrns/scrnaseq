#' Reads one or more counts tables (e.g. provided by SmartSeq-2) and converts them into Seurat objects.
# useful environment variables
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
#' @return A list of Seurat objects.
ReadCountsTable = function(counts_table, project="SeuratProject", row_name_column=1, convert_row_names=NULL, feature_type_column=NULL, columns_to_drop=NULL, drop_non_numeric_columns=TRUE, feature_type_to_assay_name=NULL, sep="\t", parse_plate_information=TRUE, plate_information_regex='^(\\S+)_(\\d+)_([A-Z])(\\d+)$', sample_name_group=1, plate_number_group=2, row_name_group=3, col_name_group=4, return_samples_as_datasets=TRUE) {
  library(magrittr)
  
  # defaults
  if (is.null(feature_type_to_assay_name)) feature_type_to_assay_name = c("Gene Expression"="RNA", 
                                                                          "Antibody Capture"="ADT", 
                                                                          "CRISPR Guide Capture"="Crispr", 
                                                                          "Custom"="Custom", 
                                                                          "ERCC"="ERCC")
  
  # read counts data
  if (!file.exists(counts_table)) futile.logger::flog.error("The counts file %s does not exist!",counts_table)
  feature_data = readr::read_delim(counts_table, delim=sep, col_names=TRUE, comment="#", progress=FALSE, col_types = readr::cols())
  
  #  row_name_column, feature_type_column and columns_to_drop: if (numeric) index, convert to column name (easier)
  if (!is.null(row_name_column) & is.numeric(row_name_column)) row_name_column = colnames(feature_data)[row_name_column]
  if (!is.null(feature_type_column) & is.numeric(feature_type_column)) feature_type_column = colnames(feature_data)[feature_type_column]
  if (!is.null(columns_to_drop) & length(columns_to_drop)>0 & is.numeric(columns_to_drop)) columns_to_drop = colnames(feature_data)[columns_to_drop]
  
  # drop unwanted columns
  if (!is.null(columns_to_drop) & length(columns_to_drop)>0) {
    feature_data = feature_data[, !colnames(feature_data) %in% columns_to_drop]
  }
  
  # check if there are other non-numeric columns
  invalid = unlist(lapply(colnames(feature_data), function(n) {
    !(n %in% c(row_name_column, feature_type_column) | is.numeric(feature_data[, n, drop=TRUE]))
  }))
  if (sum(invalid)>0) {
    if (drop_non_numeric_columns) {
      # if yes and drop_non_numeric_columns: just remove them
      feature_data = feature_data[, !invalid]
    } else {
      # otherwise print an error
      invalid = invalid[which(invalid)]
      futile.logger::flog.error("Some columns in the counts table do not contain numeric data: %s!", 
                                first_n_elements_to_string(names(invalid)))
    }
  }
  
  # create a data frame for feature id, feature name and feature type (similar to features.tsv in 10x datasets)
  # feature_type_column: indicates the type of feature; if null, assume, there is only "Gene Expression" features; 
  #   only check for ERCC controls which will be separate
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
  
  # Define a named vector for the row names: its names are the feature names in the dataset and its values are the final names in the Seurat object
  seurat_row_names = setNames(features_ids_types[, "feature_id"],features_ids_types[, "feature_id"])
  if (!is.null(convert_row_names) & length(convert_row_names)>0) {
    convert_ids = which(seurat_row_names %in% names(convert_row_names))
    seurat_row_names[convert_ids] = convert_row_names[seurat_row_names[convert_ids]]
  }
  seurat_row_names = setNames(make.unique(gsub(pattern="_", replacement="-", x=seurat_row_names, fixed=TRUE)),names(seurat_row_names))
  rownames(features_ids_types) = unname(seurat_row_names)
  
  # split by feature type and add Seurat-compatible row names
  feature_data = as.data.frame(feature_data)
  rownames(feature_data) = rownames(features_ids_types)
  feature_data = split(feature_data, as.factor(features_ids_types$feature_type))
  
  # feature type to assay name
  feature_types = names(feature_data)
  missed = feature_types[!feature_types %in% names(feature_type_to_assay_name)]
  if (length(missed)>0) futile.logger::flog.error("The 'feature_type_to_assay_name' argument misses some feature types: %s!", 
                                                  first_n_elements_to_string(missed))
  
  # sort feature types by their order in feature_type_to_assay_name
  feature_types = feature_types[order(match(feature_types, names(feature_type_to_assay_name)))]
  
  # Read metadata file (if available)
  metadata_file = file.path(dirname(path), "metadata.tsv.gz")
  if (file.exists(metadata_file)) {
    metadata_table = read.delim(metadata_file, header=TRUE, stringsAsFactors=FALSE)
    barcodes = colnames(feature_data[[1]])
    missed = barcodes[!barcodes %in% metadata_table[, 1]]
    if (length(missed)>0) futile.logger::flog.error("The 'metadata.tsv.gz' file misses some cell names: %s!", 
                                                    first_n_elements_to_string(missed))
    rownames(metadata_table) = metadata_table[, 1]
  } else {
    metadata_table = data.frame(Cells=colnames(feature_data[[1]]) ,stringsAsFactors=FALSE)
    rownames(metadata_table) = metadata_table[, 1]
  }
  
  # If a regex for parsing the plate information from the names has been provided, parse: sample name, plate number, plate row and plate column
  if (parse_plate_information) {
    plate_information = ParsePlateInformation(colnames(feature_data[[1]]), pattern=plate_information_regex, sample_name_group=sample_name_group, plate_number_group=plate_number_group, row_name_group=row_name_group, col_name_group=col_name_group)
    
    # Then add to metadata
    if (nrow(plate_information) > 0) {
      
      metadata_table = merge(metadata_table,plate_information,by="row.names", all.x=TRUE)
      rownames(metadata_table) = metadata_table$Row.names
      metadata_table$Row.names = NULL
    } else {
      futile.logger::flog.warn("The 'ReadCountsTable' method could not parse plate information from the names though requested! Check the regular expression ('plate_information_regex'). Cannot return dataset by samples ('return_samples_as_datasets')!")
      return_samples_as_datasets = FALSE
    }
  }

  # Group cells by samples for datasets
  samples_to_process = list()
  if (return_samples_as_datasets) {
    if(!parse_plate_information) futile.logger::flog.error("The 'ReadCountsTable' method cannot return datasets by samples without parsing this information from the name ('parse_plate_information')!")
    samples_to_process = split(rownames(metadata_table), metadata_table$SampleName)
  } else {
    samples_to_process[[project]] = rownames(metadata_table)
  }
  
  # Now create Seurat objects
  sc = list()
  for (s in names(samples_to_process)) {
    
    # New name
    n = paste(project, s, sep=".")
    
    # Create Seurat object with first assay
    # Include cell and feature metadata
    c = samples_to_process[[s]] # cells to include for sample
    f = feature_types[1] # feature type
    a = feature_type_to_assay_name[f] # assay name
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
    
    # check that the feature symbols are the same and add feature meta information
    nms = rownames(sc[[n]][[a]])
    missed = nms[!nms %in% rownames(features_ids_types)]
    if (length(missed)>0) futile.logger::flog.error("The 'CreateSeuratObject' method modifies feature symbols for assay %s not as expected: %s!", 
                                                  a, first_n_elements_to_string(missed))
    sc[[n]][[a]] = Seurat::AddMetaData(sc[[n]][[a]], features_ids_types[rownames(sc[[n]][[a]]),])
    
    # now add remaining assays
    for (f in feature_types[-1]) {
      a = feature_type_to_assay_name[f]
      d = as(as.matrix(feature_data[[f]][,c, drop=FALSE]), "dgCMatrix")
      sc[[n]][[a]] = Seurat::CreateAssayObject(counts=d,
                                          min.cells=0,
                                          min.features=0)
      nms = rownames(sc[[n]][[a]])
      missed = nms[!nms %in% rownames(features_ids_types)]
      if (length(missed)>0) futile.logger::flog.error("The 'CreateAssayObject' method modifies feature symbols for assay %s not as expected: %s!", 
                                                    a, first_n_elements_to_string(missed))
      sc[[n]][[a]] = Seurat::AddMetaData(sc[[n]][[a]],features_ids_types[rownames(sc[[n]][[a]]),])
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
#' @param hto_names If Antibody Capture data, is used hashtags for multiplexing, a vector with feature names to be used as hashtags. If a named vector is provided, the hashtags will be renamed accordingly.
#' @param return_samples_as_datasets If there are multiple samples in the dataset: return each as separate dataset ("TRUE").
#' @return A Seurat object.
ReadSparseMatrix = function(path, project="SeuratProject", row_name_column=2, convert_row_names=NULL, feature_type_column=3, feature_type_to_assay_name=NULL, hto_names=NULL) {
  library(magrittr)
  
  # defaults
  if (is.null(feature_type_to_assay_name)) feature_type_to_assay_name = c("Gene Expression"="RNA",
                                                                          "Antibody Capture"="ADT",
                                                                          "CRISPR Guide Capture"="Crispr",
                                                                          "Custom"="Custom",
                                                                          "ERCC"="ERCC")
  
  # Check that files exist and read barcodes.tsv.gz and features.tsv.gz separately
  if (!file.exists(file.path(path, "barcodes.tsv.gz"))) futile.logger::flog.error("Could not find file 'barcodes.tsv.gz' at directory '%s'!", path)
  barcodes = read.delim(file.path(path, "barcodes.tsv.gz"), header=FALSE, stringsAsFactors=FALSE)[, 1]
  
  if (!file.exists(file.path(path, "features.tsv.gz"))) futile.logger::flog.error("Could not find file 'features.tsv.gz' at directory '%s'!", path)
  features_ids_types = read.delim(file.path(path, "features.tsv.gz"), header=FALSE, stringsAsFactors=FALSE)
  colnames(features_ids_types) = c("feature_id", "feature_name", "feature_type")
  
  if (!file.exists(file.path(path, "matrix.mtx.gz"))) futile.logger::flog.error("Could not find file 'matrix.mtx.gz' at directory '%s'!", path)
  
  # Define a named vector for the row names: its names are the feature names in the dataset and its values are the final names in the Seurat object
  seurat_row_names = setNames(features_ids_types[, row_name_column],features_ids_types[, row_name_column])
  if (!is.null(convert_row_names) & length(convert_row_names)>0) {
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
  if ("Antibody Capture" %in% names(feature_data) & length(hto_names)>0) {
    # check to avoid special chars
    invalid = grep(pattern="[-_]", x=hto_names, v=TRUE)
    if (length(invalid)>0) futile.logger::flog.error("The 'hto_names' argument contains invalid (not allowed: -,_) names: %s!", 
                                                     first_n_elements_to_string(invalid))
    
    # if vector without names, just add values as names
    if (is.null(names(hto_names))) hto_names = setNames(hto_names, hto_names)
    
    # test if the provided names are in the assay
    missed = names(hto_names[!names(hto_names) %in% rownames(feature_data[["Antibody Capture"]])])
    if (length(missed)>0) futile.logger::flog.error("Some of names in the 'hto_names' argument are not present in the 'Antibody Capture' assay: %s!", 
                                                    first_n_elements_to_string(missed))
    
    # split assay
    is_hashtag = rownames(feature_data[["Antibody Capture"]]) %in% names(hto_names)
    feature_data[["_HashTags_"]] = feature_data[["Antibody Capture"]][is_hashtag,, drop=FALSE]
    feature_data[["Antibody Capture"]] = feature_data[["Antibody Capture"]][!is_hashtag,, drop=FALSE]
    
    # rename if neccessary
    rownames(feature_data[["_HashTags_"]]) = hto_names[rownames(feature_data[["_HashTags_"]])]
    nms = rownames(features_ids_types)
    nms[nms %in% names(hto_names)] = hto_names[nms[nms %in% names(hto_names)]]
    rownames(features_ids_types) = nms
    
    # this will be the HashTag assay
    feature_type_to_assay_name = c(feature_type_to_assay_name, c("_HashTags_" = "HTO"))
    
    # drop Antibody Capture if empty
    if (nrow(feature_data[["Antibody Capture"]])==0) feature_data[["Antibody Capture"]] = NULL
  }
  
  # Read metadata file (if available)
  metadata_table = NULL
  metadata_file = file.path(path, "metadata.tsv.gz")
  if (file.exists(metadata_file)) {
    metadata_table = read.delim(metadata_file, header=TRUE, stringsAsFactors=FALSE)
    missed = barcodes[!barcodes %in% metadata_table[, 1]]
    if (length(missed)>0) futile.logger::flog.error("The 'metadata.tsv.gz' file misses some cell names: %s!", 
                                                    first_n_elements_to_string(missed))
    rownames(metadata_table) = metadata_table[, 1]
  }
  
  # feature type to assay name
  feature_types = names(feature_data)
  missed = feature_types[!feature_types %in% names(feature_type_to_assay_name)]
  if (length(missed)>0) futile.logger::flog.error("The 'feature_type_to_assay_name' argument misses some feature types: %s!", 
                                                  first_n_elements_to_string(missed))
  
  # sort feature types by their order in feature_type_to_assay_name
  feature_types = feature_types[order(match(feature_types, names(feature_type_to_assay_name)))]
  
  # create Seurat object with first assay
  # include cell and feature metadata
  f = feature_types[1]
  a = feature_type_to_assay_name[f]
  sc = Seurat::CreateSeuratObject(counts=feature_data[[f]], 
                                  project=project, 
                                  assay=a, 
                                  min.cells=0, 
                                  min.features=0, 
                                  names.delim=NULL, 
                                  names.field=NULL, 
                                  meta.data=metadata_table)
  
  # check that the feature symbols are the same and add feature meta information
  nms = rownames(sc[[a]])
  missed = nms[!nms %in% rownames(features_ids_types)]
  if (length(missed)>0) futile.logger::flog.error("The 'CreateSeuratObject' method modifies feature symbols for assay %s not as expected: %s!", 
                                                  a, first_n_elements_to_string(missed))
  
  sc[[a]] = Seurat::AddMetaData(sc[[a]], features_ids_types[rownames(sc[[a]]),])
  
  # now add remaining assays
  for (f in feature_types[-1]) {
    a = feature_type_to_assay_name[f]
    sc[[a]] = Seurat::CreateAssayObject(counts=feature_data[[f]], min.cells=0, min.features=0)
    
    nms = rownames(sc[[a]])
    missed = nms[!nms %in% rownames(features_ids_types)]
    if (length(missed)>0) futile.logger::flog.error("The 'CreateSeuratObject' method modifies feature symbols for assay %s not as expected: %s!", 
                                                    a, first_n_elements_to_string(missed))
    
    sc[[a]] = Seurat::AddMetaData(sc[[a]], features_ids_types[rownames(sc[[a]]),])
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
#' @return The path to the data directory.
ExportSeuratAssayData = function(sc, dir="data", assays=NULL, slot="counts", assay_name_to_feature_type=NULL, include_cell_metadata_cols=NULL, include_feature_metadata_cols=NULL) {
  # defaults
  if (is.null(assay_name_to_feature_type)) assay_name_to_feature_type = c("RNA"="Gene Expression",
                                                                          "ADT"="Antibody Capture",
                                                                          "HTO"="Antibody Capture",
                                                                          "Crispr"="CRISPR Guide Capture",
                                                                          "ERCC"="ERCC",
                                                                          "Custom"="Custom")
  
  # make the directory and write data
  d = file.path(dir)
  dir.create(d, showWarnings=FALSE)
  
  # collect assay data and write
  if (is.null(assays)) assays = Seurat::Assays(sc)
  missed = assays[!assays %in% names(assay_name_to_feature_type)]
  if (length(missed)>0) futile.logger::flog.error("The 'assay_name_to_feature_type' argument misses some assays: %s!", 
                                                  first_n_elements_to_string(missed))
  
  assays = assays[order(match(assays, names(assay_name_to_feature_type)))]
  feature_data = do.call(rbind, lapply(assays, function(a) { Seurat::GetAssayData(sc, assay=a, slot=slot) }))
  mh = file.path(d, "matrix.mtx")
  Matrix::writeMM(feature_data, file=mh)
  R.utils::gzip(mh, overwrite=TRUE)
  
  # collect cell barcodes and write
  bh = gzfile(file.path(d, "barcodes.tsv.gz"), open="wb")
  write(colnames(feature_data), file=bh)
  close(bh)
  
  # collect feature names (plus additional feature meta data) and write
  # if not set, include columns 1:3
  if (is.null(include_feature_metadata_cols)) {
    include_feature_metadata_cols = c(1:3)
  }
  assay_feature_meta_data_df = dplyr::bind_rows(lapply(assays, function(a) { sc[[a]][[include_feature_metadata_cols]] }))
  fh = gzfile(file.path(d, "features.tsv.gz"), open="wb")
  write.table(assay_feature_meta_data_df, file=fh, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
  close(fh)
  
  # collect cell metadata and write
  if (!is.null(include_cell_metadata_cols) & length(include_cell_metadata_cols)>0) {
    missed = include_cell_metadata_cols[!include_cell_metadata_cols %in% colnames(sc[[]])]
    if (length(missed)>0) futile.logger::flog.error("The 'include_cell_metadata_cols' argument contains columns which are not in the metadata ob the Seurat object: %s!", 
                                                    first_n_elements_to_string(missed))
    
    metadata_table = sc[[include_cell_metadata_cols]]
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
#' @param pattern A regular expression pattern with capture groups for sample name, plate number, row or column which must match the entire cell name. Default is '^(\\S+)_(\\d+)_([A-Z])(\\d+)$'. Information for which no capture groups will be set to NA.
#' @param sample_name_group Index of the capture group which contains the sample name. Can be NULL in which case it will not be evaluated and the sample name will be NA. Default is 1.
#' @param plate_number_group Index of the capture group which contains the plate number name. Can be NULL in which case it will not be evaluated and the plate number will be NA. Default is 2.
#' @param row_name_group Index of the capture group which contains the plate row name. Can be NULL in which case it will not be evaluated and the plate row name will be NA. Default is 3.
#' @param col_name_group Index of the capture group which contains the plate col name. Can be NULL in which case it will not be evaluated and the plate col name will be NA. Default is 4.
#' @return A data frame with plate information.
ParsePlateInformation = function(cell_names, pattern='^(\\S+)_(\\d+)_([A-Z])(\\d+)$', sample_name_group=1, plate_number_group=2, row_name_group=3, col_name_group=4) {
  # split cell names
  cell_names_matched_and_split = as.data.frame(stringr::str_match(string=cell_names, pattern=pattern), stringsAsFactors=FALSE)
  cell_names_matched_and_split = merge(x=data.frame(V1=cell_names), y=cell_names_matched_and_split, by=1, all.x=TRUE)
  
  # now prepare plate information
  plate_information = data.frame(rowname=cell_names_matched_and_split$V1)
  rownames(plate_information) = plate_information$rowname
  
  # sample name
  if (!is.null(sample_name_group)) {
    plate_information$SampleName = cell_names_matched_and_split[, 1+sample_name_group]
  } else {
  }
  plate_information$SampleName = as.factor(plate_information$SampleName)
  
  # plate number
  if (!is.null(plate_number_group)) {
    plate_information$PlateNumber = cell_names_matched_and_split[, 1+plate_number_group]
  } else {
    plate_information$PlateNumber = NA
  }
  plate_information$PlateNumber = as.integer(plate_information$PlateNumber)
  
  # row
  if (!is.null(row_name_group)) {
    plate_information$PlateRow = cell_names_matched_and_split[, 1+row_name_group]
  } else {
    plate_information$PlateRow = NA
  }
  plate_information$PlateRow = as.character(plate_information$PlateRow)
  
  # col
  if (!is.null(col_name_group)) {
    plate_information$PlateCol = cell_names_matched_and_split[, 1+col_name_group]
  } else {
    plate_information$PlateCol = NA
  }
  plate_information$PlateCol = as.integer(plate_information$PlateCol)
  
  # decide on plate layout
  if ("Q" %in% plate_information$PlateRow | max(plate_information$PlateCol)>24) {
    # super plate?
    plate_information$PlateRow = factor(plate_information$PlateRow, 
                                        levels=c("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z"), 
                                        ordered=TRUE)
    plate_information$PlateCol = factor(plate_information$PlateCol, levels=1:max(plate_information$PlateCol), ordered=TRUE)
  } else if ("I" %in% plate_information$PlateRow | max(plate_information$PlateCol)>12) {
    # 384 plate
    plate_information$PlateRow = factor(plate_information$PlateRow, 
                                        levels=c("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P"), 
                                        ordered=TRUE)
    plate_information$PlateCol = factor(plate_information$PlateCol, levels=1:24, ordered=TRUE)
  } else {
    plate_information$PlateRow = factor(plate_information$PlateRow, 
                                        levels=c("A","B","C","D","E","F","G","H"),
                                        ordered=TRUE)
    plate_information$PlateCol = factor(plate_information$PlateCol, levels=1:12, ordered=TRUE)
  }
  plate_information$rowname = NULL
  
  return(plate_information)
}

#' Exports a Seurat object for visualisation with the cerebroApp (https://github.com/romanhaa/Cerebro).
#' 
#' @param sc A Seurat object.
#' @param path File path for export.
#' @param param Named list with parameters passed to the scrnaseq script.
#' @param organism Species. Can be 'Hg', 'Mm', ... .
#' @param assay Assay to export ('RNA').
#' @param column_sample Metadata column containing the sample name ('orig.ident').
#' @param column_cluster Metadata column containing the cluster identity ('seurat_clusters').
#' @param column_ccphase Metadata column containing the cell cycle phase ('Phase').
#' @param gene_lists Gene lists to include (NULL).
#' @param marker_genes Marker genes table obtained from Seurat::FindMarkers (NULL).
#' @param enriched_pathways List with enriched pathways results.
#' @return TRUE if export was successful otherwise FALSE.
ExportToCerebro = function(sc, path, param, project="scrnaseq", species, assay="RNA", column_sample="orig.ident", column_cluster="seurat_clusters", gene_lists = NULL, marker_genes=NULL, enriched_pathways=NULL, column_ccphase="Phase") {
  
  # Set to NULL if the requested CC phase is not part of the actual metadata
  if (!(column_ccphase %in% colnames(sc[[]]))) column_ccphase = NULL
  
  
  # Save the Seurat object in the cerebroApp format
  cerebroApp::exportFromSeurat(sc,
                               file=path,
                               assay=assay,
                               experiment_name=project,
                               organism=species,
                               column_sample=column_sample,
                               column_cluster=column_cluster,
                               column_nUMI=paste("nCount", assay, sep="_"),
                               column_nGene=paste("nFeature", assay, sep="_"),
                               column_cell_cycle_seurat=column_ccphase,
                               add_all_meta_data=T)
  
  # Export function does not include all information - need to add manually
  crb_obj = readRDS(path)
  
  # Add basic parameters (are hardcoded unfortunately)
  crb_obj$parameters[["gene_nomenclature"]] = "gene_name"
  crb_obj$parameters[["discard_genes_expressed_in_fewer_cells_than"]] = "NA"
  crb_obj$parameters[["keep_mitochondrial_genes"]] = TRUE
  crb_obj$parameters[["variables_to_regress_out"]] = paste(param$vars_to_regress, collapse=",")
  crb_obj$parameters[["number_PCs"]] = param$pc_n
  crb_obj$parameters[["cluster_resolution"]] = param$cluster_resolution
  crb_obj$parameters[["filtering"]] = list(UMI_min="NA", UMI_max="NA", genes_min="NA", genes_max="NA")
  crb_obj$parameters[["tSNE_perplexity"]] = "NA"
  
  # Add gene lists (G2M_phase_genes, S_phase_genes, mitochondrial_genes may be hardcoded; additional lists may be custom)
  if(!is.null(gene_lists) & length(gene_lists)>0) crb_obj$gene_lists = gene_lists
  
  # Add technical information
  crb_obj$technical_info = list(R=R.utils::captureOutput(devtools::session_info()), seurat_version=packageVersion("Seurat"))
  
  # Add sample-related information
  crb_obj$samples = list()
  crb_obj$samples[["colors"]] = param$col_samples
  crb_obj$samples[["overview"]] = data.frame(sample=levels(sc[[column_sample, drop=TRUE]]))
  
  crb_obj$samples[["by_cluster"]] = data.frame(sample=sc[[column_sample, drop=TRUE]], cluster=sc[[column_cluster, drop=TRUE]]) %>%
    dplyr::group_by(sample, cluster) %>%
    dplyr::summarise(num_cells=length(cluster)) %>%
    dplyr::group_by(sample) %>%
    dplyr::mutate(total_cell_count=sum(num_cells)) %>%
    tidyr::pivot_wider(names_from=cluster, values_from=num_cells, values_fill=list(num_cells=0))
  
  if ("Phase" %in% colnames(sc[[]])) {
    crb_obj$samples[["by_cell_cycle_seurat"]] = data.frame(sample=sc[[column_sample, drop=TRUE]], phase=sc[[column_ccphase, drop=TRUE]]) %>% 
      dplyr::group_by(sample, phase) %>%
      dplyr::summarise(num_cells=length(phase)) %>%
      dplyr::group_by(sample) %>%
      dplyr::mutate(total_cell_count=sum(num_cells)) %>%
      tidyr::pivot_wider(names_from=phase, values_from=num_cells, values_fill=list(num_cells=0))
  }
  
  # Add cluster-related information
  crb_obj$clusters = list()
  crb_obj$clusters[["colors"]] = param$col_clusters
  crb_obj$clusters[["overview"]] = data.frame(cluster=levels(sc[[column_cluster, drop=TRUE]]))
  
  crb_obj$clusters[["by_samples"]] = data.frame(sample=sc[[column_sample, drop=TRUE]], cluster=sc[[column_cluster, drop=TRUE]]) %>%
    dplyr::group_by(sample, cluster) %>%
    dplyr::summarise(num_cells=length(cluster)) %>%
    dplyr::group_by(cluster) %>%
    dplyr::mutate(total_cell_count=sum(num_cells)) %>%
    tidyr::pivot_wider(names_from=sample, values_from=num_cells, values_fill=list(num_cells=0))
  
  if ("Phase" %in% colnames(sc[[]])) {
    crb_obj$clusters[["by_cell_cycle_seurat"]] = data.frame(cluster=sc[[column_cluster, drop=TRUE]], phase=sc[[column_ccphase, drop=TRUE]]) %>% 
      dplyr::group_by(cluster, phase) %>%
      dplyr::summarise(num_cells=length(phase)) %>%
      dplyr::group_by(cluster) %>%
      dplyr::mutate(total_cell_count=sum(num_cells)) %>%
      tidyr::pivot_wider(names_from=phase, values_from=num_cells, values_fill=list(num_cells=0))
  }
  
  tree = Tool(object=sc, slot="BuildClusterTree")
  if(!is.null(tree)) crb_obj$clusters[["tree"]] = tree
  
  # Fix cell cycle information
  if ("Phase" %in% colnames(sc[[]])) {
    crb_obj$cells$cell_cycle_seurat = factor(crb_obj$cells$Phase, levels=c("G1","G2M","S"))
  }
  
  # Add marker genes
  crb_obj$marker_genes = list()
  crb_obj$marker_genes[["parameters"]] = list(only_positive=FALSE, minimum_percentage=0.25, logFC_threshold=param$log2fc, test="MAST", p_value_threshold=param$padj)
  if(!is.null(marker_genes) & nrow(marker_genes)>0) {
    crb_obj$marker_genes[["by_cluster"]] = marker_genes %>% dplyr::select(c("cluster", "gene", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj"))
  } else {
    crb_obj$marker_genes[["by_cluster"]] = "no markers found"
  }
  crb_obj$marker_genes[["by_sample"]]
  
  # Add enriched pathways
  crb_obj$enriched_pathways = list(enrichr=list())
  crb_obj$enriched_pathways$enrichr$parameters = list(databases=param$enrichr_dbs, adj_p_cutoff="NA", max_terms="NA")
  
  pathways_by_cluster = purrr::map_dfr(names(enriched_pathways), function(x){
    d = purrr::flatten_dfr(enriched_pathways[x], .id="db")
    if (nrow(d) == 0) return(d)
    d$set = x
    d$cluster = gsub("\\S+_cluster_","",d$set)
    return(d[,c(ncol(d), ncol(d)-1, 1, 2:(ncol(d)-2))])
  })
  pathways_by_cluster$cluster = factor(pathways_by_cluster$cluster, levels=unique(pathways_by_cluster$cluster))
  pathways_by_cluster$db = factor(pathways_by_cluster$db, levels=unique(pathways_by_cluster$db))
  pathways_by_cluster$Term = paste(pathways_by_cluster$Term, ifelse(grepl("DEG_up", pathways_by_cluster$set), "UP", "DOWN"))
  pathways_by_cluster$set = NULL
  
  crb_obj$enriched_pathways$enrichr[["by_cluster"]] = pathways_by_cluster
  crb_obj$enriched_pathways$enrichr[["by_sample"]] = "no enrichment found"
  
  # Save modified object
  saveRDS(crb_obj, file=path)
  
  return(TRUE)
}
