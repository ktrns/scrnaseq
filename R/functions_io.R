#' Reads a counts table (e.g. provided by smartseq2) and converts it into a Seurat object.
#' 
#' @param path Path to a counts table.
#' @param project A project name for the dataset ("SeuratProject").
#' @param feature_name_column Name or index of the column which should be used for feature names (1).
#' @param convert_feature_names Named vector for converting feature names (e.g. from Ensembl ids to gene symbols). Does not have to contain all feature names. Can be NULL in which case feature names are not converted.
#' @param feature_type_column Name or index of the column which should be used for feature types (NULL). Can be NULL in which case feature types are "Gene expression" and "ERCC".
#' @param columns_to_drop Names or indices of the columns which should be dropped (excluding feature_name_column or feature_type_column).
#' @param drop_non_numeric_columns After dropping columns manually with columns_to_drop: Should any remaining non-numeric columns (excluding feature_name_column or feature_type_column) be droppped automatically? If TRUE, they will be dropped silently. If FALSE, then an error message will be printed. 
#' @param feature_type_to_assay_name How should the assays for the different feature type be named? Default is: "Gene expression" = "RNA","Antibody Capture" = "ADT","CRISPR Guide Capture" = "Crispr", "ERCC" = "ERCC" and "Custom" = "Custom". Also sets the order in which the assays are loaded.
#' @param col_sep Separator used in the counts table (tab).
#' @return A Seurat object.
ReadSmartseq2Dataset = function(counts_table, project="SeuratProject", feature_name_column=1, convert_feature_names=NULL, feature_type_column=NULL, columns_to_drop=NULL, drop_non_numeric_columns=TRUE, feature_type_to_assay_name=NULL, sep="\t") {
  library(magrittr)

  # defaults
  if (is.null(feature_type_to_assay_name)) feature_type_to_assay_name = c("Gene Expression"="RNA", 
                                                                          "Antibody Capture"="ADT", 
                                                                          "CRISPR Guide Capture"="Crispr", 
                                                                          "Custom"="Custom", 
                                                                          "ERCC"="ERCC")
  
  # read counts data
  if (!file.exists(counts_table)) futile.logger::flog.error("The counts file %s does not exist!",counts_table)
  feature_data = readr::read_delim(counts_table, delim=sep, col_names=TRUE, comment="#", progress=FALSE)
  
  #  feature_name_column, feature_type_column and feature_columns_to_drop: if (numeric) index, convert to column name (easier)
  if (!is.null(feature_name_column) & is.numeric(feature_name_column)) feature_name_column = colnames(feature_data)[feature_name_column]
  if (!is.null(feature_type_column) & is.numeric(feature_type_column)) feature_type_column = colnames(feature_data)[feature_type_column]
  if (!is.null(columns_to_drop) & length(columns_to_drop)>0 & is.numeric(columns_to_drop)) columns_to_drop = colnames(feature_data)[columns_to_drop]
  
  # drop unwanted columns
  if (!is.null(columns_to_drop) & length(columns_to_drop)>0) {
    feature_data = feature_data[, !colnames(feature_data) %in% columns_to_drop]
  }
  
  # check if there are other non-numeric columns
  invalid = unlist(lapply(colnames(feature_data), function(n) {
    !(n %in% c(feature_name_column, feature_type_column) | is.numeric(feature_data[, n, drop=TRUE]))
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
  features_ids_types = data.frame(feature_id = feature_data[, feature_name_column, drop=TRUE], 
                                  feature_name = feature_data[, feature_name_column, drop=TRUE], 
                                  stringsAsFactors=FALSE)
  feature_data[, feature_name_column] = NULL
  if (is.null(feature_type_column)) {
    features_ids_types$feature_type = ifelse(grepl(pattern="^ERCC-", x=features_ids_types$feature_name), "ERCC", "Gene Expression")
  } else{
    features_ids_types$feature_type = feature_data[, feature_type_column]
    feature_data[, feature_type_column] = NULL
  }
  
  # rename feature names if requested
  if (!is.null(convert_feature_names) & length(convert_feature_names)>0) {
    convert_ids = which(features_ids_types[, "feature_name"] %in% names(convert_feature_names))
    features_ids_types[convert_ids, "feature_name"] = convert_feature_names[features_ids_types[convert_ids, "feature_name"]]
  }
  
  # now make feature names Seurat compatible
  features_ids_types$rowname = gsub(pattern="_", replacement="-", x=features_ids_types$feature_name, fixed=TRUE)
  features_ids_types = features_ids_types %>% 
    dplyr::group_by_at(feature_type_column) %>%
    dplyr::mutate(rowname=make.unique(rowname)) %>% as.data.frame()
  rownames(features_ids_types) = features_ids_types$rowname
  features_ids_types$rowname = NULL
  
  # split by feature type
  feature_data = as.data.frame(feature_data)
  rownames(feature_data) = rownames(features_ids_types)
  feature_data = split(feature_data, as.factor(features_ids_types$feature_type))
  
  # Read metadata file (if available)
  metadata_table = NULL
  
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
  sc = Seurat::CreateSeuratObject(counts=as(as.matrix(feature_data[[f]]), "dgCMatrix"), 
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
  feature_types = feature_types[-1]
  for (f in feature_types) {
    a = feature_type_to_assay_name[f]
    sc[[a]] = Seurat::CreateAssayObject(counts=as(as.matrix(feature_data[[f]]), "dgCMatrix"), min.cells=0, min.features=0)
    
    nms = rownames(sc[[a]])
    missed = nms[!nms %in% rownames(features_ids_types)]
    if (length(missed)>0) futile.logger::flog.error("The 'CreateSeuratObject' method modifies feature symbols for assay %s not as expected: %s!", 
                                                    a, first_n_elements_to_string(missed))
    
    sc[[a]] = Seurat::AddMetaData(sc[[a]],features_ids_types[rownames(sc[[a]]),])
  }
  
  return(sc)
}

#' Reads a sparse matrix dataset (e.g. provided by 10X) and converts it into a Seurat object.
#' 
#' @param path Path to the directory with sparse matrix dataset. Must contain matrix.mtx.gz, features.tsv.gz and barcodes.tsv.gz. Can contain additional metadata in a file metadata.tsv.gz where the first column is the cell name.
#' @param project A project name for the dataset.
#' @param feature_name_column Name or index of the column which should be used for feature names (2).
#' @param convert_feature_names Named vector for converting feature names (e.g. from Ensembl ids to gene symbols). Does not have to contain all feature names. Can be NULL in which case feature names are not converted.
#' @param feature_type_column Name or index of the column which should be used for feature types (3).
#' @param feature_type_to_assay_name How should the assays for the different feature type be named? Default is: "Gene expression" = "RNA","Antibody Capture" = "ADT","CRISPR Guide Capture" = "Crispr", "ERCC" = "ERCC" and "Custom" = "Custom". Also sets the order in which the assays are loaded.
#' @param hto_names If Antibody Capture data, is used hashtags for multiplexing, a vector with feature names to be used as hashtags. If a named vector is provided, the hashtags will be renamed accordingly.
#' @return A Seurat object.
Read10xDataset = function(path, project="SeuratProject", feature_name_column=2, convert_feature_names=NULL, feature_type_column=3, feature_type_to_assay_name=NULL, hto_names=NULL) {
  library(magrittr)
  
  # defaults
  if (is.null(feature_type_to_assay_name)) feature_type_to_assay_name = c("Gene Expression"="RNA",
                                                                          "Antibody Capture"="ADT",
                                                                          "CRISPR Guide Capture"="Crispr",
                                                                          "Custom"="Custom",
                                                                          "ERCC"="ERCC")
  
  # Read barcodes.tsv.gz and features.tsv.gz separately
  if (!file.exists(file.path(path, "barcodes.tsv.gz"))) futile.logger::flog.error("Could not find file 'barcodes.tsv.gz' at directory '%s'!", path)
  barcodes = read.delim(file.path(path, "barcodes.tsv.gz"), header=FALSE, stringsAsFactors=FALSE)[, 1]
  
  if (!file.exists(file.path(path, "features.tsv.gz"))) futile.logger::flog.error("Could not find file 'features.tsv.gz' at directory '%s'!", path)
  features_ids_types = read.delim(file.path(path, "features.tsv.gz"), header=FALSE, stringsAsFactors=FALSE)
  colnames(features_ids_types) = c("feature_id", "feature_name", "feature_type")
  
  if (!file.exists(file.path(path, "matrix.mtx.gz"))) futile.logger::flog.error("Could not find file 'matrix.mtx.gz' at directory '%s'!", path)
  
  # rename feature names if requested
  if (!is.null(convert_feature_names) & length(convert_feature_names)>0) {
    convert_ids = which(features_ids_types[, "feature_name"] %in% names(convert_feature_names))
    features_ids_types[convert_ids, "feature_name"] = convert_feature_names[features_ids_types[convert_ids, "feature_name"]]
  }
  
  # this will make the feature symbols Seurat-compatible
  features_ids_types$rowname = gsub(pattern="_", replacement="-", x=features_ids_types[, feature_name_column], fixed=TRUE)
  features_ids_types = features_ids_types %>% 
    dplyr::group_by_at(feature_type_column) %>%
    dplyr::mutate(rowname=make.unique(rowname)) %>% 
    as.data.frame()
  rownames(features_ids_types) = features_ids_types$rowname
  features_ids_types$rowname = NULL
  
  # Read feature data in a list
  feature_data = Seurat::Read10X(path, gene.column=feature_name_column)
  if (!is.list(feature_data)) {
    # only one feature type: need to know which and convert into list
    feature_types = unique(features_ids_types[, feature_type_column]) 
    feature_data = list(feature_data)
    names(feature_data) = feature_types[1]
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
  feature_types = feature_types[-1]
  for (f in feature_types) {
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
parse_plate_information = function(cell_names, pattern='^(\\S+)_(\\d+)_([A-Z])(\\d+)$', sample_name_group=1, plate_number_group=2, row_name_group=3, col_name_group=4) {
  # split cell names
  cell_names_matched_and_split = as.data.frame(stringr::str_match(string=cell_names, pattern=pattern), stringsAsFactors=FALSE)
  cell_names_matched_and_split = merge(x=data.frame(V1=cell_names), y=cell_names_matched_and_split, by=1, all.x=TRUE)
  
  # now prepare plate information
  plate_information = data.frame(rowname = cell_names_matched_and_split$V1)
  rownames(plate_information) = plate_information$rowname
  
  # sample name
  if (!is.null(sample_name_group)) {
    plate_information$SampleName = cell_names_matched_and_split[, 1+sample_name_group]
  } else {
    plate_information$SampleName = NA
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