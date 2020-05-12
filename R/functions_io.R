#' Reads a sparse matrix dataset (e.g. provided by 10X) and converts it into a Seurat object.
#' 
#' @param path Path to the directory with sparse matrix dataset. Must contain matrix.mtx.gz, features.tsv.gz and barcodes.tsv.gz. Can contain additional metadata in a file metadata.tsv.gz where the first column is the cell name.
#' @param project A project name for the dataset.
#' @param gene_column Name or index of the column which should be used for feature names.
#' @param unique_features Make feature names unique.
#' @param feature_type_to_assay_name How should the assays for the different feature type be named? Default is: "Gene expression" = "RNA","Antibody Capture" = "ADT","CRISPR Guide Capture" = "Crispr" and "Custom" = "Custom". Also sets the order in which the assays are loaded.
#' @param hto_names If Antibody Capture data, is used hashtags for multiplexing, a vector with feature names to be used as hashtags. If a named vector is provided, the hashtags will be renamed accordingly.
#' @return A Seurat object.
read_10x_dataset = function(path, project = "SeuratProject", feature_column = 2, feature_type_column = 3, unique_features = TRUE, feature_type_to_assay_name = NULL, hto_names = NULL){
  # dplyr, Seurat, futile.logger; there may be a better solution
  library(dplyr)
  library(Seurat)
  library(futile.logger)
  library(magrittr)
  
  # defaults
  if(is.null(feature_type_to_assay_name)) feature_type_to_assay_name = c("Gene Expression" = "RNA","Antibody Capture" = "ADT","CRISPR Guide Capture" = "Crispr","Custom" = "Custom")
  
  # Read barcodes.tsv.gz and features.tsv.gz separately
  barcodes = read.delim(file.path(path, "barcodes.tsv.gz"), header = FALSE,stringsAsFactors = FALSE)[,1]
  features_ids_types = read.delim(file.path(path, "features.tsv.gz"), header = FALSE,stringsAsFactors = FALSE)
  colnames(features_ids_types) = c("orig_feature_id","orig_feature_symbol","orig_feature_type")
  
  # this will make the feature symbols Seurat-compatible
  features_ids_types$rowname = gsub(pattern="_",replacement = "-",features_ids_types[,feature_column],fixed=TRUE)
  if(unique_features){
    features_ids_types = features_ids_types %>% 
      group_by_at(feature_type_column) %>%
      mutate(rowname=make.unique(rowname)) %>% as.data.frame()
  }
  rownames(features_ids_types) = features_ids_types$rowname
  features_ids_types$rowname = NULL
  
  # Read feature data in a list
  feature_data = Read10X(path,gene.column = feature_column,unique.features = unique_features)
  if(!is.list(feature_data)){
    # only one feature type: need to know which and convert into list
    feature_types = unique(features_ids_types[,feature_type_column]) 
    feature_data = list(feature_data)
    names(feature_data) = feature_types[1]
  }
  
  # special case: if hashtags are specified and if there is a Antibody Capture data, then split and add as separate assay
  if("Antibody Capture" %in% names(feature_data) & length(hto_names)>0){
    # check to avoid special chars
    invalid = grep("[-_]",hto_names,v=T)
    if(length(invalid)>0) flog.error("The 'hto_names' argument contains invalid (not allowed: -,_) names: %s!",first_n_elements_to_string(invalid))
    
    # if vector without names, just add values as names
    if(is.null(names(hto_names))) hto_names = setNames(hto_names,hto_names)
    
    # test if the provided names are in the assay
    missed = names(hto_names[!names(hto_names) %in% rownames(feature_data[["Antibody Capture"]])])
    if(length(missed)>0) flog.error("Some of names in the 'hto_names' argument are not present in the 'Antibody Capture' assay: %s!",first_n_elements_to_string(missed))
    
    # split assay
    is_hashtag = rownames(feature_data[["Antibody Capture"]]) %in% names(hto_names)
    feature_data[["_HashTags_"]] = feature_data[["Antibody Capture"]][is_hashtag,,drop=F]
    feature_data[["Antibody Capture"]] = feature_data[["Antibody Capture"]][!is_hashtag,,drop=F]
    
    # rename if neccessary
    rownames(feature_data[["_HashTags_"]]) = hto_names[rownames(feature_data[["_HashTags_"]])]
    nms = rownames(features_ids_types)
    nms[nms %in% names(hto_names)] = hto_names[nms[nms %in% names(hto_names)]]
    rownames(features_ids_types) = nms
    
    # this will be the HashTag assay
    feature_type_to_assay_name = c(feature_type_to_assay_name,c("_HashTags_" = "HTO"))
    
    # drop Antibody Capture if empty
    if(nrow(feature_data[["Antibody Capture"]])==0) feature_data[["Antibody Capture"]] = NULL
  }
  
  # Read metadata file (if available)
  metadata_table = NULL
  metadata_file = file.path(path,"metadata.tsv.gz")
  if(file.exists(metadata_file)){
    metadata_table = read.delim(metadata_file,header = TRUE,stringsAsFactors = FALSE)
    missed = barcodes[!barcodes %in% metadata_table[,1]]
    if(length(missed)>0) flog.error("The 'metadata.tsv.gz' file misses some cell names: %s!",first_n_elements_to_string(missed))
    rownames(metadata_table) = metadata_table[,1]
  }
  
  # feature type to assay name
  feature_types = names(feature_data)
  missed = feature_types[!feature_types %in% names(feature_type_to_assay_name)]
  if(length(missed)>0) flog.error("The 'feature_type_to_assay_name' argument misses some feature types: %s!",first_n_elements_to_string(missed))
  
  feature_types = feature_types[order(match(feature_types,names(feature_type_to_assay_name)))]
  
  # create Seurat object with first assay
  # include cell and feature metadata
  f = feature_types[1]
  a = feature_type_to_assay_name[f]
  sc = CreateSeuratObject(counts = feature_data[[f]], project = project,assay = a,min.cells = 1,min.features = 0,names.delim = NULL, names.field = NULL,meta.data = metadata_table)
  
  # check that the feature symbols are the same and add feature meta information
  nms = rownames(sc[[a]])
  missed = nms[!nms %in% rownames(features_ids_types)]
  if(length(missed)>0) flog.error("The 'CreateSeuratObject' method modifies feature symbols for assay %s not as expected: %s!",a,first_n_elements_to_string(missed))
  
  sc[[a]] = AddMetaData(sc[[a]],features_ids_types[rownames(sc[[a]]),])
  
  # now add remaining assays
  feature_types = feature_types[-1]
  for(f in feature_types){
    a = feature_type_to_assay_name[f]
    sc[[a]] = CreateAssayObject(counts = feature_data[[f]],min.cells = 0,min.features = 0)
    
    nms = rownames(sc[[a]])
    missed = nms[!nms %in% rownames(features_ids_types)]
    if(length(missed)>0) flog.error("The 'CreateSeuratObject' method modifies feature symbols for assay %s not as expected: %s!",a,first_n_elements_to_string(missed))
    
    sc[[a]] = AddMetaData(sc[[a]],features_ids_types[rownames(sc[[a]]),])
  }
  
  sc
}

#' Saves the assay data of a Seurat object as a sparse matrix dataset in a format similar to that of 10X.
#' 
#' @param object A seurat object.
#' @param dir Name of the directory.
#' @param assays Which assays to export. Default is all but can be specified.
#' @param slot Which slot to use. Default is counts.
#' @param feature_type_to_assay_name How should the feature types for the different assays be named? Default is: "RNA" = "Gene expression","ADT" = "Antibody Capture" ,"HTO" = "Antibody Capture","Crispr" = "CRISPR Guide Capture" and "Custom" = "Custom". Also sets the order in which the feature types are written.
#' @param include_cell_metadata_cols Vector with names of cell metadata columns to be written to the metadata.tsv.gz. Can be NULL in which case no file is created.
#' @param include_feature_metadata_cols Vector with names of feature metadata columns to be written. Can be NULL in which case only the first three columns (feature id, feature symbol and feature type) are written.
#' @return The path to the data directory.
export_seurat_assay_data_to_dir = function(sc,dir="data",assays=NULL,slot="counts",assay_name_to_feature_type = NULL,include_cell_metadata_cols = NULL, include_feature_metadata_cols = NULL){
  library(Seurat)
  library(Matrix)
  library(futile.logger)
  library(R.utils)
  library(dplyr)
  
  # defaults
  if(is.null(assay_name_to_feature_type)) assay_name_to_feature_type = c("RNA" = "Gene Expression","ADT" = "Antibody Capture","HTO" = "Antibody Capture","Crispr" = "CRISPR Guide Capture","Custom" = "Custom")
  
  # make the directory and write data
  d = file.path(dir)
  dir.create(d,showWarnings = F)
  
  # collect assay data and write
  if(is.null(assays)) assays = Assays(sc)
  missed = assays[!assays %in% names(assay_name_to_feature_type)]
  if(length(missed)>0) flog.error("The 'assay_name_to_feature_type' argument misses some assays: %s!",first_n_elements_to_string(missed))
  
  assays = assays[order(match(assays,names(assay_name_to_feature_type)))]
  feature_data = do.call(rbind,lapply(assays,function(a){GetAssayData(sc,assay=a,slot=slot)}))
  mh = file.path(d,"matrix.mtx")
  Matrix::writeMM(feature_data,file=mh)
  gzip(mh,overwrite=T)
  
  # collect cell barcodes and write
  bh = gzfile(file.path(d, "barcodes.tsv.gz"), open="wb")
  write(colnames(feature_data), file=bh)
  close(bh)
  
  # collect feature names (plus additional feature meta data) and write
  # if not set, include columns 1:3
  if(is.null(include_feature_metadata_cols)){
    include_feature_metadata_cols = c(1:3)
  }
  assay_feature_meta_data_df = bind_rows(lapply(assays,function(a){sc[[a]][[include_feature_metadata_cols]]}))
  fh = gzfile(file.path(d, "features.tsv.gz"), open="wb")
  write.table(assay_feature_meta_data_df, file=fh, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
  close(fh)
  
  # collect cell metadata and write
  if(!is.null(include_cell_metadata_cols) & length(include_cell_metadata_cols)>0){
    missed = include_cell_metadata_cols[!include_cell_metadata_cols %in% colnames(sc[[]])]
    if(length(missed)>0) flog.error("The 'include_cell_metadata_cols' argument contains columns which are not in the metadata ob the Seurat object: %s!",first_n_elements_to_string(missed))
    
    metadata_table = sc[[include_cell_metadata_cols]]
    metadata_table_nms = colnames(metadata_table)
    metadata_table$CELL_ID_METADATA = rownames(metadata_table)
    metadata_table = metadata_table[,c("CELL_ID_METADATA",metadata_table_nms)]
    
    metah = gzfile(file.path(d, "metadata.tsv.gz"), open="wb")
    write.table(metadata_table, file = metah, row.names=FALSE, col.names=TRUE, quote=TRUE, sep="\t")
    close(metah)
  }  
  
  return(d)
}