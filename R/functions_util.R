
#' Transformer for the glue package:
#' - if quote=TRUE, quotes all variables in the glue string
#' - if a variable has multiple values and its glue string contains a '*' (e.g. {variable*}, include all values separated by sep. If quote=TRUE, they will be quoted as well.
#' 
#' https://glue.tidyverse.org/articles/transformers.html
#' 
#' @param sep When a variable contains multiple values and is marked with a '*', which separator to use to collapse the values
#' @param quote Whether to quote variables
#' @return A function to transform the glue string
GlueTransformer_quote_collapse = function(sep=", ", quote=TRUE, ...) {
  function(text, envir) {
    collapse = grepl("[*]$", text)
    if (collapse) {
      text = sub("[*]$", "", text)
    }
    res = glue::identity_transformer(text, envir)
    if (quote) {
      res = glue::single_quote(res)
    }
    if (collapse) {
      glue::glue_collapse(res, sep, ...)  
    } else {
      res
    }
  }
}

#' Formats messages. 
#' 
#' - Variables can be automatically inserted with '{variable}'.
#'  - If quote=TRUE, quotes all variables.
#' - If a variable has multiple values, use {variable*} to include all values separated by sep. If quote=TRUE, they will be quoted as well.
#' 
#' @param sep When a variable contains multiple values and is marked with a '*', which separator to use to collapse the values
#' @param quote Whether to quote variables
#' @return The formatted message
FormatMessage = function(msg, quote=TRUE, sep=", ") {
  return(glue::glue(msg, .transformer=GlueTransformer_quote_collapse(), .envir=parent.frame()))
}

#' Returns the content of the profile yaml.
#' 
#' Note: The current profile must be set via the environment variable 'QUARTO_PROFILE'.
#' 
#' @return The content as nested list.
GetProfileYaml = function() {
  profile = Sys.getenv("QUARTO_PROFILE")
  assertthat::assert_that(nchar(profile) > 0,
                          msg=FormatMessage("Environment variable 'QUARTO_PROFILE' must be set to the current profile."))
  
  profile_yml = yaml::read_yaml(paste0("_quarto-", profile, ".yml"), eval.expr=TRUE)
  return(profile_yml)
}

#' Given a module directory, returns the module directory that was run before this module according to the workflow ('chapters') defined in the current profile.
#' 
#' Note: The current profile must be set via the environment variable 'QUARTO_PROFILE'.
#' 
#' @return The previous module directory or NULL if there is none.
PreviousModuleDir = function(current_module_dir) {
  current_module_dir = module_dir
  
  # Get profile yaml parameter 'chapters' in section 'book'
  profile_yml = GetProfileYaml()
  assertthat::assert_that("book" %in% names(profile_yml),
                          msg=FormatMessage("Profile yaml does not contain the parameter 'book'."))
  assertthat::assert_that("chapters" %in% names(profile_yml$book),
                          msg=FormatMessage("Profile yaml does not contain the parameter 'chapters' in the section 'book'."))
  chapters = dirname(profile_yml$book$chapters)
  
  # Find position
  idx = match(current_module_dir, chapters)
  idx = idx[1]
  assertthat::assert_that(!is.na(idx),
                          msg=FormatMessage("Module is not part of parameter 'chapters' in the section 'book'."))
  
  if (idx == 1) {
    return(NULL)
  } else {
    return(chapters[(idx-1)])
  }
}

#' Access a workflow parameters set in the profile and the module yaml.
#' 
#' This function merged parameters defined in the module yaml head general with general and module-specific parameters from the profile yaml (in this order). This is currently the only way to work with profile and module 
#' parameters a) interactively in rstudio as well as b) during rendering by quarto.
#' 
#' Note: The current profile must be set via the environment variable 'QUARTO_PROFILE'.
#' 
#' @param p Parameter to access. If NULL, returns all parameters.
#' @return One or more parameters as list
param = function(p=NULL) {

  # Read module parameter (document params yaml) and get module name
  assertthat::assert_that("module" %in% names(params),
                          msg=FormatMessage("Module does not contain the document yaml parameter 'module' with the module name."))
  module_name = params[["module"]]
  param_set = params
  
  # Read profile params from the profile yaml file
  profile_yml = GetProfileYaml()
  if ("params" %in% names(profile_yml)) {
    profile_params = profile_yml[["params"]]
    
    # Get general parameter if available
    if ("general" %in% names(profile_params)) {
      param_set = purrr::list_modify(param_set, !!!profile_params[["general"]])
    }
    
    # Get module-specific parameter if available
    if ("modules" %in% names(profile_params)) {
      profile_module_params = profile_params[["modules"]]
      
      if (module_name %in% names(profile_module_params)) {
        param_set = purrr::list_modify(param_set, !!!profile_module_params[[module_name]])
      }
    }
  }
  
  if(is.null(p)) {
    return(param_set)
  } else {
    return(param_set[[p]])
  }
}

#' Returns an Ensembl biomaRt for a species and Ensembl version.
#' 
#' @param species Latin species name in format genus_species (for example homo_sapiens or mus_musculus).
#' @param ensembl_version Ensembl version (for example 98).
#' @return A biomaRt object.
GetBiomaRt = function(species, ensembl_version) {
  # Check if we can find find an Ensembl database for this species and annotation version
  ensembl_archives = biomaRt::listEnsemblArchives()
  assertthat::assert_that(
    ensembl_version %in% ensembl_archives$version,
    msg = FormatMessage("Could not find or access Ensembl version {ensembl_version}.")
  )
  
  # Get mart and check if the species is part of ensembl
  ensembl_mart = biomaRt::useMart(biomart = "ensembl", host = ensembl_archives[match(ensembl_version, ensembl_archives$version), "url"])
  ensembl_datasets = biomaRt::listDatasets(ensembl_mart)
  
  if (species == "heterocephalus_glaber") {
    # Hack for heterocephalus_glaber which does not fit in the usual Ensembl way of naming datasets
    species_dataset_name = paste0("hgfemale", "_gene_ensembl")
  } else {
    species_dataset_name = paste0(gsub("^(.)[a-z]+_", "\\1", species), "_gene_ensembl")
  }
  
  idx = which(ensembl_datasets$dataset == species_dataset_name)
  assertthat::assert_that(
    length(idx) == 1,
    msg = FormatMessage(
      "Could not find species {species} dataset (name: {species_dataset_name}) for Ensembl version {ensembl_annotation_version}."
    )
  )
  
  # Get dataset
  ensembl_mart = biomaRt::useDataset(ensembl_datasets[idx[1], "dataset"], ensembl_mart)
  
  return(ensembl_mart)
}

#' Fetches gene information from Ensembl using biomaRt.
#' 
#' @param ids A list of Ensembl ids or - if symbols is TRUE - gene symbols.
#' @param symbols If TRUE, ids are interpreted as gene symbols.
#' @param species Species.
#' @param ensembl_version Ensembl version.
#' @param mart_attributes Ensembl attributes to fetch. Can be a character vector or a named character vector. Defaults to: c(ensembl_id="ensembl_gene_id, ensembl_symbol="external_gene_name", ensembl_biotype="gene_biotype", ensembl_description="description", ensembl_chr="chromosome_name", ensembl_start_position="start_position", ensembl_end_position="end_position", ensembl_strand="strand").
#' @param useCache Use local cache for faster querying. Default is TRUE. Set to FALSE if there are problems.
#' @return A table with gene information. Ids that were not found are included but most of the information will be NA.
EnsemblFetchGeneInfo = function(ids, symbols=FALSE, species, ensembl_version, mart_attributes=c(ensembl_id="ensembl_gene_id", ensembl_symbol="external_gene_name", ensembl_biotype="gene_biotype", ensembl_description="description", ensembl_chr="chromosome_name", ensembl_start_position="start_position", ensembl_end_position="end_position", ensembl_strand="strand"), useCache=TRUE) {
  # Get species mart
  species_mart = GetBiomaRt(species, ensembl_version)
  
  # Check that we have Ensembl ids
  if (!symbols) {
    assertthat::assert_that(any(grepl("^ENS", ids)),
                          msg = FormatMessage(
                            "None of the ids in this dataset is Ensembl. Cannot fetch gene information with this method."))
    id_column = "ensembl_gene_id"
  } else {
    id_column = "external_gene_name"
  }
  
  if (!id_column %in% mart_attributes) {
    mart_attributes = c(setNames(id_column, id_column), mart_attributes)
  }
  
  # Fetch attributes from Ensembl (use cache to allow multiple fetches)
  species_annotation = biomaRt::getBM(mart=species_mart, 
                                      filters=id_column, 
                                      values=ids, 
                                      attributes=mart_attributes, 
                                      useCache=useCache)
  if (!is.null(names(mart_attributes))) {
    colnames(species_annotation) = names(mart_attributes)
  }
  assertthat::assert_that(nrow(species_annotation)>0,
    msg = FormatMessage(
      "Could not find fetch any gene information for this dataset. Is the species {species} correct?"))
  
  # Add rows for ids that were not found
  idx = which(mart_attributes == id_column)
  first_y_col = names(mart_attributes)[idx]
  x_df = data.frame(ids)
  colnames(x_df) = first_y_col
  annotation_ensembl = dplyr::left_join(x_df, species_annotation, by=first_y_col)
  
  return(annotation_ensembl)
}

#' Fetches orthologues information between two species from Ensembl using biomaRt.
#' 
#' @param ids A list of Ensembl ids or - if symbols is TRUE - gene symbols.
#' @param symbols If TRUE, ids are interpreted as gene symbols.
#' @param species1 Species 1.
#' @param species1 Species 2.
#' @param ensembl_version Ensembl version.
#' @param mart_attributes1 Ensembl attributes to fetch fo species 1. Can be a character vector or a named character vector. Defaults to: c(ensembl_id="ensembl_gene_id, ensembl_symbol="external_gene_name"). Can be empty.
#' @param mart_attributes1 Ensembl attributes to fetch fo species 2. Can be a character vector or a named character vector. Defaults to: c(ensembl_id="ensembl_gene_id, ensembl_symbol="external_gene_name"). Cannot be empty.
#' @param useCache Use local cache for faster querying. Default is TRUE. Set to FALSE if there are problems.
#' @return A table with orthologues between species 1 and species 2. Ids that were not found are included but most of the information will be NA.
EnsemblFetchOrthologues = function(ids, symbols=FALSE, species1, species2, ensembl_version, mart_attributes1=c(ensembl_id1="ensembl_gene_id", ensembl_symbol1="external_gene_name"), mart_attributes2=c(ensembl_id2="ensembl_gene_id", ensembl_symbol2="external_gene_name"), useCache=TRUE) {
  # Get species marts
  species1_mart = GetBiomaRt(species1, ensembl_version)
  species2_mart = GetBiomaRt(species2, ensembl_version)
  
  # Check that we have Ensembl ids
  if (!symbols) {
    assertthat::assert_that(any(grepl("^ENS", ids)),
                            msg = FormatMessage(
                              "None of the ids in this dataset is Ensembl. Cannot fetch gene information with this method."))
    id_column = "ensembl_gene_id"
  } else {
    id_column = "external_gene_name"
  }
  
  if (!id_column %in% mart_attributes1) {
    mart_attributes1 = c(setNames(id_column, id_column), mart_attributes1)
  }
  
  ortholog_annotation = biomaRt::getLDS(mart=species1_mart,
                            martL=species2_mart,
                            filters=id_column,
                            values=ids,
                            attributes=mart_attributes1,
                            attributesL=mart_attributes2, 
                            uniqueRows=TRUE)
  if (!is.null(names(mart_attributes1)) & !is.null(names(mart_attributes2))) {
    colnames(ortholog_annotation) = c(names(mart_attributes1), names(mart_attributes2))
  }
  assertthat::assert_that(nrow(ortholog_annotation)>0,
                          msg = FormatMessage(
                            "Could not find fetch any ortholog information for this dataset. Is the species {species1} correct?"))

  
  # Identify one-to-one orthologues
  idx = which(mart_attributes1 == id_column)
  first_y_col = names(mart_attributes1)[idx]
  one_to_one = ortholog_annotation %>% 
    dplyr::group_by(!!sym(first_y_col)) %>%
    dplyr::summarise(n=dplyr::n()) %>%
    dplyr::filter(n==1) %>%
    dplyr::pull(!!sym(first_y_col))
  
  ortholog_annotation = ortholog_annotation %>%
    dplyr::mutate(one_to_one = !!sym(first_y_col) %in% one_to_one)
  
  # Add rows for ids that were not found
  idx = which(mart_attributes1 == id_column)
  first_y_col = names(mart_attributes1)[idx]
  x_df = data.frame(ids)
  colnames(x_df) = first_y_col
  ortholog_annotation = dplyr::left_join(x_df, ortholog_annotation, by=first_y_col)
  
  return(ortholog_annotation)
  
}  

#' Adds feature metadata to an Seurat object or an Assay object.
#' 
#' Note: Seurat::AddMetaData does not seem to work for features in Seurat v5.
#' 
#' @param obj A Seurat or Assay object.
#' @param assay The assay to which to add the feature metadata. Will be ignored if adding to an Assay object.
#' @param metadata A table with feature metadata with feature names being row names.
#' @return The Seurat (v5) or Assay object with updated feature metadata.
AddFeatureMetadata = function(obj, assay=NULL, metadata) {
  # Checks
  valid_objs = c("Seurat", "Assay5", "Assay")
  assertthat::assert_that(class(obj) %in% valid_objs,
                          msg = FormatMessage(
                            "AddFeatureMetadata works only for 'Seurat v5' or 'Assay5' objects"))
  
  # Prepare metadata tables for join
  metadata = metadata %>% tibble::rownames_to_column()
  
  if (is(obj, "Seurat")) {
    feature_metadata = obj[[assay]]@meta.data
  } else if (is(obj, "Assay5") | is(obj, "Assay")) {
    feature_metadata = obj@meta.data
  }
  feature_metadata = feature_metadata %>% tibble::rownames_to_column()
 
  # Join
  feature_metadata = dplyr::left_join(feature_metadata, metadata, by="rowname") %>% as.data.frame()
  rownames(feature_metadata) = feature_metadata$rowname
  feature_metadata = feature_metadata %>% dplyr::select(-rowname)
  
  # Add again
  if (is(obj, "Seurat")) {
    feature_metadata = feature_metadata[rownames(obj[[assay]]), , drop=FALSE]
    obj[[assay]]@meta.data = feature_metadata
  } else if (is(obj, "Assay5") | is(obj, "Assay")) {
    feature_metadata = feature_metadata[rownames(obj), , drop=FALSE]
    obj@meta.data = feature_metadata
  }
  
  return(obj)
}

#' Evaluates a knitr chunk as R code.
#' 
#' @param x Knitr chunk code
#' @return The last return value of the code
EvalKnitrChunk = function(x) {
  chunks = unlist(stringr::str_extract_all(string=x, 
                           pattern=stringr::regex(pattern="```\\s*\\{r\\}.*?\\n```", dotall=TRUE)
                           ))
  chunks = gsub("```", "#", chunks)
  r_code = paste(chunks, collapse="\n")
  r_code = parse(text=r_code)
  return(eval(r_code))
}

#' Prepares the barcode filter.
#' 
#' @param filter Filter from yaml configuration
#' @param orig_idents The samples in the analysis
#' @return A filter with entries for each sample
PrepareBarcodeFilter = function(filter, orig_idents) {
  if (is.null(filter)) {
    filter = rep(list(NULL), length(orig_idents))
    names(filter) = orig_idents
    return(filter)
  }
  
  sample_specific_filter = filter[names(filter) %in% orig_idents]
  general_filter = filter[!names(filter) %in% orig_idents]
  filter = list()
  
  # Apply general filter to all samples
  if (length(general_filter) > 0) {
    filter = rep(list(general_filter), length(orig_idents))
    names(filter) = orig_idents
  }
  
  # Apply sample_specific filter (overwrite or add)
  filter = purrr::list_modify(filter, !!!sample_specific_filter)
  
  return(filter)
}

#' Prepares the feature filter.
#' 
#' @param filter Filter from yaml configuration
#' @param orig_idents The samples in the analysis
#' @return A filter with entries for each sample
PrepareFeatureFilter = function(filter, orig_idents) {
  if (is.null(filter)) {
    filter = list(min_counts=1, min_cells=1)
    filter = rep(filter, length(orig_idents))
    names(filter) = orig_idents
  }
  
  sample_specific_filter = filter[names(filter) %in% orig_idents]
  general_filter = filter[!names(filter) %in% orig_idents]
  filter = list()
  
  # Apply general filter to all samples
  if (length(general_filter) > 0) {
    filter = rep(list(general_filter), length(orig_idents))
    names(filter) = orig_idents
  }
  
  # Apply sample_specific filter (overwrite or add)
  filter = purrr::list_modify(filter, !!!sample_specific_filter)
  
  return(filter)
}

#' Adds one or more lists to the misc slot of the Seurat object.
#' 
#' @param sc A Seurat sc object.
#' @param lists A named list with one or more lists (named or unnamed vectors only).
#' @param lists_slot In which slot of the Seurat misc slot should the list(s) be stored.
#' @param add_to_list When a list with this name already exists, add to the list instead of overwriting the list. Default is FALSE.
#' @param make_unique Make lists unique (after they were stored in the misc slot). Default is FALSE.
#' @return A Seurat sc object with updated list(s).
ScAddLists = function(sc, lists, lists_slot='gene_lists', add_to_list=FALSE, make_unique=FALSE) {
  stored_lists = Seurat::Misc(sc, slot=lists_slot)
  if (is.null(stored_lists)) stored_lists = list()
  
  # Check that the lists have names, otherwise set to 'list_1', 'list_2' etc
  lists_names = purrr::map_chr(seq(lists), function(i) {
    n = names(lists[i])
    if (is.null(n) || nchar(n) == 0) return(paste0("list_",as.character(i))) else return(names(lists[i]))
  })
  names(lists) = lists_names
  
  for(n in names(lists)) {
    if (n %in% names(stored_lists) & add_to_list == TRUE) {
      stored_lists[[n]] = c(stored_lists[[n]], unlist(lists[[n]]))
    } else {
      stored_lists[[n]] = unlist(lists[[n]])
    }
    
    if (make_unique) stored_lists[[n]] = unique(stored_lists[[n]])
  }
  
  # Add to Seurat object
  suppressWarnings({Seurat::Misc(sc, slot=lists_slot) = stored_lists})
  return(sc)
}

#' Get one or more lists from the misc slot of the Seurat object.
#' 
#' @param sc A Seurat sc object.
#' @param lists One or more list names.
#' @param lists_slot From which slot of the Seurat misc slot should the lists be pulled. If NULL, pull from the top level of the Seurat misc slot.
#' @return Lists saved in the misc slot of the Seurat object.
ScLists = function(sc, lists, lists_slot=NULL) {
  stored_lists = Seurat::Misc(sc, slot=lists_slot)
  assertthat::assert_that(!is.null(stored_lists), 
                          msg=FormatMessage("No lists found in misc slot of Seurat object (list slot: {{lists_slot}})."))
  assertthat::assert_that(all(lists %in% names(stored_lists)), 
                          msg=FormatMessage("List(s) {{lists}} not found in misc slot of Seurat object (list slot: {{lists_slot}})."))
  
  if (length(lists) > 1) {
    return(stored_lists[lists])
  } else {
    return(stored_lists[[lists]])
  }
}

#' Adds colours for one or more categories to the misc slot of the Seurat object.
#' 
#' @param sc A Seurat sc object.
#' @param colours A named list with one or more lists with colours for categories (named or unnamed vectors only).
#' @param colours_slot Name of the misc slot which stores the colours. Default is 'colour_lists'.
#' @param add_to_list When a colour list with this name already exists, add to the colour list instead of overwriting the list. Default is FALSE.
#' @param make_unique Make colour lists unique (after they were stored in the misc slot). Default is FALSE.
#' @return A Seurat sc object with updated colour list(s).
ScAddColours = function(sc, colours, colours_slot='colour_lists', add_to_list=FALSE, make_unique=FALSE) {
  return(ScAddLists(sc, colours, lists_slot=colours_slot, add_to_list=add_to_list, make_unique=make_unique))
}

#' Gets colours for one or more categories from the misc slot of the Seurat object.
#' 
#' @param sc A Seurat sc object.
#' @param categories One or more category names.
#' @param colours_slot Name of the misc slot which stores the colours. Default is 'colour_lists'.
#' @return Colour lists for categories.
ScColours = function(sc, categories, colours_slot="colour_lists") {
  return(ScLists(sc, categories, lists_slot=colours_slot))
}

######################
######################

#' Tests if values can be converted to numbers.
#' 
#' @param x A vector with values.
#' @return TRUE if they can be converted otherwise FALSE.
converts_to_number = function(x) {
  return(suppressWarnings(!is.na(as.numeric(na.omit(x)))))
}

#' Tests if values can be converted to logical.
#' 
#' @param x A vector with values.
#' @return TRUE if they can be converted otherwise FALSE.
converts_to_logical = function(x) {
  return(suppressWarnings(!is.na(as.logical(na.omit(x)))))
}


#' Given a vector, report at most n elements as concatenated string.
#' 
#' @param x A vector.
#' @param n Number of elements to report at most.
#' @param sep Separator for string concatenation.
#' @return A string with at most n elements to concatenated.
first_n_elements_to_string = function(x, n=5, sep=",") {
  s = paste(x[1:min(n,length(x))], collapse=sep)
  if (length(x) > n) s = paste(s, "...", sep=sep)
  return(s)
}

#' Get scrnaseq git repository version
#' 
#' @param path_to_git: Path to git repository.
#' @return The git repository version.
GitRepositoryVersion = function(path_to_git) {
  repo = tryCatch({system(paste0("git --git-dir=", path_to_git, "/.git rev-parse HEAD"), intern=TRUE)},
                  warning = function(war) {return("NA")})
  return(repo)
}

#' Get container information (if available)
#' 
#' @return A string with container git name, container git commit id and container build date.
ContainerVersion = function(path_to_git) {
  if (nchar(Sys.getenv("CONTAINER_GIT_NAME")) > 0) {
    container_info = paste(Sys.getenv("CONTAINER_GIT_NAME"), Sys.getenv("CONTAINER_GIT_COMMIT_ID"), Sys.getenv("CONTAINER_BUILD_DATE"), sep=", ")
  } else {
    container_info = c("NA")
  }
  return(container_info)
}

#' Report session info in a table
#' 
#' @param path_to_git: Path to git repository.
#' @return The session info as table.
ScrnaseqSessionInfo = function(path_to_git=".") {
  out = matrix(NA, nrow=0, ncol=2)
  colnames(out) = c("Name", "Value")
  
  # Run on
  out = rbind(out, c("Run on:", format(Sys.time(), "%a %b %d %X %Y")))
  
  # Git
  repo = GitRepositoryVersion(path_to_git)
  out = rbind(out, c("ktrns/scrnaseq", repo))
  
  # Container (if available)
  out = rbind(out, c("Container", ContainerVersion()))
  
  # System
  info_session = sessionInfo()
  out = rbind(out, c("R", info_session$R.version$version.string))
  out = rbind(out, c("Platform", info_session$platform))
  out = rbind(out, c("Operating system", info_session$running))
  out = rbind(out, c("Host name", unname(Sys.info()["nodename"])))
  out = rbind(out, c("Host OS", paste0(unname(Sys.info()["version"]), " (", unname(Sys.info()["release"]), ") ")))
  
  info_pkgs = sessioninfo::package_info()
  out = rbind(out, c("Packages", paste(paste(info_pkgs$package, info_pkgs$loadedversion, sep=""), collapse=", ")))
  
  return(out)
}

#' Subsamples barcodes of a Seurat v5 object.
#' 
#' @param sc A Seurat v5 object.
#' @param n Total number of barcodes to subsample. Default is 500.
#' @param seed Seed for sampling. Default is 1.
#' @param group If not NULL, sample the same number of barcodes from each group defined by this barcode metadata column. The number of barcodes per group is then the total number of barcodes divided by the number of groups. Default is NULL.
#' @return Sampled barcodes.
SubsampleSC = function(sc, n=500, seed=1, group=NULL) {
  barcode_metadata = sc[[]]
  barcodes = rownames(barcode_metadata)
  if (!is.null(group)) {
    barcodes = split(barcodes, barcode_metadata[, group, drop=TRUE])
  } else {
    barcodes = list(barcodes)
  }
  
  n = round(n/length(barcodes))
  
  barcodes = purrr::map(barcodes, function(x) {
    set.seed(seed)
    return(sample(x, min(n, length(x))))
  })
  
  return(unlist(purrr::flatten(barcodes)))
}



#' Returns the names of an object.
#' 
#' @param x A list or vector with names.
#' @return A named vector with names as names and names as values.
list_names = function(x) {
  return(setNames(names(x), names(x)))
}

#' Returns a vector with its values as names.
#' 
#' @param x A vector.
#' @return A vector with its values as names.
values_to_names = function(x) {
  return(setNames(x,x))
}

#' Returns the indices of an object.
#' 
#' @param x A list or vector with names.
#' @return A named vector with names as names and indices as values.
list_indices = function(x) {
  return(setNames(seq(x), names(x)))
}

#' Generate colours based on a palette. If the requested number exceeds the number of colours in the palette, then the palette is reused but with a different alpha.
#' 
#' @param num_colours The number of colours to generate.
#' @param names A character vector with names to be assigned to the colour values. If NULL, no names. 
#' @param palette A palette function for generating the colours.
#' @param palette_options List of additional arguments (beside alpha) to pass on to the palette function.
#' @param alphas Alpha value(s) to use. If the number of colours exceeds the palette, multiple alpha value are used to generate more colours.
#' @return The generated colours.
GenerateColours = function(num_colours, names=NULL, palette="ggsci::pal_igv", alphas=c(1,0.7,0.3), palette_options=list()) {
  palette = tryCatch({eval(parse(text=palette))}, error=function(cond) return(NULL))
  if (is.null(palette)) stop("GenerateColours: Could not find specified palette!")
  
  colours = purrr::flatten_chr(purrr::map(alphas, function(a) {
    palette_options[["alpha"]] = a
    cols = suppressWarnings(do.call(do.call(palette, palette_options), list(100)))
    cols[!is.na(cols)]
  }))
  
  if (num_colours > length(colours)) {
    stop("GenerateColours: Cannot generate the requested number of colours. Please change palette or add alpha values.")
  }
  
  colours = colours[1:num_colours]
  if (!is.null(names)) colours = setNames(colours, names)
  return(colours)
}

#' Returns a message formatted for markdown. 
#' See: https://www.w3schools.com/bootstrap/bootstrap_alerts.asp and https://bookdown.org/yihui/rmarkdown-cookbook/output-hooks.html
#' 
#' @param x The message.
#' @param options Further options.
#' @return The message formatted for markdown.
format_message = function(x, options){
  x = gsub('##', '<br/>', gsub('^## Message:','',x))
  msg = paste(c('\n\n:::{class="alert alert-info alert-dismissible"}',
                '<style> .alert-info { background-color: #abd9c6; color: black; word-wrap: break-word; } </style>', 
                '<a href="#" class="close" data-dismiss="alert" aria-label="close">&times;</a>',
                '<strong>(Message)</strong>',
                x,
                ':::\n\n'), collapse = '\n')
  return(msg)
}

#' Prints a message so that it will be included in the markdown document.
#' Note that cat is used since print will not work.
#' 
#' @param x The message.
#' @param options Further options.
#' @return No return value.
Message = function(x, options){
  # Function asis_output: prints something in mode results="asis"; normal_print: prints something in mode results="markup"
  knitr::asis_output(format_message(x, options))
}

#' Returns a warning message formatted for markdown. 
#' See: https://www.w3schools.com/bootstrap/bootstrap_alerts.asp and https://bookdown.org/yihui/rmarkdown-cookbook/output-hooks.html.
#' @param x The message.
#' @param options Further options.
#' @return The message formatted for markdown.
format_warning = function(x, options){
  x = gsub('##', '<br/>', gsub('^## Warning:','',x))
  warn = paste(c('\n\n:::{class="alert alert-warning alert-dismissible"}',
                 '<style> .alert-warning { background-color: #fae39c; color: black; word-wrap: break-word; } </style>', 
                 '<a href="#" class="close" data-dismiss="alert" aria-label="close">&times;</a>',
                 '<strong>(Warning)</strong>',
                 x,
                 ':::\n\n'), collapse = '\n')
  return(warn)
}

#' Prints a warning so that it will be included in the markdown document.
#' Note that cat is used since print will not work. Will only work with chunk option results="asis".
#' 
#' @param x The message.
#' @param options Further options.
#' @return No return value.
Warning = function(x, options){
  # Function asis_output: prints something in mode results="asis"; normal_print: prints something in mode results="markup".
  knitr::asis_output(format_warning(x, options))
}

#' Returns a error formatted for markdown.
#' See: https://www.w3schools.com/bootstrap/bootstrap_alerts.asp and https://bookdown.org/yihui/rmarkdown-cookbook/output-hooks.html.
#' @param x The message.
#' @param options Further options.
#' @return The message formatted for markdown.
format_error = function(x, options){
  x = gsub('^##','',x)
  x = gsub("Error in eval(expr, envir, enclos):", "", x, fixed = TRUE)
  err = paste(c('\n\n:::{class="alert alert-danger"}',
                 '<style> .alert-danger { background-color: #ffb6c1; color: black; word-wrap: break-word; } </style>', 
                 '<strong>(Error)</strong>',
                 x,
                 ':::\n\n'), collapse = '\n')
  return(err)
}

#' Prints an error so that it will be included in the markdown document.
#' 
#' Note that cat is used since print will not work.
#' @param x The message.
#' @param options Further options.
#' @return No return value.
Error = function(x, options){
  # asis_output: prints something in mode results="asis"; normal_print: prints something in mode results="markup"
  knitr::asis_output(format_error(x, options))
}

#' Report parameters in a table.
#' @param params The parameter list.
#' @return A table with parameters for printing.
ScrnaseqParamsInfo = function(params) { 
  
  # Initialize output table
  out = matrix(NA, nrow=0, ncol=2)
  colnames(out) = c("Name", "Value")
  
  # Convert a (named) vector to a string
  VectorToString = function(x) { 
    if (is.null(names(x))) return(toString(x))
    if (sum(names(x) != "") != length(x)) return(toString(x))
    return(paste(names(x), x, sep="=", collapse=", "))   
  }
  
  # Convert list to table
  for (i in seq(params)) {
    x = params[[i]]
    # List, function or simple vector?
    if(is.list(x)) {
      y = paste(names(x), sapply(names(x), function(j) VectorToString(x[[j]])), sep=":", collapse="; ")
    } else if (is.function(x)) {
      y = "function"
    } else {
      y = VectorToString(x)
    }
    out = rbind(out, c(names(params)[i], y))
  }
  
  return(out)
}

# Checks if the parameters of the scrnaseq workflow are valid and converts them if needed.
#'
#' @param param The parameter list.
#' @return Returns a list with error messages.
check_parameters_scrnaseq = function(param) {
  error_messages = c()
  param[["error_messages"]] = NULL
  
  # Check project_id ###
  if (!"project_id" %in% names(param)) {
    error_messages = c(error_messages, "The parameter 'project_id' is missing!")
  } else {
    param$project_id = as.character(param$project_id)
  }
  
  # Check path_data ###
  if (!"path_data" %in% names(param) | 
      !is.data.frame(param$path_data) | 
      !ncol(param$path_data) >= 4 | 
      !nrow(param$path_data) > 0 | 
      any(!c("name", "type", "path", "stats") %in% colnames(param$path_data))) {
    
    error_messages = c(error_messages, "The parameter 'path_data' needs to be a non-empty data.frame with at least four columns - 'name' (dataset name), 'type' (10x or smartseq2), 'path' (path to counts directory or file), 'stats' (path to 10x metrics summary file; can be NA) - and an optional fifth column 'suffix' (suffix to cell names)")
  } else {
    # name
    param$path_data$name = as.character(param$path_data$name)
    if (any(duplicated(param$path_data$name))) {
      error_messages = c(error_messages, "The column 'name' of the parameter 'path_data' must contain unique values!")
    }
    
    # type
    param$path_data$type = as.character(param$path_data$type)
    if (any(!param$path_data$type %in% c("10x","smartseq2"))) {
      error_messages = c(error_messages, "The column 'type' of the parameter 'path_data' should be either '10x' or 'smartseq2'!")
    }
    
    # path
    param$path_data$path = file.path(param$path_data$path)
    datasets_10x = param$path_data %>% dplyr::filter(type=="10x")
    is_valid = purrr::map_lgl(datasets_10x$path, function(p){
      return(dir.exists(p) & file.exists(file.path(p,"barcodes.tsv.gz")) & file.exists(file.path(p,"features.tsv.gz")) & file.exists(file.path(p,"matrix.mtx.gz")))
    })
    if (length(is_valid) > 0 & any(!is_valid)) {
      error_messages = c(error_messages, "At least one 10x dataset 'path' of the parameter 'path_data' is not a directory or misses at least one of the following files: 'barcodes.tsv', 'features.tsv.gz', 'matrix.mtx.gz'!")
    }
  
    datasets_smartseq2 = param$path_data %>% dplyr::filter(type=="smartseq2")
    is_valid = purrr::map_lgl(datasets_smartseq2$path, file.exists)
    if (length(is_valid) > 0 & any(!is_valid)) {
      error_messages = c(error_messages, "For at least one smartseq2 dataset, the 'path' of the parameter 'path_data' could not be found!")
    }
  
    # stats
    param$path_data$stats = ifelse(is.na(param$path_data$stats), NA, file.path(param$path_data$stats))
    is_valid = purrr::map_lgl(param$path_data$stats, function(p){
      if (is.na(p)) return(TRUE) else return(file.exists(p))
    })
    if (length(is_valid) > 0 && any(!is_valid)) {
      error_messages = c(error_messages, "At least one 'stats' file of the parameter 'path_data' could not be found. If not available, please set to NA!")
    }
    
    # suffix
    if (!"suffix" %in% colnames(param$path_data)) {
      param$path_data$suffix = paste0("-", 1:nrow(param$path_data))
    }
  }
  
  # Check downsample_cells_n ###
  if (("downsample_cells_n" %in% names(param)) & (!is.null(param$downsample_cells_n))) {
    if (!converts_to_number(param$downsample_cells_n)) {
      error_messages = c(error_messages, "The parameter 'path_out' (path for output) is missing!")
    } else {
      param$downsample_cells_n = as.numeric(param$downsample_cells_n)
    }
  } else {
    param["downsample_cells_n"] = list(NULL)
  }
  
  # Check path_out ###
  if (!"path_out" %in% names(param)) {
    error_messages = c(error_messages, "The parameter 'path_out' (path for output) is missing!")
  } else {
    param$path_out = file.path(param$path_out)
  }
  
  # Check file_known_markers ###
  if (!is.null(param$file_known_markers)){
    param$file_known_markers = file.path(param$file_known_markers)
    if (!file.exists(param$file_known_markers)) error_messages = c(error_messages, "The parameter 'file_known_markers' (Excel file with markers) is set but the file cannot be found. If not available, please set to NULL!")
  }
  
  # Check mart_dataset ###
  if (!"mart_dataset" %in% names(param)) {
    error_messages = c(error_messages, "The parameter 'mart_dataset' (Biomart dataset name) is missing!")
  } else {
    param$mart_dataset = as.character(param$mart_dataset)
  }
  
  # Check annot_version ###
  if (!"annot_version" %in% names(param)) {
    error_messages = c(error_messages, "The parameter 'annot_version' (Ensembl version) is missing!")
  } else {
    param$annot_version = as.character(param$annot_version)
  }
  
  # Check annot_main ###
  if (!"annot_main" %in% names(param) |
      any(!c("ensembl", "symbol", "entrez") %in% names(param$annot_main))) {
    error_messages = c(error_messages, "The parameter 'annot_main' is missing or is not a named vector with names 'ensembl', 'symbol' and 'entrez' as well as corresponding values!")
  } else {
    param$annot_main = setNames(as.character(param$annot_main), names(param$annot_main))
  }
  
  # Check file_annot ###
  if (!is.null(param$file_annot)) {
    param$file_annot = file.path(param$file_annot)
  }
  
  # Check mart_attributes ###
  if (!"mart_attributes" %in% names(param) || !is.vector(param$mart_attributes) || length(param$mart_attributes)==0) {
    error_messages = c(error_messages, "The parameter 'mart_attributes' (Biomart attributes) is missing or is not a non-empty vector!")
  } else {
    param$mart_attributes = as.character(param$mart_attributes)
  }
  
  # Check biomart_mirror ###
  if (!is.null(param$biomart_mirror)) {
    if (!param$biomart_mirror %in% c("www", "useast", "uswest", "asia")) error_messages = c(error_messages, "The parameter 'biomart_mirror' (Biomart mirror) is set but does not contain one of the following values: 'www', 'uswest', 'useast' ,'asia'!")
    param$biomart_mirror = as.character(param$biomart_mirror)
  } else {
    param$biomart_mirror = "www"
  }
  
  # Check mt ###
  if (!"mt" %in% names(param)) {
    error_messages = c(error_messages, "The parameter 'mt' (prefix of mitochondrial genes) is missing!")
  } else {
    param$mt = as.character(param$mt)
  }
  
  # Check cell_filter: can be numeric with length 2 or character or factor with any length; can also contain sublists per sample with the same criteria ###
  if ("cell_filter" %in% names(param) && length(param$cell_filter) > 0) {
    is_valid = TRUE
    for (i in seq(param$cell_filter)) {
      f = param$cell_filter[[i]]
      
      if (is.list(f)) {
        # Sample-specific values
        for (j in seq(f)) {
          f_s = f[[j]]
          if (length(f_s) == 2 & all(converts_to_number(f_s))) {
            f_s = as.numeric(f_s)
          } else if ( (is.character(f_s) | is.factor(f_s)) & length(f_s) > 0) {
            f_s = as.character(f_s)
          } else {
            is_valid = FALSE
          }
          f[[j]] = f_s
        }
      } else {
        # Global values
        if (length(f) == 2 & all(converts_to_number(f))) {
          f = as.numeric(f)
        } else if ( (is.character(f) | is.factor(f)) & length(f) > 0) {
          f = as.character(f)
        } else {
          is_valid = FALSE
        }
      }
      param$cell_filter[[i]] = f
    }
    
    if (!is_valid) {
      error_messages = c(error_messages, "The parameter 'cell_filter' should contain filters for cell properties with the following structures: a) lists of length 2 with minimum and maxium for numeric properties (set NA if not applicable/no min/no max), b) character/factor vectors for categorial properties. This can also specified by sample using sublists with the sample names as used by the script!")
    }
  } else {
    param$cell_filter = list()
  }
  
  # Check feature_filter: so far only contains: min_counts and min_cells; can also contain sublists per sample with the same criteria  ###
  if ("feature_filter" %in% names(param) && length(param$feature_filter) > 0) {
    valid = TRUE
    
    # Always set a global default for the minimum number of counts; but sample-specific values will overrule it
    if (!"min_counts" %in% names(param$feature_filter)) {
      param$feature_filter[["min_counts"]] = 1
    } else {
      if (converts_to_number(param$feature_filter[["min_counts"]])) {
        param$feature_filter[["min_counts"]] = as.numeric(param$feature_filter[["min_counts"]])
      } else {
        valid = FALSE
      }
    }
    
    # Always set a global default for the minimum number of cells; but sample-specific values will overrule it
    if (!"min_cells" %in% names(param$feature_filter)) {
      param$feature_filter[["min_cells"]] = 1
    } else {
      if (converts_to_number(param$feature_filter[["min_cells"]])) {
        param$feature_filter[["min_cells"]] = as.numeric(param$feature_filter[["min_cells"]])
      } else {
        valid = FALSE
      }
    }
    
    # Check sample-specific values
    valid = TRUE
    for (n in setdiff(names(param$feature_filter), c("min_counts", "min_cells"))) {
      f = param$feature_filter[[n]]
      
      if ("min_counts" %in% names(f)) {
        if (converts_to_number(f[["min_counts"]])) {
          f[["min_counts"]] = as.numeric(f[["min_counts"]])
        } else {
          valid = FALSE
        }
      }
      
      if ("min_cells" %in% names(f)) {
        if (converts_to_number(f[["min_cells"]])) {
          f[["min_cells"]] = as.numeric(f[["min_cells"]])
        } else {
          valid = FALSE
        }
      }
      
      param$feature_filter[[n]] = f
    }
    
    if (!valid) error_messages =  c(error_messages, paste("The parameter 'feature_filter' can contain: a) 'min_counts' for the minimum counts for a gene to be considered expressed,",
                   "b) 'min_cells' for the minimum number of cells in which a gene must be expressed. This can also specified by sample using sublists with the sample names as used by the script!"))
    
  } else {
    param$feature_filter = list(min_counts = 1, min_cells = 1)
  }
  
  # Check samples_to_drop  ###
  if ("samples_to_drop" %in% names(param)) {
    param$samples_to_drop = as.character(param$samples_to_drop)
  } else {
    param$samples_to_drop = as.character(NULL)
  }
  
  # Check samples_min_cells  ###
  if ("samples_min_cells" %in% names(param)) {
    if (converts_to_number(param$samples_min_cells)) {
      param$samples_min_cells = as.numeric(param$samples_min_cells)
    } else{
      error_messages = c(error_messages, "The parameter 'samples_min_cells' must be a number specifying the minimum number of cells a sample must have. Please set to NULL if there is no minimum!")
    }
  } else {
    param$samples_min_cells = 10
  }
  
  # Check cc_remove  ###
  if ("cc_remove" %in% names(param)) {
    if (converts_to_logical(param$cc_remove)) {
      param$cc_remove = as.logical(param$cc_remove)
    } else {
      error_messages = c(error_messages, "The parameter 'cc_remove' (correct for cell cycle) is missing or is not a logical value!")
    }
  } else {
    param$cc_remove = FALSE
  }

  # Check cc_remove_all  ###
  if ("cc_remove_all" %in% names(param)) {
    if (converts_to_logical(param$cc_remove_all)) {
      param$cc_remove_all = as.logical(param$cc_remove_all)
    } else {
      error_messages = c(error_messages, "The parameter 'cc_remove_all' (remove all cell cycle) is missing or is not a logical value!")
    }
    
    if (param$cc_remove_all && !param$cc_remove) {
      error_messages = c(error_messages, "The parameter 'cc_remove_all' (remove all cell cycle)  cannot be set to TRUE while the parameter 'cc_remove' is set to FALSE!")
    }
    
  } else {
    param$cc_remove_all = FALSE
  }
  
  # Check cc_rescore_after_merge  ###
  if ("cc_rescore_after_merge" %in% names(param)) {
    if (converts_to_logical(param$cc_rescore_after_merge)) {
      param$cc_rescore_after_merge = as.logical(param$cc_rescore_after_merge)
    } else {
      error_messages = c(error_messages, "The parameter 'cc_rescore_after_merge' (rescore after merging/integrating multiple samples) is missing or is not a logical value!")
    }
    
    if (param$cc_remove_all && !param$cc_remove) {
      error_messages = c(error_messages, "The parameter 'cc_rescore_after_merge' (rescore after merging/integrating multiple samples)  cannot be set to TRUE while the parameter 'cc_remove' is set to FALSE!")
    }
  } else {
    param$cc_rescore_after_merge = FALSE
  }
  
  # Check vars_to_regress  ###
  if ("vars_to_regress" %in% names(param)) {
    param$vars_to_regress = as.character(param$vars_to_regress)
  } else {
    param["vars_to_regress"] = NULL
  }
  
  # Check latent_vars  ###
  if ("latent_vars" %in% names(param)) {
    param$latent_vars = as.character(param$latent_vars)
  } else {
    param["latent_vars"] = NULL
  }
  
  # Check integrate_samples
  if ("integrate_samples" %in% names(param) && is.list(param$integrate_samples)) {
    if (!"method" %in% names(param$integrate_samples) || !param$integrate_samples$method %in% c("single", "merge", "integrate")) {
      error_messages = c(error_messages, "The parameter 'integrate_samples' misses a 'method' entry that is one of: 'single' (only one dataset), 'merge' (just merge) or 'integrate' (integrate)!")
    }
    
    if ("method" %in% names(param$integrate_samples) && param$integrate_samples$method=="integrate") {
      if (!"dimensions" %in% names(param$integrate_samples)) {
        error_messages = c(error_messages, "The parameter 'integrate_samples' misses a 'dimensions' entry. Please specify the number of dimensions to include for integration!")
      } else if (!converts_to_number(param$integrate_samples$dimensions)) {
        error_messages = c(error_messages, "The 'dimensions' entry of the parameter 'integrate_samples' must be numeric!")
      } else {
        param$integrate_samples$dimensions = as.numeric(param$integrate_samples$dimensions)
      }
      
      if (!"reference" %in% names(param$integrate_samples)) {
        param$integrate_samples["reference"] = NULL
      }
      
      if (!"use_reciprocal_pca" %in% names(param$integrate_samples)) {
        param$integrate_samples[["use_reciprocal_pca"]] = FALSE
      }
    }
  } else {
    param$integrate_samples = list(method="merge")
  }

  # Check norm (normalisation)
  if ("norm" %in% names(param)) {
    if (!param$norm %in% c("RNA", "SCT")) {
      error_messages = c(error_messages, "The parameter 'norm' (normalisation method to use) must be one of: 'RNA', 'SCT'!")
    }
  } else {
    param$norm = "RNA"
  }

  # Check pc_n
  if ("pc_n" %in% names(param)) {
    if (converts_to_number(param$pc_n)) {
      param$pc_n = as.numeric(param$pc_n)
    } else {
      error_messages = c(error_messages, "The parameter 'pc_n' is not a numeric value!")
    }
  } else {
    param$pc_n = 10
  }
  
  # Check cluster_k
  if ("cluster_k" %in% names(param)) {
    if (converts_to_number(param$cluster_k)) {
      param$cluster_k = as.numeric(param$cluster_k)
    } else {
      error_messages = c(error_messages, "The parameter 'cluster_k' is not a numeric value!")
    }
  } else {
    param$cluster_k = 20
  }
  
  # Check umap_k
  if ("umap_k" %in% names(param)) {
    if (converts_to_number(param$umap_k)) {
      param$umap_k = as.numeric(param$umap_k)
    } else {
      error_messages = c(error_messages, "The parameter 'umap_k' is not a numeric value!")
    }
  } else {
    param$umap_k = 30
  }
  
  # Check cluster_resolution_test
  if ("cluster_resolution_test" %in% names(param)) {
    if (all(purrr::map_lgl(param$cluster_resolution_test, converts_to_number))) {
      param$cluster_resolution_test = as.numeric(param$cluster_resolution_test)
    } else {
      error_messages = c(error_messages, "The parameter 'cluster_resolution_test' does not contain numeric values!")
    }
  } else {
    param$cluster_resolution_test = c()
  }
  
  # Check cluster_resolution
  if ("cluster_resolution_test" %in% names(param)) {
    if (all(purrr::map_lgl(param$cluster_resolution_test, converts_to_number))) {
      param$cluster_resolution_test = as.numeric(param$cluster_resolution_test)
    } else {
      error_messages = c(error_messages, "The parameter 'cluster_resolution_test' does not contain numeric values!")
    }
  } else {
    param$cluster_resolution_test = c()
  }
  
  if ("cluster_resolution" %in% names(param)) {
    if (converts_to_number(param$cluster_resolution)) {
      param$cluster_resolution = as.numeric(param$cluster_resolution)
    } else {
      error_messages = c(error_messages, "The parameter 'cluster_resolution' does not contain a numeric value!")
    }
  } else {
    param$cluster_resolution = 0.5
  }

  # Check marker_padj
  if ("marker_padj" %in% names(param)) {
    if (converts_to_number(param$marker_padj)) {
      param$marker_padj = as.numeric(param$marker_padj)
    } else {
      error_messages = c(error_messages, "The parameter 'marker_padj' is not a numeric value!")
    }
  } else {
    param$marker_padj = 0.05
  }

  # Check marker_log2FC
  if ("marker_log2FC" %in% names(param)) {
    if (converts_to_number(param$marker_log2FC)) {
      param$marker_log2FC = as.numeric(param$marker_log2FC)
    } else {
      error_messages = c(error_messages, "The parameter 'log2fc' is not a numeric value!")
    }
  }
  
  # Check deg_contrasts
  if ("deg_contrasts" %in% names(param) & !is.null(param$deg_contrasts)) {
    if (is.character(param$deg_contrasts) && file.exists(param$deg_contrasts)) {
      deg_contrasts_table = openxlsx::read.xlsx(param$deg_contrasts)
    } else {
      deg_contrasts_table = param$deg_contrasts
    }
    
    if (!is.data.frame(deg_contrasts_table)) {
      error_messages = c(error_messages, "The parameter 'deg_contrasts' must be either a data.frame or an existing Excel file with a valid table in the first sheet!")
    } else {
      if (!all(c("condition_column","condition_group1", "condition_group2") %in% colnames(deg_contrasts_table))){
        error_messages = c(error_messages, "The table specified by parameter 'deg_contrasts' must contain the following columns: 'condition_column', 'condition_group1' and 'condition_group2'!")
      }
      
    }
  }

  # Check enrichr_padj
  if (!"enrichr_padj" %in% names(param)) {
    error_messages = c(error_messages, "The parameter 'p_enrichr' is missing!")
  }

  # Check col
  if (!"col" %in% names(param) || !param$col %in% colors()) {
    error_messages = c(error_messages, "The parameter 'col' is missing or not a valid colour!")
  }

  # Check col_palette_samples
  if (!"col_palette_samples" %in% names(param)) {
    error_messages = c(error_messages, "The parameter 'col_palette_samples' is missing!")
  }
  
  # Check col_palette_clusters
  if (!"col_palette_clusters" %in% names(param)) {
    error_messages = c(error_messages, "The parameter 'col_palette_clusters' is missing!")
  }
  
  # Add to param object and return
  if (length(error_messages) > 0) param[["error_messages"]] = error_messages
  return(param)
  
  # Check
}

# Checks if python is valid.
#'
#' @return Returns a list with error messages.
check_python = function() {
  error_messages = c()
  
  if (!reticulate::py_available(initialize = TRUE) || is.null(reticulate::py_config())) {
    return("Python is not installed on this system or not found in the specified path!")
  }
  
  python_modules = c("leidenalg", "anndata", "scipy")
  is_available = purrr::map_lgl(python_modules, reticulate::py_module_available)
  if (any(!is_available)) {
    error_messages = c(error_messages, paste0("The following python packages are missing: ", paste(python_modules[!is_available], sep=", "),"!"))
  }

  return(error_messages)
}

# Checks if pandoc is valid.
#'
#' @return Returns a list with error messages.
check_pandoc = function() {
  error_messages = c()
  
  if (length(rmarkdown::find_pandoc()) == 0) {
    return("Pandoc is not installed on this system or not found in the specified path!")
  }

  return(error_messages)
}

#' Checks if enrichR is live.
#'
#' @param databases The enrichR databases to use.
#' @return Returns a list with error messages.
check_enrichr = function(databases, site="Enrichr") {
  if(is.null(databases) || length(databases)==0) return(c())
  
  # Is enrichR live
  if (is.null(options("enrichR.live")) || !options("enrichR.live")[[1]]) {
    return("EnrichR is not available or cannot connect to the databases")
  }
  
  # Set Enrichr site
  suppressMessages(enrichR::setEnrichrSite(param$enrichr_site))

  # Are databases available at all
  available_databases = tryCatch({ enrichR::listEnrichrDbs()[,"libraryName"] }, error=function(e) {return(NULL) })
  if (is.null(available_databases)) return("Could not list databases available at enrichR! Please check the enrichR vignette!")
  
  
  # Are the requested databases available
  are_valid = databases %in% available_databases
  if (any(!are_valid)) {
    return(paste0("The following enrichR databases are not available: ",paste(databases[!are_valid],collapse=", "),"!"))
  }
  
  return(c())
}

#' Checks if packages are installed.
#'
#' @param packages A character vector with package names.
#' @return Logical vector with TRUE for installed and FALSE for not installed
packages_installed = function(packages) {
  return(packages %in% installed.packages()[ , "Package"])
}


#' Checks if all packages required for the scrnaseq workflow are installed.
#'
#' @return Returns a list with error messages.
check_installed_packages_scrnaseq = function() {
  required_packages = c("Seurat", "ggplot2", "patchwork", "magrittr",
                        "reticulate", "enrichR", "future", "knitr",
                        "dplyr", "tidyr", "purrr", "stringr", "sctransform", 
                        "Matrix", "kableExtra", "DT", "ggsci",
                        "openxlsx", "readr", "R.utils", "biomaRt",
                        "MAST", "enrichR", "sessioninfo", "cerebroApp",
                        "knitcitations", "sceasy")
  
  is_installed = packages_installed(packages=required_packages)
  if(any(!is_installed)) {
    return(paste0("The R packages '", required_packages[!is_installed],"' are not installed!"))
  } else {
    return(c())
  }
}

#' Checks if Ensembl is available so that annotation and cell cycle markers can be downloaded.
#'
#' @param biomart Biomart database name.
#' @param dataset Dataset name.
#' @param mirror Ensembl mirror.
#' @param version Ensembl version.
#' @param attributes The attributes to download.
#' @param file_annot File with existing annotation to use.
#' @param file_cc_markers File with existing cell cycle markers to use.
#' @return Returns a list with error messages.
check_ensembl = function(biomart, dataset, mirror, version, attributes, file_annot=NULL, file_cc_markers=NULL) {
  error_messages = c()
  
  if (is.null(file_annot) || !file.exists(file_annot) || is.null(file_cc_markers) || !file.exists(file_cc_markers)) {
    # See if mart is available
    annot_mart = suppressWarnings(GetBiomaRt(biomart, dataset, mirror, version))
    if (is.null(annot_mart)) {
      if (is.null(mirror)) {
        return(paste0("Cannot download Ensembl annotation for dataset '",dataset,"', version '",version ,"' using biomaRt!"))
      } else {
        return(paste0("Cannot download Ensembl annotation for dataset '",dataset,"', version '",version ,"' at mirror '",mirror,"' using biomaRt!"))
      }
    }
  
    # See if attributes are valid
    attributes = unique(attributes)
    available_attributes = biomaRt::listAttributes(annot_mart)[,1]
    is_available = attributes %in% available_attributes
    if (any(!is_available)) {
      error_messages = c(error_messages, paste0("The following Ensembl attributes could not be found using biomaRt: ", paste(attributes[!is_available], sep=", "),"!"))
    }
  } else {
    file_annot_tbl = read.delim(file_annot, nrows=50)
    is_available = attributes %in% colnames(file_annot_tbl)
    if (any(!is_available)) {
      error_messages = c(error_messages, paste0("The existing annotation file '", file_annot,"' misses the following Ensembl attributes: ", paste(attributes[!is_available], sep=", "),"!"))
    }
  }
  
  
  return(error_messages)
}

#' On error, R will not start a debugger but just print a traceback. For non-interactive use.
#' @return None.
on_error_just_print_traceback = function(x) {
  options(rlang_trace_top_env = rlang::caller_env())
  options(error = function() {
    sink()
    print(rlang::trace_back(bottom = sys.frame(-1)), simplify = "none")
  })
  return(invisible(NULL))
}

#' On error, R will start a debugger on the terminal. For interactive use without X11.
#' @return None.
on_error_start_terminal_debugger = function(x) {
  options(error = function() {
    sink()
    recover()
  })
  return(invisible(NULL))
}

#' On error, R will run the default debugging process. Default.
#' @return None.
on_error_default_debugging = function(x) {
  return(invisible(NULL))
}

#' Wrapper around citep and citet. Takes care of connection problems.
#'
#' @param reference Argument for citep or citet.
#' @param type Use 'citet' or 'citep'.
#' @return Returns the output of citep or citet.
Cite = function(reference, type="citet") {
  formatted = tryCatch({
    if (type=="citet") knitcitations::citet(reference) else knitcitations::citep(reference)
  },
  error=function(cond) {
    return(NULL)
  })
  
  if (is.null(formatted)) formatted = reference
  return(formatted)
}

