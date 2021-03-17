#' Given a vector, report at most n elements as concatenated string
#' 
#' @param x A vector
#' @param n Number of elements to report at most
#' @param sep Separator for string concatenation
#' @return A string with at most n elements concatenated
FirstElementsToString = function(x, n=5, sep=",") {
  s = paste(x[1:min(n,length(x))], collapse=sep)
  if (length(x) > n) s = paste(s, "...", sep=sep)
  return(s)
}

#' Report session info in a table
#' 
#' @param path_to_git: Path to git repository
#' @return Table containing the session info
SessionInfo = function(path_to_git=".") {
  out = matrix(NA, nrow=0, ncol=2)
  colnames(out) = c("Name", "Version")
  
  repo = tryCatch({system(paste0("git --git-dir ", path_to_git, "/.git log --format='%H' -n 1"), intern=TRUE)},
    warning = function(war) {return("Unknown")})
  out = rbind(out, c("ktrns/scrnaseq", repo))
  
  info_session = sessionInfo()
  out = rbind(out, c("R", info_session$R.version$version.string))
  out = rbind(out, c("Platform", info_session$platform))
  out = rbind(out, c("Operating system", info_session$running))
  
  info_pkgs = sessioninfo::package_info()
  out = rbind(out, c("Packages", paste(paste(info_pkgs$package, info_pkgs$loadedversion, sep=""), collapse=", ")))
  
  return(out)
}

#' Subsample cells
#' 
#' @param sc Seurat sc object
#' @param n Number of cells to subsample
#' @param seed Seed for sampling (default: 1)
#' @param group Metadata column group to consider for sampling with group_proportional
#' @param group_proportional TRUE, if the number of cells sampled from each group should be proportional, FALSE otherwise (default: TRUE; requires group to be defined)
#' @return Names of sampled cells
ScSampleCells = function(sc, n, seed=1, group=NULL, group_proportional=TRUE) {
  set.seed(seed)

  # Set n
  if (n > ncol(sc)) n = ncol(sc)
  
  # Sample cells
  cell_names = sc[[]] %>% tibble::rownames_to_column() %>% dplyr::select(dplyr::all_of(c("rowname", group)))
  colnames(cell_names) = ifelse(is.null(group), "rowname", c("rowname", "group"))
  
  if (!is.null(group) && !group_proportional) {
    num_groups = length(unique(cell_names$group))
    n_per_group = round(n/num_groups)
    
    cell_names_sample = lapply(split(cell_names$rowname, cell_names$group), function(l) {
      if (length(l) < n_per_group) {
        return(sample(l, length(l)))
      } else {
        return(sample(l, n_per_group))
      }
    }) %>% unlist() %>% unname() 
    
  } else {
    cell_names_sample = sample(cell_names$rowname, n)
  }
  
  return(cell_names_sample)
}


#' Take a named vector and return its names as a new vector
#' 
#' @param x A list or vector with names
#' @return A named vector with names of x as names and values
ListNames = function(x) {
  return(setNames(names(x), names(x)))
}

#' Take a vector and return its values as names in a new vector
#' 
#' @param x A vector
#' @return A named vector with values of x as names and values
ValuesToNames = function(x) {
  return(setNames(x,x))
}

#' Take a named vector and return the indices of each element in a new vector
#' 
#' @param x A list or vector with names
#' @return A named vector with names of x as names and indices as values
ListIndices = function(x) {
  return(setNames(seq(x), names(x)))
}

#' Wrapper around the biomaRt::useEnsembl function to cope with unavailable Ensembl mirrors. Tries different Ensembl mirrors and returns a mart object with the mirror that works.
#' 
#' @param biomart A biomaRt database name. Possible database names can be retrieved with the function listEnsembl().
#' @param dataset Dataset you want to use. Possible dataset names can be retrieved with the function listDatasets(mart_obj).
#' @param mirror Specify an Ensembl mirror to connect to. The valid options here are 'www', 'uswest', 'useast', 'asia'. If no mirror is specified, then the first mirror that works will be used. Will be ignored if a host is specified.
#' @param version Ensembl version to connect to when wanting to connect to an archived Ensembl version.
#' @param host Host to connect to. Only needs to be specified if different from www.ensembl.org. 
#' @return A biomaRt object.
GetBiomaRt = function(biomart, dataset, mirror=NULL, version=NULL, host=NULL) {
  
  # Which mirrors to test
  if (is.null(mirror)) {
    mirrors_to_test = c("www", "uswest", "useast", "asia")
  } else {
    mirrors_to_test = c(mirror)
  }
  
  mart_obj = NULL
  if(is.null(host)) {
    # Test and if a mirror is not available, check the next one
    for(m in mirrors_to_test) {
      mart_obj = tryCatch({
        biomaRt::useEnsembl(biomart=biomart, dataset=dataset, mirror=m, version=version)
      },
      error=function(cond) {
        return(NULL)
      })
    
      if(!is.null(mart_obj)) break
    }

  } else {
    # Use specific host
    mart_obj = tryCatch({
      biomaRt::useEnsembl(biomart=biomart, dataset=dataset, host=host, version=version)
    },
    error=function(cond) {
      return(NULL)
    })
  }
  
  return(mart_obj)
}


#' Returns the mirror of a biomaRt object.
#' 
#' @param mart_obj A biomaRt object obtained by GetBiomaRt or useEnsembl name.
#' @return The mirror of the biomaRt object. Can be 'www', 'uswest', 'useast' or 'asia'.
GetBiomaRtMirror = function(mart_obj) {
  mirrors_to_test = c("uswest", "useast", "asia")
  mirror = "www"
  
  for(m in mirrors_to_test){
    if(grepl(pattern=m, x=mart_obj@host)){
      mirror = m
      break
    }
  }
  
  return(mirror)
}

#' Generate colours based on a palette. If the requested number exceeds the number of colours in the palette, then the palette is reused but with a different alpha.
#' 
#' @param num_colours The number of colours to generate.
#' @param palette A palette function for generating the colours.
#' @param palette_options List of additional arguments (beside alpha) to pass on to the palette function.
#' @param alphas Alpha value(s) to use. If the number of colours exceeds the palette, multiple alpha value can be provided to generate more colours.
#' @return The generated colours.
GenerateColours = function(num_colours, palette="ggsci::pal_igv", alphas=c(1,0.7,0.3), palette_options=list()) {
  palette = eval(parse(text=palette))
  colours = purrr::flatten_chr(purrr::map(alphas, function(a) {
    palette_options[["alpha"]] = a
    cols = suppressWarnings(do.call(do.call(palette, palette_options), list(100)))
    cols[!is.na(cols)]
  }))
  
  if (num_colours > length(colours)) {
    stop("GenerateColours: Cannot generate the requested number of colours. Please change palette or add alpha values.")
  }
  
  return(colours[1:num_colours])
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
                '<style> .alert-info { background-color: #abd9c6; color: black; } </style>', 
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
                 '<style> .alert-warning { background-color: #fae39c; color: black; } </style>', 
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
                 '<style> .alert-danger { background-color: #ffb6c1; color: black; } </style>', 
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

# Report parameters in a table.
#' @param params The parameter list.
#' @return A table with parameters for printing.
scrnaseq_params_info = function(params) { 
  
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

# Checks if the parameters are valid.
#'
#' @param params The parameter list.
#' @param 
#' @return Returns a list with error messages.
check_parameters = function(param) {
  error_messages = c()
  
  # Check project_id
  if (!"project_id" %in% names(param)) {
    error_messages = c(error_messages, "The parameter 'project_id' is missing!")
  }
  
  # Check path_data
  if (!"path_data" %in% names(param) | 
      !is.data.frame(param$path_data) | 
      !ncol(param$path_data) == 4 | 
      !nrow(param$path_data) > 0 | 
      any(!c("name", "type", "path", "stats") %in% colnames(param$path_data))) {
    
    error_messages = c(error_messages, "The parameter 'path_data' (table with count data information) is missing or is not a data.frame with three columns and at least one row or misses at least one of the following columns: 'name' (dataset name), 'type' (10x or smartseq2), 'path' (path to counts directory or file), 'stats' (path to 10x metrics summary file; can be NA)!")
  } else {
    if (any(duplicated(param$path_data$name))) {
      error_messages = c(error_messages, "The column 'name' of the parameter 'path_data' must contain unique values!")
    }
    
    if (any(!param$path_data$type %in% c("10x","smartseq2"))) {
      error_messages = c(error_messages, "The column 'type' of the parameter 'path_data' should be either '10x' or 'smartseq2'!")
    }
    
    datasets_10x = param$path_data %>% dplyr::filter(type=="10x")
    is_valid = purrr::map_lgl(datasets_10x$path, function(p){
      return(dir.exists(p) & file.exists(paste(p,"barcodes.tsv.gz", sep="/")) & file.exists(paste(p,"features.tsv.gz", sep="/")) & file.exists(paste(p,"matrix.mtx.gz", sep="/")))
    })
    if (length(is_valid) > 0 & any(!is_valid)) {
      error_messages = c(error_messages, "At least one 10x dataset 'path' of the parameter 'path_data' is not a directory or misses at least one of the following files: 'barcodes.tsv', 'features.tsv.gz', 'matrix.mtx.gz'!")
    }
  
    datasets_smartseq2 = param$path_data %>% dplyr::filter(type=="smartseq2")
    is_valid = purrr::map_lgl(datasets_smartseq2$path, file.exists)
    if (length(is_valid) > 0 & any(!is_valid)) {
      error_messages = c(error_messages, "For at least one smartseq2 dataset, the 'path' of the parameter 'path_data' could not be found!")
    }
  
    is_valid = purrr::map_lgl(param$path_data$stats, function(p){
      if (is.na(p)) return(TRUE) else return(file.exists(p))
    })
    if (length(is_valid) > 0 && any(!is_valid)) {
      error_messages = c(error_messages, "At least one 'stats' file of the parameter 'path_data' could not be found. If not available, please set to NA!")
    }
  }
  # Check path_out
  if (!"path_out" %in% names(param) || !dir.exists(dirname(dirname(dirname("."))))) {
    error_messages = c(error_messages, "The parameter 'path_out' (path for output) is missing or cannot be accessed!")
  }
  
  # Check file_known_markers
  if (!is.null(param$file_known_markers) && !file.exists(param$file_known_markers)) {
    error_messages = c(error_messages, "The parameter 'file_known_markers' (Excel file with markers) is set but the file cannot be found. If not available, please set to NULL!")
  }
  
  # Check mart_dataset
  if (!"mart_dataset" %in% names(param)) {
    error_messages = c(error_messages, "The parameter 'mart_dataset' (Biomart dataset name) is missing!")
  }
  
  # Check annot_version
  if (!"annot_version" %in% names(param)) {
    error_messages = c(error_messages, "The parameter 'annot_version' (Ensembl version) is missing!")
  }
  
  # Check annot_main
  if (!"annot_main" %in% names(param) |
      any(!c("ensembl", "symbol", "entrez") %in% names(param$annot_main))) {
    error_messages = c(error_messages, "The parameter 'annot_main' is missing or is not a named vector with names 'ensembl', 'smybol' and 'entrez' as well as corresponding values!")
  }
  
  # Check file_annot
  if (!is.null(param$file_annot) && !file.exists(param$file_annot)) {
    error_messages = c(error_messages, "The parameter 'file_annot' (annotation file) is set but the file cannot be found. If not available, please set to NULL!")
  }
  
  # Check mart_attributes
  if (!"mart_attributes" %in% names(param) || !is.vector(param$mart_attributes) || length(param$mart_attributes)==0) {
    error_messages = c(error_messages, "The parameter 'mart_attributes' (Biomart attributes) is missing or is not a non-empty vector!")
  }
  
  # Check biomart_mirror
  if (!is.null(param$biomart_mirror) && !param$biomart_mirror %in% c("www", "useast", "uswest", "asia")) {
    error_messages = c(error_messages, "The parameter 'biomart_mirror' (Biomart mirror) is set but does not contain one of the following values: 'www', 'uswest', 'useast' ,'asia'!")
  }
  
  # Check mt
  if (!"mt" %in% names(param)) {
    error_messages = c(error_messages, "The parameter 'mt' (prefix of mitochondrial genes) is missing!")
  }
  
  # Check cell_filter: can be numeric with length 2 or character or factor with any length; can also contain sublists per sample with the same criteria
  if ("cell_filter" %in% names(param) && length(param$cell_filter) > 0) {
    is_valid = purrr::map_lgl(param$cell_filter, function(f) {
      if (is.list(f)) {
        return(any(!purrr::map_lgl(f, function(g){
          return((is.numeric(g) & length(g)==2) | is.character(f) | is.factor(f))
        })))
      } else {
        return((is.numeric(f) & length(f)==2) | is.character(f) | is.factor(f))
      }
    })
    if (any(!is_valid)) {
      error_messages = c(error_messages, "The parameter 'cell_filter' should contain filters for cell properties with the following structures: a) lists of length 2 with minimum and maxium for numeric properties (set NA if not applicable), b) character/factor vectors for categorial properties. This can also specified by sample using sublists with the sample names as used by the script!")
    }
  }
  
  # Check feature_filter: so far only contains: min_counts and min_cells; can also contain sublists per sample with the same criteria
  if (!"feature_filter" %in% names(param) || !is.list(param$feature_filter)) {
    error_messages = c(error_messages, "The parameter 'feature_filter' is missing or not a list!")
  } else {
    if ("feature_filter" %in% names(param) && length(param$feature_filter) > 0) {
      is_valid = sum(c("min_counts", "min_cells") %in% names(param$feature_filter))==2
      if (!is_valid) {
        is_valid = purrr::map_lgl(param$feature_filter, function(f) {
          return(sum(c("min_count", "min_cells") %in% names(f))==2)
        })
      }
      if (any(!is_valid)) {
        error_messages = c(error_messages, "The parameter 'feature_filter' should contain global entries for 'min_counts' (minimum count for gene to be expressed) and 'min_cells' (minimum cells for a gene to be considered). This can also specified by sample using sublists with the sample names as used by the script!")
      }
    }
  }
  
  # Check samples_to_drop
  if ("samples_to_drop" %in% names(param) && !is.character(param$samples_to_drop)) {
    error_messages = c(error_messages, "The parameter 'samples_to_drop' is not a character vector. Please set to NULL if there are no samples to drop!")
  }
  
  # Check samples_min_cells
  if ("samples_min_cells" %in% names(param) && !is.numeric(param$samples_min_cells)) {
    error_messages = c(error_messages, "The parameter 'samples_min_cells' must be a number specifying the minimum number of cells a sample must have. Please set to NULL if there is no minimum!")
  }
  
  # Check cc_remove
  if (!"cc_remove" %in% names(param) || !is.logical(param$cc_remove)) {
    error_messages = c(error_messages, "The parameter 'cc_remove' (correct for cell cycle) is missing or is not a logical value!")
  }

  # Check cc_remove_all
  if (!"cc_remove_all" %in% names(param) || !is.logical(param$cc_remove_all)) {
    error_messages = c(error_messages, "The parameter 'cc_remove_all' (remove all cell cycle) is missing or is not a logical value!")
  }
  
  # Check vars_to_regress
  if ("vars_to_regress" %in% names(param) && !is.character(param$vars_to_regress)) {
    error_messages = c(error_messages, "The parameter 'vars_to_regress' is not a character vector. Please set to NULL if there are no variables to regress out!")
  }
  
  # Check latent_vars
  if ("latent_vars" %in% names(param) && !is.character(param$latent_vars)) {
    error_messages = c(error_messages, "The parameter 'latent_vars' is not a character vector. Please set to NULL if there are no variables to included in statistical tests!")
  }
  
  # Check integrate_samples
  if (!"integrate_samples" %in% names(param) || !is.list(param$integrate_samples)) {
    error_messages = c(error_messages, "The parameter 'integrate_samples' is missing or not a list!")
  } else {
    if (!"method" %in% names(param$integrate_samples) || !param$integrate_samples$method %in% c("single", "merge", "standard", "reference", "reciprocal")) {
      error_messages = c(error_messages, "The sub parameter 'method' of parameter 'integrate_samples' is missing or not one of: 'single' (only one dataset), 'merge' (just merge), 'standard' (integrate), 'reference' (integrate based on reference), 'reciprocal' (fast integrate in PCA space)!")
    }
    
    if ("method" %in% names(param$integrate_samples) && param$integrate_samples$method %in% c("standard", "reference", "reciprocal") && !"dimensions" %in% names(param$integrate_samples)) {
      error_messages = c(error_messages, "The sub parameter 'dimensions' of parameter 'integrate_samples' is missing. Please specify the number of dimensions to include for integration!")
    }
    
    if ("method" %in% names(param$integrate_samples) && param$integrate_samples$method %in% c("reference") && !"reference" %in% names(param$integrate_samples)) {
      error_messages = c(error_messages, "The sub parameter 'reference' of parameter 'integrate_samples' is missing. Please specify the reference sample to use for integration!")
    }

    if ("method" %in% names(param$integrate_samples) && param$integrate_samples$method %in% c("single") && nrow(param$path_data) > 1) {
      error_messages = c(error_messages, "The sub parameter 'method' of parameter 'integrate_samples' cannot be 'single' when there are multiple datasets! Please set to another method. If after filtering only one dataset remains, the script will automatically set the 'method' parameter to 'single'.")
    }
    
  }

  # Check normalisation_default
  if (!"norm" %in% names(param) || !param$norm %in% c("RNA", "SCT")) {
    error_messages = c(error_messages, "The parameter 'norm' (normalisation method to use) is missing or not one of: 'RNA', 'SCT'!")
  }

  # Check pc_n
  if (!"pc_n" %in% names(param)) {
    error_messages = c(error_messages, "The parameter 'pc_n' is missing!")
  }

  # Check cluster_resolution
  if (!"cluster_resolution" %in% names(param)) {
    error_messages = c(error_messages, "The parameter 'cluster_resolution' is missing!")
  }

  # Check marker_padj
  if (!"marker_padj" %in% names(param)) {
    error_messages = c(error_messages, "The parameter 'padj' is missing!")
  }

  # Check marker_log2FC
  if (!"marker_log2FC" %in% names(param)) {
    error_messages = c(error_messages, "The parameter 'log2fc' is missing!")
  }
  
  # Check deg_contrasts
  if ("deg_contrasts" %in% names(param) & !is.null(param$deg_contrasts)) {
    if (is.character(param$deg_contrasts) && file.exists(param$deg_contrasts)) {
      deg_contrasts_table = openxlsx::read.xlsx(param$deg_contrasts)
    } else {
      deg_contrasts_table = param$deg_contrasts
    }
    
    if (!is.data.frame(param$deg_contrasts)) {
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
  
  return(error_messages)
}

# Checks if python is valid.
#'
#' @return Returns a list with error messages.
check_python = function() {
  error_messages = c()
  
  if (!reticulate::py_available(initialize = TRUE) || is.null(reticulate::py_config())) {
    return("Python is not installed on this system or not found in the specified path!")
  }
  
  if (!reticulate::py_module_available("leidenalg")) {
    error_messages = c(error_messages, "The python package 'leidenalg' is missing!")
  }

  return(error_messages)
}

# Checks if enrichR is live.
#'
#' @param databases The enrichR databases to use.
#' @return Returns a list with error messages.
check_enrichr = function(databases) {
  if(is.null(databases) || length(databases)==0) return(c())
  
  # Is enrichR live
  if (is.null(options("enrichR.live")) || !options("enrichR.live")[[1]]) {
    return("EnrichR is not available or cannot connect to the databases")
  }
  
  # Are databases available at all
  available_databases = tryCatch({ enrichR::listEnrichrDbs()[,"libraryName"] }, error=function(e) {return(NULL) })
  if (is.null(available_databases)) return("Could not list databases available at enrichR! Please check the enrichR vignette!")
  
  
  # Are the requested databases available
  are_valid = databases %in% available_databases
  if (any(!are_valid)) {
    return(paste("The following enrichR databases are not available:",paste(databases[!are_valid],collapse=", "),"!"))
  }
  
  return(c())
}

# Checks if all required packages are installed.
#'
#' @return Returns a list with error messages.
check_installed_packages = function() {
  required_packages = c("Seurat", "ggplot2", "patchwork", "magrittr",
                        "reticulate", "enrichR", "future", "knitr",
                        "dplyr", "tidyr", "purrr", "stringr", "sctransform", 
                        "Matrix", "kableExtra", "DT", "ggsci", "ggpubr",
                        "openxlsx", "readr", "R.utils", "biomaRt",
                        "MAST", "enrichR", "sessioninfo", "cerebroApp",
                        "knitcitations")
  
  is_installed = required_packages %in% installed.packages()[,"Package"]
  
  if(any(!is_installed)) {
    return(paste0("The R package '", required_packages[!is_installed],"' is not installed!"))
  } else {
    return(c())
  }
}

# Checks if Ensembl annotation is available.
#'
#' @param biomart Biomart database name.
#' @param dataset Dataset name.
#' @param mirror Ensembl mirror.
#' @param version Ensembl version.
#' @param attributes The attributes to download.
#' @return Returns a list with error messages.
check_ensembl = function(biomart, dataset, mirror, version, attributes) {
  error_messages = c()

  # See if mart is available
  annot_mart = suppressWarnings(GetBiomaRt(biomart, dataset, mirror, version))
  if (is.null(annot_mart)) {
    if (is.null(mirror)) {
      return(paste0("Cannot download Ensembl annotation for dataset '",dataset,"', version '",version ," using biomaRt'!"))
    } else {
      return(paste0("Cannot download Ensembl annotation for dataset '",dataset,"', version '",version ,"' at mirror '",mirror,"' using biomaRt!"))
    }
  }
  
  # See if attributes are valid
  attributes = unique(attributes)
  available_attributes = biomaRt::listAttributes(annot_mart)[,1]
  is_available = attributes %in% available_attributes
  if (any(!is_available)) {
    error_messages = c(error_messages, paste("The following Ensembl attributes could not be found using biomaRt:",paste(attributes[!is_available], sep=", "),"!"))
  }

  return(error_messages)
}

#' On error, R will not start a debugger but just print a traceback. For non-interactive use.
#' @return None.
on_error_just_print_traceback = function(x) {
  options(rlang_trace_top_env = rlang::current_env())
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

# Wrapper around citep and citet. Takes care of connection problems.
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
  
  if (is.null(formatted)) formatted = paste0("'", reference, "' (Citation server error)")
  return(formatted)
}
