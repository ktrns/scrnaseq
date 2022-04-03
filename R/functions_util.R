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

#' Returns a subsample of cells.
#' 
#' @param sc A Seurat sc object.
#' @param n Number of cells to subsample.
#' @param seed Seed for sampling. Default is 1.
#' @param group Metadata column group to consider for sampling with group_proportional.
#' @param group_proportional Should the number of cells sampled from each group be proportional (TRUE) or should the number of cells be the same for each group (FALSE)? Only works if group is not NULL.
#' @return Sampled cell names.
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

#' Adds one more lists to the misc slot of the Seurat object.
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
    stop("GenerateColours: Cannot generate the requested number of colours since the palette does not provide enough colours. Please change palette or use alpha values.")
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
#' @param parameter The parameter list used by the scrnaseq script.
#' @return Returns a list with error messages.
check_parameters_scrnaseq = function(parameter) {
  error_messages = c()

  # Check author ###
  if (!"author" %in% names(parameter)) {
    error_messages = c(error_messages, "The parameter 'author' is missing!")
  } else {
    parameter$author = as.character(parameter$author)
  }
  
  # Check project_id ###
  if (!"project_id" %in% names(parameter)) {
    error_messages = c(error_messages, "The parameter 'project_id' is missing!")
  } else {
    parameter$project_id = as.character(parameter$project_id)
  }
  
  # Check path_data ###
  if (!"path_data" %in% names(parameter) | 
      !is.data.frame(parameter$path_data) | 
      !ncol(parameter$path_data) >= 4 | 
      !nrow(parameter$path_data) > 0 | 
      any(!c("name", "type", "path", "stats") %in% colnames(parameter$path_data))) {
    
    error_messages = c(error_messages, "The parameter 'path_data' needs to be a non-empty data.frame with at least four columns - 'name' (dataset name), 'type' (10x or smartseq2), 'path' (path to counts directory or file), 'stats' (path to 10x metrics summary file; can be NA) - and an optional fifth column 'suffix' (suffix to cell names)")
  } else {
    # name
    parameter$path_data$name = as.character(parameter$path_data$name)
    if (any(duplicated(parameter$path_data$name))) {
      error_messages = c(error_messages, "The column 'name' of the parameter 'path_data' must contain unique values!")
    }
    
    # type
    parameter$path_data$type = as.character(parameter$path_data$type)
    if (any(!parameter$path_data$type %in% c("10x","smartseq2"))) {
      error_messages = c(error_messages, "The column 'type' of the parameter 'path_data' should be either '10x' or 'smartseq2'!")
    }
    
    # path
    parameter$path_data$path = file.path(parameter$path_data$path)
    datasets_10x = parameter$path_data %>% dplyr::filter(type=="10x")
    is_valid = purrr::map_lgl(datasets_10x$path, function(p){
      return(dir.exists(p) & file.exists(file.path(p,"barcodes.tsv.gz")) & file.exists(file.path(p,"features.tsv.gz")) & file.exists(file.path(p,"matrix.mtx.gz")))
    })
    if (length(is_valid) > 0 & any(!is_valid)) {
      error_messages = c(error_messages, "At least one 10x dataset 'path' of the parameter 'path_data' is not a directory or misses at least one of the following files: 'barcodes.tsv', 'features.tsv.gz', 'matrix.mtx.gz'!")
    }
  
    datasets_smartseq2 = parameter$path_data %>% dplyr::filter(type=="smartseq2")
    is_valid = purrr::map_lgl(datasets_smartseq2$path, file.exists)
    if (length(is_valid) > 0 & any(!is_valid)) {
      error_messages = c(error_messages, "For at least one smartseq2 dataset, the 'path' of the parameter 'path_data' could not be found!")
    }
  
    # stats
    parameter$path_data$stats = ifelse(is.na(parameter$path_data$stats), NA, file.path(parameter$path_data$stats))
    is_valid = purrr::map_lgl(parameter$path_data$stats, function(p){
      if (is.na(p)) return(TRUE) else return(file.exists(p))
    })
    if (length(is_valid) > 0 && any(!is_valid)) {
      error_messages = c(error_messages, "At least one 'stats' file of the parameter 'path_data' could not be found. If not available, please set to NA!")
    }
    
    # suffix
    if (!"suffix" %in% colnames(parameter$path_data)) {
      parameter$path_data$suffix = paste0("-", 1:nrow(parameter$path_data))
    }
  }
  
  # Check assay_raw (type of input)
  if ("assay_raw" %in% names(parameter)) {
    if (!parameter$assay_raw %in% c("RNA", "Spatial")) {
      error_messages = c(error_messages, "The parameter 'assay_raw' (type of input) must be one of: 'RNA', 'Spatial'!")
    }
  } else {
    parameter$norm = "RNA"
  }
  
  # Check downsample_cells_n ###
  if (("downsample_cells_n" %in% names(parameter)) & (!is.null(parameter$downsample_cells_n))) {
    if (!converts_to_number(parameter$downsample_cells_n)) {
      error_messages = c(error_messages, "The parameter 'downsample_cells_n' is not a number!")
    } else {
      parameter$downsample_cells_n = as.numeric(parameter$downsample_cells_n)
    }
  } else {
    parameter["downsample_cells_n"] = list(NULL)
  }
  
  # Check path_out ###
  if (!"path_out" %in% names(parameter)) {
    error_messages = c(error_messages, "The parameter 'path_out' (path for output) is missing!")
  } else {
    parameter$path_out = file.path(parameter$path_out)
  }
  
  # Check overwrite ###
  if ("overwrite" %in% names(parameter)) {
    if (converts_to_logical(parameter$overwrite)) {
      parameter$overwrite = as.logical(parameter$overwrite)
    } else {
      error_messages = c(error_messages, "The parameter 'overwrite' is not a logical value!")
    }
  } else {
    parameter$overwrite = FALSE
  }
  
  # Check file_known_markers ###
  if (!is.null(parameter$file_known_markers)){
    parameter$file_known_markers = file.path(parameter$file_known_markers)
    if (!file.exists(parameter$file_known_markers)) error_messages = c(error_messages, "The parameter 'file_known_markers' (Excel file with markers) is set but the file cannot be found. If not available, please set to NULL!")
  } else {
    parameter["file_known_markers"] = list(NULL)
  }
  
  # Check mart_dataset ###
  if (!"mart_dataset" %in% names(parameter)) {
    error_messages = c(error_messages, "The parameter 'mart_dataset' (Biomart dataset name) is missing!")
  } else {
    parameter$mart_dataset = as.character(parameter$mart_dataset)
  }
  
  # Check annot_version ###
  if (!"annot_version" %in% names(parameter)) {
    error_messages = c(error_messages, "The parameter 'annot_version' (Ensembl version) is missing!")
  } else {
    parameter$annot_version = as.character(parameter$annot_version)
  }
  
  # Check annot_main ###
  if (!"annot_main" %in% names(parameter) |
      any(!c("ensembl", "symbol", "entrez") %in% names(parameter$annot_main))) {
    error_messages = c(error_messages, "The parameter 'annot_main' is missing or is not a named vector with names 'ensembl', 'symbol' and 'entrez' as well as corresponding values!")
  } else {
    parameter$annot_main = setNames(as.character(parameter$annot_main), names(parameter$annot_main))
  }
  
  # Check file_annot (no check done) ###
  parameter$file_annot = file.path(parameter$file_annot)
  
  # Check mart_attributes ###
  if (!"mart_attributes" %in% names(parameter) || !is.vector(parameter$mart_attributes) || length(parameter$mart_attributes)==0) {
    error_messages = c(error_messages, "The parameter 'mart_attributes' (Biomart attributes) is missing or is not a non-empty vector!")
  } else {
    parameter$mart_attributes = as.character(parameter$mart_attributes)
  }
  
  # Check biomart_mirror ###
  if (!is.null(parameter$biomart_mirror)) {
    if (!parameter$biomart_mirror %in% c("www", "useast", "uswest", "asia")) error_messages = c(error_messages, "The parameter 'biomart_mirror' (Biomart mirror) is set but does not contain one of the following values: 'www', 'uswest', 'useast' ,'asia'!")
    parameter$biomart_mirror = as.character(parameter$biomart_mirror)
  } else {
    parameter$biomart_mirror = "www"
  }
  
  # Check mt ###
  if (!"mt" %in% names(parameter)) {
    error_messages = c(error_messages, "The parameter 'mt' (prefix of mitochondrial genes) is missing!")
  } else {
    parameter$mt = as.character(parameter$mt)
  }
  
  # Check cell_filter: can be numeric with length 2 or character or factor with any length; can also contain sublists per sample with the same criteria ###
  if ("cell_filter" %in% names(parameter) && length(parameter$cell_filter) > 0) {
    is_valid = TRUE
    for (i in seq(parameter$cell_filter)) {
      f = parameter$cell_filter[[i]]
      
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
      parameter$cell_filter[[i]] = f
    }
    
    if (!is_valid) {
      error_messages = c(error_messages, "The parameter 'cell_filter' should contain filters for cell properties with the following structures: a) lists of length 2 with minimum and maxium for numeric properties (set NA if not applicable/no min/no max), b) character/factor vectors for categorial properties. This can also specified by sample using sublists with the sample names as used by the script!")
    }
  } else {
    parameter$cell_filter = list()
  }
  
  # Check feature_filter: so far only contains: min_counts and min_cells; can also contain sublists per sample with the same criteria  ###
  if ("feature_filter" %in% names(parameter) && length(parameter$feature_filter) > 0) {
    valid = TRUE
    
    # Always set a global default for the minimum number of counts; but sample-specific values will overrule it
    if (!"min_counts" %in% names(parameter$feature_filter)) {
      parameter$feature_filter[["min_counts"]] = 1
    } else {
      if (converts_to_number(parameter$feature_filter[["min_counts"]])) {
        parameter$feature_filter[["min_counts"]] = as.numeric(parameter$feature_filter[["min_counts"]])
      } else {
        valid = FALSE
      }
    }
    
    # Always set a global default for the minimum number of cells; but sample-specific values will overrule it
    if (!"min_cells" %in% names(parameter$feature_filter)) {
      parameter$feature_filter[["min_cells"]] = 1
    } else {
      if (converts_to_number(parameter$feature_filter[["min_cells"]])) {
        parameter$feature_filter[["min_cells"]] = as.numeric(parameter$feature_filter[["min_cells"]])
      } else {
        valid = FALSE
      }
    }
    
    # Check sample-specific values
    valid = TRUE
    for (n in setdiff(names(parameter$feature_filter), c("min_counts", "min_cells"))) {
      f = parameter$feature_filter[[n]]
      
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
      
      parameter$feature_filter[[n]] = f
    }
    
    if (!valid) error_messages =  c(error_messages, paste("The parameter 'feature_filter' can contain: a) 'min_counts' for the minimum counts for a gene to be considered expressed,",
                   "b) 'min_cells' for the minimum number of cells in which a gene must be expressed. This can also specified by sample using sublists with the sample names as used by the script!"))
    
  } else {
    parameter$feature_filter = list(min_counts = 1, min_cells = 1)
  }
  
  # Check samples_to_drop  ###
  if ("samples_to_drop" %in% names(parameter) & !is.null(parameter$samples_to_drop)) {
    parameter$samples_to_drop = as.character(parameter$samples_to_drop)
  } else {
    parameter["samples_to_drop"] = list(NULL)
  }
  
  # Check samples_min_cells  ###
  if ("samples_min_cells" %in% names(parameter)) {
    if (converts_to_number(parameter$samples_min_cells)) {
      parameter$samples_min_cells = as.numeric(parameter$samples_min_cells)
    } else{
      error_messages = c(error_messages, "The parameter 'samples_min_cells' must be a number specifying the minimum number of cells a sample must have. Please set to NULL if there is no minimum!")
    }
  } else {
    parameter$samples_min_cells = 10
  }
  
  # Check norm (normalisation)  ###
  if ("norm" %in% names(parameter)) {
    if (!parameter$norm %in% c("RNA", "SCT")) {
      error_messages = c(error_messages, "The parameter 'norm' (normalisation method to use) must be one of: 'RNA', 'SCT'!")
    }
  } else {
    parameter$norm = "RNA"
  }
  
  # Check cc_remove  ###
  if ("cc_remove" %in% names(parameter)) {
    if (converts_to_logical(parameter$cc_remove)) {
      parameter$cc_remove = as.logical(parameter$cc_remove)
    } else {
      error_messages = c(error_messages, "The parameter 'cc_remove' (correct for cell cycle) is missing or is not a logical value!")
    }
  } else {
    parameter$cc_remove = FALSE
  }

  # Check cc_remove_all  ###
  if ("cc_remove_all" %in% names(parameter)) {
    if (converts_to_logical(parameter$cc_remove_all)) {
      parameter$cc_remove_all = as.logical(parameter$cc_remove_all)
    } else {
      error_messages = c(error_messages, "The parameter 'cc_remove_all' (remove all cell cycle) is missing or is not a logical value!")
    }
    
    if (parameter$cc_remove_all && !parameter$cc_remove) {
      error_messages = c(error_messages, "The parameter 'cc_remove_all' (remove all cell cycle)  cannot be set to TRUE while the parameter 'cc_remove' is set to FALSE!")
    }
    
  } else {
    parameter$cc_remove_all = FALSE
  }
  
  # Check cc_rescore_after_merge  ###
  if ("cc_rescore_after_merge" %in% names(parameter)) {
    if (converts_to_logical(parameter$cc_rescore_after_merge)) {
      parameter$cc_rescore_after_merge = as.logical(parameter$cc_rescore_after_merge)
    } else {
      error_messages = c(error_messages, "The parameter 'cc_rescore_after_merge' (rescore after merging/integrating multiple samples) is missing or is not a logical value!")
    }
    
    if (parameter$cc_remove_all && !parameter$cc_remove) {
      error_messages = c(error_messages, "The parameter 'cc_rescore_after_merge' (rescore after merging/integrating multiple samples)  cannot be set to TRUE while the parameter 'cc_remove' is set to FALSE!")
    }
  } else {
    parameter$cc_rescore_after_merge = FALSE
  }
  
  # Check vars_to_regress  ###
  if ("vars_to_regress" %in% names(parameter) & !is.null(parameter$samples_to_drop)) {
    parameter$vars_to_regress = as.character(parameter$vars_to_regress)
  } else {
    parameter["vars_to_regress"] = list(NULL)
  }
  
  # Check latent_vars  ###
  if ("latent_vars" %in% names(parameter) & !is.null(parameter$samples_to_drop)) {
    parameter$latent_vars = as.character(parameter$latent_vars)
  } else {
    parameter["latent_vars"] = list(NULL)
  }
  
  # Check integrate_samples  ###
  if ("integrate_samples" %in% names(parameter) && is.list(parameter$integrate_samples)) {
    if (!"method" %in% names(parameter$integrate_samples) || !parameter$integrate_samples$method %in% c("single", "merge", "integrate")) {
      error_messages = c(error_messages, "The parameter 'integrate_samples' misses a 'method' entry that is one of: 'single' (only one dataset), 'merge' (just merge) or 'integrate' (integrate)!")
    }
    
    if ("method" %in% names(parameter$integrate_samples) && parameter$integrate_samples$method=="integrate") {
      if (!"dimensions" %in% names(parameter$integrate_samples)) {
        error_messages = c(error_messages, "The parameter 'integrate_samples' misses a 'dimensions' entry. Please specify the number of dimensions to include for integration!")
      } else if (!converts_to_number(parameter$integrate_samples$dimensions)) {
        error_messages = c(error_messages, "The 'dimensions' entry of the parameter 'integrate_samples' must be numeric!")
      } else {
        parameter$integrate_samples$dimensions = as.numeric(parameter$integrate_samples$dimensions)
      }
      
      if (!"reference" %in% names(parameter$integrate_samples)) {
        parameter$integrate_samples["reference"] = NULL
      }
      
      if (!"use_reciprocal_pca" %in% names(parameter$integrate_samples)) {
        parameter$integrate_samples[["use_reciprocal_pca"]] = FALSE
      }
    }
  } else {
    parameter$integrate_samples = list(method="merge")
  }

  # Check pc_n  ###
  if ("pc_n" %in% names(parameter)) {
    if (converts_to_number(parameter$pc_n)) {
      parameter$pc_n = as.numeric(parameter$pc_n)
    } else {
      error_messages = c(error_messages, "The parameter 'pc_n' is not a numeric value!")
    }
  } else {
    parameter$pc_n = 10
  }
  
  # Check cluster_k  ###
  if ("cluster_k" %in% names(parameter)) {
    if (converts_to_number(parameter$cluster_k)) {
      parameter$cluster_k = as.numeric(parameter$cluster_k)
    } else {
      error_messages = c(error_messages, "The parameter 'cluster_k' is not a numeric value!")
    }
  } else {
    parameter$cluster_k = 20
  }
  
  # Check umap_k  ###
  if ("umap_k" %in% names(parameter)) {
    if (converts_to_number(parameter$umap_k)) {
      parameter$umap_k = as.numeric(parameter$umap_k)
    } else {
      error_messages = c(error_messages, "The parameter 'umap_k' is not a numeric value!")
    }
  } else {
    parameter$umap_k = 30
  }
  
  # Check cluster_resolution_test  ###
  if ("cluster_resolution_test" %in% names(parameter) & !is.null(parameter$samples_to_drop)) {
    if (all(purrr::map_lgl(parameter$cluster_resolution_test, converts_to_number))) {
      parameter$cluster_resolution_test = as.numeric(parameter$cluster_resolution_test)
    } else {
      error_messages = c(error_messages, "The parameter 'cluster_resolution_test' does not contain numeric values!")
    }
  } else {
    parameter["cluster_resolution_test"] = list(NULL)
  }
  
  # Check cluster_resolution  ###
  if ("cluster_resolution" %in% names(parameter)) {
    if (converts_to_number(parameter$cluster_resolution)) {
      parameter$cluster_resolution = as.numeric(parameter$cluster_resolution)
    } else {
      error_messages = c(error_messages, "The parameter 'cluster_resolution' does not contain a numeric value!")
    }
  } else {
    parameter$cluster_resolution = 0.5
  }

  # Check marker_padj  ###
  if ("marker_padj" %in% names(parameter)) {
    if (converts_to_number(parameter$marker_padj)) {
      parameter$marker_padj = as.numeric(parameter$marker_padj)
    } else {
      error_messages = c(error_messages, "The parameter 'marker_padj' is not a numeric value!")
    }
  } else {
    parameter$marker_padj = 0.05
  }

  # Check marker_log2FC  ###
  if ("marker_log2FC" %in% names(parameter)) {
    if (converts_to_number(parameter$marker_log2FC)) {
      parameter$marker_log2FC = as.numeric(parameter$marker_log2FC)
    } else {
      error_messages = c(error_messages, "The parameter 'log2fc' is not a numeric value!")
    }
  }
  
  # Check deg_contrasts  ###
  if ("deg_contrasts" %in% names(parameter) & !is.null(parameter$deg_contrasts)) {
    if (is.character(parameter$deg_contrasts) && file.exists(parameter$deg_contrasts)) {
      deg_contrasts_table = openxlsx::read.xlsx(parameter$deg_contrasts)
    } else {
      deg_contrasts_table = parameter$deg_contrasts
    }
    
    if (!is.data.frame(deg_contrasts_table)) {
      error_messages = c(error_messages, "The parameter 'deg_contrasts' must be either a data.frame or an existing Excel file with a valid table in the first sheet!")
    } else {
      if (!all(c("condition_column","condition_group1", "condition_group2") %in% colnames(deg_contrasts_table))){
        error_messages = c(error_messages, "The table specified by parameter 'deg_contrasts' must contain the following columns: 'condition_column', 'condition_group1' and 'condition_group2'!")
      }
      
    }
  }
  
  # Check Enrichr site ("Enrichr", "FlyEnrichr", "WormEnrichr", "YeastEnrichr", "FishEnrichr")  ###
  if ("enrichr_site" %in% names(parameter)) {
    if (!parameter$enrichr_site %in% c("Enrichr", "FlyEnrichr", "WormEnrichr", "YeastEnrichr", "FishEnrichr")) {
      error_messages = c(error_messages, "The parameter 'norm' (normalisation method to use) must be one of: 'Enrichr', 'FlyEnrichr', 'WormEnrichr', 'YeastEnrichr', 'FishEnrichr'!")
    }
  } else {
    parameter$norm = "Enrichr"
  }
  
  # Check enrichr_padj  ###
  if ("enrichr_padj" %in% names(parameter)) {
    if (converts_to_number(parameter$enrichr_padj)) {
      parameter$enrichr_padj = as.numeric(parameter$enrichr_padj)
    } else {
      error_messages = c(error_messages, "The parameter 'enrichr_padj' is not a numeric value!")
    }
  } else {
    parameter$enrichr_padj = 0.05
  }
  
  # Check Enrichr databases of interest  ###
  if ("enrichr_dbs" %in% names(parameter)) {
    parameter$enrichr_dbs = as.character(parameter$enrichr_dbs)
  } else {
    parameter$enrichr_dbs = c("GO_Molecular_Function_2018", "GO_Biological_Process_2018", "GO_Cellular_Component_2018", "Azimuth_Cell_Types_2021", "CellMarker_Augmented_2021", "Descartes_Cell_Types_and_Tissue_2021")
  }

  # Check col  ###
  if ("col" %in% names(parameter)) {
    if (!parameter$col %in% colors() & any(!grepl("#[A-F0-9]+", parameter$col, ignore.case=TRUE))) {
      error_messages = c(error_messages, "The parameter 'col' is missing or not a valid colour! Specify either a valid R color name or a color in hex color code.")
    }
  } else {
    parameter$col = "palevioletred"
  }

  # Check col_palette_samples  ###
  if (!"col_palette_samples" %in% names(parameter)) {
    error_messages = c(error_messages, "The parameter 'col_palette_samples' is missing!")
  }
  
  # Check col_palette_clusters  ###
  if (!"col_palette_clusters" %in% names(parameter)) {
    error_messages = c(error_messages, "The parameter 'col_palette_clusters' is missing!")
  }
  
  # Check debugging_mode   ###
  if ("debugging_mode" %in% names(parameter)) {
    if (!parameter$debugging_mode %in% c("default_debugging", "terminal_debugger", "print_traceback")) {
      error_messages = c(error_messages, "The parameter 'debugging_mode' must be one of: 'default_debugging', 'terminal_debugger', 'print_traceback'!")
    }
  } else {
    parameter$debugging_mode = "default_debugging"
  }
  
  # Add to param object and return
  parameter["error_messages"] = list(NULL)
  if (length(error_messages) > 0) parameter[["error_messages"]] = error_messages
  return(parameter)
}

# Checks if the parameters of the scrnaseq hto workflow are valid and converts them if needed.
#'
#' @param parameter The parameter list used by the scrnaseq_hto script.
#' @return Returns a list with error messages.
check_parameters_scrnaseq_hto = function(parameter) {
  error_messages = c()
  
  # Check author ###
  if (!"author" %in% names(parameter)) {
    error_messages = c(error_messages, "The parameter 'author' is missing!")
  } else {
    parameter$author = as.character(parameter$author)
  }
  
  # Check project_id ###
  if (!"project_id" %in% names(parameter)) {
    error_messages = c(error_messages, "The parameter 'project_id' is missing!")
  } else {
    parameter$project_id = as.character(parameter$project_id)
  }
  
  # Check path_data ###
  if (!"path_data" %in% names(parameter) | !is.character(parameter$path_data) | !dir.exists(parameter$path_data)) {
    error_messages = c(error_messages, "The parameter 'path_data' is missing or not a valid path for 10x counts matrix directory!")
  } else {
    parameter$path_data = file.path(parameter$path_data)
  }
  
  # Check stats ###
  if ("stats" %in% names(parameter)) {
    if (!is.character(parameter$stats) | !file.exists(parameter$stats)) {
      error_messages = c(error_messages, "The parameter 'stats' is not a valid path for a 10x metrics summary file!")
    } else {
      parameter$stats = file.path(parameter$stats)
    }
  }
  
  # Check path_out ###
  if (!"path_out" %in% names(parameter)) {
    error_messages = c(error_messages, "The parameter 'path_out' (path for output) is missing!")
  } else {
    parameter$path_out = file.path(parameter$path_out)
  }
  
  
  # Check downsample_cells_n ###
  if (("downsample_cells_n" %in% names(parameter)) & (!is.null(parameter$downsample_cells_n))) {
    if (!converts_to_number(parameter$downsample_cells_n)) {
      error_messages = c(error_messages, "The parameter 'downsample_cells_n' is not a number!")
    } else {
      parameter$downsample_cells_n = as.numeric(parameter$downsample_cells_n)
    }
  } else {
    parameter["downsample_cells_n"] = list(NULL)
  }
  
  # Check overwrite ###
  if ("overwrite" %in% names(parameter)) {
    if (converts_to_logical(parameter$overwrite)) {
      parameter$overwrite = as.logical(parameter$overwrite)
    } else {
      error_messages = c(error_messages, "The parameter 'overwrite' is not a logical value!")
    }
  } else {
    parameter$overwrite = FALSE
  }
  
  # Check norm (normalisation)  ###
  if ("norm" %in% names(parameter)) {
    if (!parameter$norm %in% c("CLR", "LogNormalize")) {
      error_messages = c(error_messages, "The parameter 'norm' (normalisation method to use) must be one of: 'CLR', 'LogNormalize'!")
    }
  } else {
    parameter$norm = "CLR"
  }
  
  # Check col  ###
  if ("col" %in% names(parameter)) {
    if (!parameter$col %in% colors() & any(!grepl("#[A-F0-9]+", parameter$col, ignore.case=TRUE))) {
      error_messages = c(error_messages, "The parameter 'col' is missing or not a valid colour! Specify either a valid R color name or a color in hex color code.")
    }
  } else {
    parameter$col = "palevioletred"
  }
  
  # Check debugging_mode   ###
  if ("debugging_mode" %in% names(parameter)) {
    if (!parameter$debugging_mode %in% c("default_debugging", "terminal_debugger", "print_traceback")) {
      error_messages = c(error_messages, "The parameter 'debugging_mode' must be one of: 'default_debugging', 'terminal_debugger', 'print_traceback'!")
    }
  } else {
    parameter$debugging_mode = "default_debugging"
  }
  
  # Add to param object and return
  parameter["error_messages"] = list(NULL)
  if (length(error_messages) > 0) parameter[["error_messages"]] = error_messages
  return(parameter)
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
#' @param site The enrichR site (sever) to use.
#' @return Returns a list with error messages.
check_enrichr = function(databases, site="Enrichr") {
  if(is.null(databases) || length(databases)==0) return(c())
  
  # Is enrichR live
  if (is.null(options("enrichR.live")) || !options("enrichR.live")[[1]]) {
    return("EnrichR is not available or cannot connect to the databases")
  }
  
  # Set Enrichr site
  suppressMessages(enrichR::setEnrichrSite(site))

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
                        "knitcitations", "sceasy", "ROpenSci")
  
  is_installed = packages_installed(packages=required_packages)
  if(any(!is_installed)) {
    return(paste0("The R package '", required_packages[!is_installed],"' is not installed!"))
  } else {
    return(c())
  }
}

#' Checks if all packages required for the scrnaseq hto workflow are installed.
#'
#' @return Returns a list with error messages.
check_installed_packages_scrnaseq_hto = function() {
  required_packages = c("Seurat", "ggplot2", "patchwork", "magrittr",
                        "reticulate", "future", "knitr",
                        "dplyr", "tidyr", "purrr", "stringr", 
                        "Matrix", "kableExtra", "DT", "ggsci",
                        "openxlsx", "readr", "R.utils",
                        "sessioninfo", "knitcitations", "sceasy")
  
  is_installed = packages_installed(packages=required_packages)
  if(any(!is_installed)) {
    return(paste0("The R package '", required_packages[!is_installed],"' is not installed!"))
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

