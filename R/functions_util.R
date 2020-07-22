# 'Calculate enrichment of cells per sample per cluster
#' 
#' @param sc Seurat object
cells_fisher = function(sc) {
  cell_samples = sc[[]] %>% dplyr::pull(orig.ident) %>% unique() %>% sort()
  cell_clusters = sc[[]] %>% dplyr::pull(seurat_clusters) %>% unique() %>% sort()
  out = matrix(0+NA, nrow=length(cell_clusters), ncol=0) %>% as.data.frame()
  for(s in cell_samples) {
    ft.list = lapply(cell_clusters, function(cl) { 
      a = sc[[]] %>% dplyr::filter(orig.ident==s, seurat_clusters==cl) %>% dplyr::count() %>% as.numeric()
      b = sc[[]] %>% dplyr::filter(orig.ident!=s, seurat_clusters==cl) %>% dplyr::count() %>% as.numeric()
      c = sc[[]] %>% dplyr::filter(orig.ident==s, seurat_clusters!=cl) %>% dplyr::count() %>% as.numeric()
      d = sc[[]] %>% dplyr::filter(orig.ident!=s, seurat_clusters!=cl) %>% dplyr::count() %>% as.numeric()
      tbl.2by2 = matrix(c(a, b, c, d), ncol=2, nrow=2, byrow=TRUE)
      ft = fisher.test(tbl.2by2, alternative="greater")
      return(c(oddsRatio=round(as.numeric(ft$estimate), 2),
               p=round(as.numeric(ft$p.value), 2)))
    })
    ft.matrix = purrr::reduce(ft.list, .f=rbind)
    colnames(ft.matrix) = paste0(s, ".", colnames(ft.matrix))
    out = cbind(out, ft.matrix)
  }
  return(out)
}

#' Given a vector, report at most n elements as concatenated string.
#' 
#' @param x A vector.
#' @param n Number of elements to report at most.
#' @param sep Separator for string concatenation.
#' @return A string with at most n elements to concatenated.
first_n_elements_to_string = function(x, n=5, sep=",") {
  s = paste(x[1:min(n,length(x))], collapse=sep)
  if (length(x)>n) s = paste(s, "...", sep=sep)
  return(s)
}

#' Report session info in a table
#' 
#' @return The session info as table.
scrnaseq_session_info = function(path_to_git=".") {
  out=matrix(NA, nrow=0, ncol=2)
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

#' Returns the names of an object
#' @param x A list or vector with names.
#' @return A named vector with names as names and names as values.
list_names = function(x) {
  return(setNames(names(x), names(x)))
}

#' Returns a vector with its values as names.
#' @param x A vector.
#' @return A vector with its values as names.
values_to_names = function(x) {
  return(setNames(x,x))
}

#' Returns the indices of an object
#' @param x A list or vector with names.
#' @return A named vector with names as names and indices as values.
list_indices = function(x) {
  return(setNames(seq(x), names(x)))
}

#' Wrapper around the biomaRt::useEnsembl function to cope with unavailable Ensembl mirrors. Tries different Ensembl mirrors and returns a mart object with the mirror that works.
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
  
  mart_obj = NA
  if(is.null(host)) {
    # Test and if a mirror is not available, check the next one
    for(m in mirrors_to_test) {
      mart_obj = tryCatch({
        biomaRt::useEnsembl(biomart=biomart, dataset=dataset, mirror=m, version=version)
      },
      error=function(cond) {
        return(NA)
      })
    
      if(!is.na(mart_obj)) break
    }
    if(is.na(mart_obj)) stop("The requested Ensembl mirror(s) are not available.")
    
  } else {
    # Use specific host
    mart_obj = tryCatch({
      biomaRt::useEnsembl(biomart=biomart, dataset=dataset, host=host, version=version)
    },
    error=function(cond) {
      return(NA)
    })
    if(is.na(mart_obj)) stop("The requested Ensembl host is not available.")
  }
  
  return(mart_obj)
}


#' Returns the mirror of a biomaRt object.
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
#' @param num_colours The number of colours to generate.
#' @param palette A palette function for generating the colours.
#' @param palette_options List of additional arguments (beside alpha) to pass on to the palette function.
#' @param alphas Alpha value(s) to use. If the number of colours exceeds the palette, multiple alpha value can be provided to generate more colours.
#' @return The generated colours.
GenerateColours = function(num_colours, palette=ggsci::pal_igv, alphas=c(1,0.7,0.3), palette_options=list()) {
  colours = purrr::flatten_chr(purrr::map(alphas, function(a) {
    palette_options[["alpha"]] = a
    cols = suppressWarnings(do.call(do.call(ggsci::pal_d3,palette_options),list(100)))
    cols[!is.na(cols)]
  }))
  
  if (num_colours>length(colours)) {
    stop("GenerateColours: Cannot generate the requested number of colours. Please change palette or add alpha values.")
  }
  
  return(colours[1:num_colours])
}

#' Prints a message formatted for markdown. 
#' See: https://www.w3schools.com/bootstrap/bootstrap_alerts.asp and https://bookdown.org/yihui/rmarkdown-cookbook/output-hooks.html
#' @param x The message.
#' @param options Further options.
#' @return The message formatted for markdown.
Message = function(x, options){
  x = gsub('^##','',x)
  msg = paste(c('\n\n:::{class="alert alert-info alert-dismissible"}',
          '<a href="#" class="close" data-dismiss="alert" aria-label="close">&times;</a>',
          '<strong>Information:</strong>',
          x,
          ':::\n'), collapse = '\n')
  return(msg)
}

#' Prints a warning formatted for markdown. 
#' See: https://www.w3schools.com/bootstrap/bootstrap_alerts.asp and https://bookdown.org/yihui/rmarkdown-cookbook/output-hooks.html
#' @param x The message.
#' @param options Further options.
#' @return The message formatted for markdown.
Warning = function(x, options){
  x = gsub('^##','',x)
  warn = paste(c('\n\n:::{class="alert alert-warning alert-dismissible"}',
                '<a href="#" class="close" data-dismiss="alert" aria-label="close">&times;</a>',
                '<strong>Warning:</strong>',
                x,
                ':::\n'), collapse = '\n')
  return(warn)
}
