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
scrnaseq_session_info = function() {
  
  out=matrix(NA, nrow=0, ncol=2)
  colnames(out) = c("Name", "Version")
  
  repo = system("git log --format='%H' -n 1", intern=TRUE)
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

#' Wrapper around the biomaRt::useEnsembl function to cope with unavailable Ensembl mirrors. Tries different Ensembl mirrors and returns a mart object with the mirror that works.
#' @param biomart A biomaRt database name. Possible database names can be retrieved with the function listEnsembl().
#' @param dataset Dataset you want to use. Possible dataset names can be retrieved with the function listDatasets(mart_obj).
#' @param mirror Specify an Ensembl mirror to connect to. The valid options here are 'www', 'uswest', 'useast', 'asia'. If no mirror is specified, then the first mirror that works will be used. Will be ignored if a host is specified.
#' @param version Ensembl version to connect to when wanting to connect to an archived Ensembl version.
#' @return A biomaRt object.
#' @param host Host to connect to. Only needs to be specified if different from www.ensembl.org. 
GetBiomart = function(biomart, dataset, mirror=NULL, version=NULL, host=NULL) {
  
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
