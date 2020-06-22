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

# Report session info in a table
scrnaseq_session_info = function(path_to_git) {
  
  out=matrix(NA, nrow=0, ncol=2)
  colnames(out) = c("Name", "Version")
  
  repo = system(paste0("git --git-dir ", path_to_git, "/.git log --format='%H' -n 1"), intern=TRUE)
  out = rbind(out, c("ktrns/scrnaseq", repo))
  
  info_session = sessionInfo()
  out = rbind(out, c("R", info_session$R.version$version.string))
  out = rbind(out, c("Platform", info_session$platform))
  out = rbind(out, c("Operating system", info_session$running))
  
  info_pkgs = sessioninfo::package_info()
  out = rbind(out, c("Packages", paste(paste(info_pkgs$package, info_pkgs$loadedversion, sep=""), collapse=", ")))
  
  return(out)
}
