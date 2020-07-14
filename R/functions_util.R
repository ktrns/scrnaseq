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

# Report parameters of hto report in a table
scrnaseq_hto_parameter_info = function(project, path_data, path_out, hto_names, mt, color, sample_cells) {
  out=matrix(NA, nrow=0, ncol=2)
  colnames(out) = c("Parameter", "Value")
  
  out = rbind(out, c("Project ID", project))
  out = rbind(out, c("Input data path in case Cell Ranger was run", path_data))
  out = rbind(out, c("Output path", path_out))
  out = rbind(out, c("HTO names", toString(hto_names)))
  out = rbind(out, c("Prefix of mitochondrial genes", mt))
  out = rbind(out, c("Main color to use for plots", toString(color)))
  out = rbind(out, c("Sample data", toString(sample_cells)))
  
  return(out)
}