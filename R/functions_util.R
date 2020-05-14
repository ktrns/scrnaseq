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