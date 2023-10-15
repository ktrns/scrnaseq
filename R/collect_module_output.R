#!/usr/bin/env Rscript

#QUARTO_PROJECT_RENDER_ALL	Set to “1” if this is a render of all files in the project (as opposed to an incremental render or a render for preview). This unset if Quarto is not rendering all files.
#QUARTO_PROJECT_OUTPUT_DIR	Output directory
#QUARTO_PROJECT_INPUT_FILES	Newline separated list of all input files being rendered (passed only to pre-render)
#QUARTO_PROJECT_OUTPUT_FILES	Newline separated list of all output files rendered (passed only to post-render).

# Libraries
library(magrittr)

# Splits a path
split_path = function(x) {
  if (dirname(x)==x) return(x) else return(c(split_path(dirname(x)), basename(x)))
}

# The environment variable QUARTO_PROJECT_OUTPUT_FILES should contains a newline-separated list of all quarto output files rendered 
quarto_output_files = Sys.getenv("QUARTO_PROJECT_OUTPUT_FILES")
quarto_project_output_dir = Sys.getenv("QUARTO_PROJECT_OUTPUT_DIR")

#quarto_output_files = "_book\\index.html\n_book\\modules\\read_data\\read_data.html\n_book\\modules\\visualisation\\visualisation.html"


if (nchar(quarto_output_files) > 0 & nchar(quarto_project_output_dir) > 0) {
  quarto_output_files = stringi::stri_split_lines(quarto_output_files, omit_empty=TRUE) %>% unlist()

  # Find all modules that were rendered
  modules_rendered = purrr::map(quarto_output_files, function(f) {
    m = NULL
    
    f = split_path(f)
    i = which(f == "modules")
    if (length(i) > 0) {
      i = i[1]
      if (length(f) > i) {
        d = do.call(file.path, as.list(f[1:(i+1)]))
        if (dir.exists(d)) {
          m = f[i+1]
        }
      }
    }
    
    return(m)
  }) %>% purrr::flatten_chr() %>% unique()
  
  print(paste("Module", modules_rendered, collapse = "\n"))

}