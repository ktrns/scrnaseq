# Do not convert character vectors automatically into factors
options(stringsAsFactors=FALSE)

# Get rid of info of dplyr when grouping: `summarise()` regrouping output by 'species' (override with `.groups` argument)
options(dplyr.summarise.inform=FALSE)

# The default number of cores
options(mc.cores=8)

# The total size of all global objects that need to be exported for the future expression for parallel computations
options(future.globals.maxSize=+Inf)

# The strategy used for parallel computations with future
options(future.plan="multisession")

# Use v5 assays in Seurat
options(Seurat.object.assay.version="v5")

# Python3 path needed for clustering, umap, other python packages
reticulate_python3_path = unname(Sys.which("python3"))
Sys.setenv(RETICULATE_PYTHON=reticulate_python3_path)

# Buffer for reading large text files
Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 2)

# Hook to measure the computation time of a chunk
# Activate with '#| timeit: true'. Deactivate with '#| timeit: null'
knitr::knit_hooks$set(timeit = function(before, options, envir) {
  if(before) {
    ## code to be run before a chunk
    tictoc::tic()
  } else {
    ## code to be run after a chunk
    elapsed = tictoc::toc()$toc
    print(paste0("Execution took ", elapsed, " seconds"))
  }
})

# Set seed
options(random_seed=11)
