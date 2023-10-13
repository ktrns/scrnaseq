# Dataset provided by 10x: fresh frozen mouse brain for xenium (small subset)
# https://www.10xgenomics.com/resources/datasets/fresh-frozen-mouse-brain-for-xenium-explorer-demo-1-standard

# Important: Run this script in its directory
files = c("analysis_summary.html", "analysis.zarr.zip", "cell_boundaries.csv.gz", "cell_boundaries.parquet", "cell_feature_matrix.h5",
          "cell_feature_matrix.zarr.zip", "cells.csv.gz", "cells.parquet", "cells.zarr.zip", "experiment.xenium", "gene_panel.json",
          "metrics_summary.csv", "morphology_focus.ome.tif", "morphology_mip.ome.tif", "morphology.ome.tif", "nucleus_boundaries.csv.gz",
          "nucleus_boundaries.parquet", "transcripts.csv.gz", "transcripts.parquet", "transcripts.zarr.zip")
unlink(files)
dirs = c("analysis", "cell_feature_matrix") 
unlink(dirs, recursive=TRUE)

# download and unzip
url = 'https://cf.10xgenomics.com/samples/xenium/1.0.2/Xenium_V1_FF_Mouse_Brain_Coronal_Subset_CTX_HP/Xenium_V1_FF_Mouse_Brain_Coronal_Subset_CTX_HP_outs.zip'
curl::curl_download(url=url, destfile=basename(path=url))
utils::unzip(basename(path=url))
unlink(basename(path=url))
