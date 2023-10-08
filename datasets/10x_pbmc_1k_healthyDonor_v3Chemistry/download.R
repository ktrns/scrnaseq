# Dataset provided by 10x: 1k Peripheral blood mononuclear cells (PBMCs) from a healthy donor
# https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.0/pbmc_1k_v3

# Important: Run this script in its directory

unlink("filtered_feature_bc_matrix", recursive=T)

# download and untar
url = 'https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_filtered_feature_bc_matrix.tar.gz'
curl::curl_download(url=url, destfile=basename(path=url))
untar(tarfile = basename(path=url))
unlink(basename(path=url))

# get metrics summary
url = 'https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_metrics_summary.csv'
curl::curl_download(url=url, destfile='metrics_summary.csv')

# get some analysis results and use "analysis/clustering/graphclust/clusters.csv" as barcode metadata
url = 'https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_analysis.tar.gz'
curl::curl_download(url=url, destfile=basename(path=url))
untar(tarfile = basename(path=url))
unlink(basename(path=url))
