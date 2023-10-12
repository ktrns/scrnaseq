# Dataset provided by 10x: 5k Peripheral blood mononuclear cells (PBMCs) from a healthy donor with cell surface proteins (v3 chemistry)
# https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.2/5k_pbmc_protein_v3

# Important: Run this script in its directory

unlink("filtered_feature_bc_matrix.h5", recursive=T)

# download and untar
url = 'http://cf.10xgenomics.com/samples/cell-exp/3.0.2/5k_pbmc_protein_v3/5k_pbmc_protein_v3_filtered_feature_bc_matrix.h5'
curl::curl_download(url=url, destfile="filtered_feature_bc_matrix.h5")

# get metrics summary
url = 'http://cf.10xgenomics.com/samples/cell-exp/3.0.2/5k_pbmc_protein_v3/5k_pbmc_protein_v3_metrics_summary.csv'
curl::curl_download(url=url, destfile='metrics_summary.csv')

# get some analysis results and use "analysis/clustering/graphclust/clusters.csv" as barcode metadata
url = 'https://cf.10xgenomics.com/samples/cell-exp/3.1.0/5k_pbmc_protein_v3/5k_pbmc_protein_v3_analysis.tar.gz'
curl::curl_download(url=url, destfile=basename(path=url))
untar(tarfile = basename(path=url))
unlink(basename(path=url))

