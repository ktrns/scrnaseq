# Dataset provided by 10x: 5k Peripheral blood mononuclear cells (PBMCs), Chromatin accessibility data.
# https://www.10xgenomics.com/resources/datasets/5-k-peripheral-blood-mononuclear-cells-pbm-cs-from-a-healthy-donor-next-gem-v-1-1-1-1-standard-2-0-0

# Important: Run this script in its directory

# download and untar
url = 'https://cf.10xgenomics.com/samples/cell-atac/2.0.0/atac_pbmc_5k_nextgem/atac_pbmc_5k_nextgem_filtered_peak_bc_matrix.h5'
curl::curl_download(url=url, destfile="filtered_peak_bc_matrix.h5")

url = 'https://cf.10xgenomics.com/samples/cell-atac/2.0.0/atac_pbmc_5k_nextgem/atac_pbmc_5k_nextgem_filtered_tf_bc_matrix.h5'
curl::curl_download(url=url, destfile="filtered_tf_bc_matrix.h5")

# get metadata
url = 'https://cf.10xgenomics.com/samples/cell-atac/2.0.0/atac_pbmc_5k_nextgem/atac_pbmc_5k_nextgem_singlecell.csv'
curl::curl_download(url=url, destfile="single_cell.csv")

# get fragments
url = 'https://cf.10xgenomics.com/samples/cell-atac/2.0.0/atac_pbmc_5k_nextgem/atac_pbmc_5k_nextgem_fragments.tsv.gz'
curl::curl_download(url=url, destfile="fragments.tsv.gz")
url = 'https://cf.10xgenomics.com/samples/cell-atac/2.0.0/atac_pbmc_5k_nextgem/atac_pbmc_5k_nextgem_fragments.tsv.gz.tbi'
curl::curl_download(url=url, destfile="fragments.tsv.gz.tbi")

# get metrics summary
url = 'https://cf.10xgenomics.com/samples/cell-atac/2.0.0/atac_pbmc_5k_nextgem/atac_pbmc_5k_nextgem_summary.csv'
curl::curl_download(url=url, destfile='metrics_summary.csv')
