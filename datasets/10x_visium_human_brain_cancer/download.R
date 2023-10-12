# Dataset provided by 10x: Visium spatial gene expression dataset, human brain cancer, 11 mm capture area (FFPE)
# https://www.10xgenomics.com/resources/datasets/human-brain-cancer-11-mm-capture-area-ffpe-2-standard

# Important: Run this script in its directory

unlink("filtered_feature_bc_matrix.h5", recursive=T)

# download and untar
url = 'https://cf.10xgenomics.com/samples/spatial-exp/2.0.1/CytAssist_11mm_FFPE_Human_Glioblastoma/CytAssist_11mm_FFPE_Human_Glioblastoma_filtered_feature_bc_matrix.h5'
curl::curl_download(url=url, destfile="filtered_feature_bc_matrix.h5")

# get metrics summary
url = 'https://cf.10xgenomics.com/samples/spatial-exp/2.0.1/CytAssist_11mm_FFPE_Human_Glioblastoma/CytAssist_11mm_FFPE_Human_Glioblastoma_metrics_summary.csv'
curl::curl_download(url=url, destfile='metrics_summary.csv')

# get spatial data
url = 'https://cf.10xgenomics.com/samples/spatial-exp/2.0.1/CytAssist_11mm_FFPE_Human_Glioblastoma/CytAssist_11mm_FFPE_Human_Glioblastoma_spatial.tar.gz'
curl::curl_download(url=url, destfile="spatial.tar.gz")
untar(tarfile = basename(path="spatial.tar.gz"))
unlink(basename(path="spatial.tar.gz"))

