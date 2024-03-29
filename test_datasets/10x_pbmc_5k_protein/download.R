# Dataset provided by 10: 5k Peripheral blood mononuclear cells (PBMCs) from a healthy donor with cell surface proteins (v3 chemistry)
# https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.2/5k_pbmc_protein_v3

# Important: Run this script in its directory

unlink("filtered_feature_bc_matrix", recursive=T)
unlink("counts", recursive=T)

# download and untar
url = 'http://cf.10xgenomics.com/samples/cell-exp/3.0.2/5k_pbmc_protein_v3/5k_pbmc_protein_v3_filtered_feature_bc_matrix.tar.gz'
curl::curl_download(url=url, destfile=basename(path=url))
untar(tarfile = basename(path=url))
unlink(basename(path=url))

file.rename("filtered_feature_bc_matrix","counts")
