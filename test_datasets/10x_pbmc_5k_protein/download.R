# Dataset provided by 10x: 5k Peripheral blood mononuclear cells (PBMCs) from a healthy donor with cell surface proteins (v3 chemistry)
# https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.2/5k_pbmc_protein_v3

# Important: Run this script in its directory

unlink("counts", recursive=T)

# download and untar
url = 'http://cf.10xgenomics.com/samples/cell-exp/3.0.2/5k_pbmc_protein_v3/5k_pbmc_protein_v3_filtered_feature_bc_matrix.tar.gz'
curl::curl_download(url=url, destfile=basename(path=url))
untar(tarfile = basename(path=url))
unlink(basename(path=url))
file.rename("filtered_feature_bc_matrix","counts")

# get metrics summary
url = 'http://cf.10xgenomics.com/samples/cell-exp/3.0.2/5k_pbmc_protein_v3/5k_pbmc_protein_v3_metrics_summary.csv'
curl::curl_download(url=url, destfile='metrics_summary.csv')

# write marker file
openxlsx::write.xlsx(data.frame(bcell=c("ENSG00000105369", "ENSG00000156738"),
                                tcell=c("ENSG00000167286", NA),
                                tcell_cd8p=c("ENSG00000153563", "ENSG00000172116"),
                                nk=c("ENSG00000115523", "ENSG00000105374"),
                                myeloid=c("ENSG00000101439", "ENSG00000090382"),
                                monocytes=c("ENSG00000203747", NA),
                                dendritic=c("ENSG00000179639", NA)),
                     "known_markers.xlsx")