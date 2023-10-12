# Minimal dataset for testing.
#
# The 10x dataset "1k Peripheral blood mononuclear cells (PBMCs)" is used to create two artifical samples with 100 cells and 500 genes each.
# Dataset: https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.0/pbmc_1k_v3

# Important: Run this script in its directory

unlink("sample1", recursive=T)
dir.create("sample1", showWarnings=FALSE)
dir.create("sample1/filtered_feature_bc_matrix", showWarnings=FALSE)

unlink("sample2", recursive=T)
dir.create("sample2", showWarnings=FALSE)
dir.create("sample2/filtered_feature_bc_matrix", showWarnings=FALSE)

# Download dataset
url = 'https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_1k_v3/pbmc_1k_v3_filtered_feature_bc_matrix.tar.gz'
curl::curl_download(url=url, destfile=basename(path=url))
untar(tarfile = basename(path=url))
unlink(basename(path=url))

# Read 10x data
tenx_data = Seurat::Read10X("filtered_feature_bc_matrix")

# Create two samples with 300 cells
set.seed(11)
cells = sample(colnames(tenx_data), 600)
cells_sample1 = cells[1:100]
cells_sample2 = cells[101:200]

sample1_data = tenx_data[, cells_sample1]
sample2_data = tenx_data[, cells_sample2]

# Write to disk
mh = file.path("sample1", "filtered_feature_bc_matrix", "matrix.mtx")
Matrix::writeMM(sample1_data, file=mh)
R.utils::gzip(mh, overwrite=TRUE)

bh = gzfile(file.path("sample1", "filtered_feature_bc_matrix", "barcodes.tsv.gz"), open="wb")
write(colnames(sample1_data), file=bh)
close(bh)

fh = gzfile(file.path("sample1", "filtered_feature_bc_matrix", "features.tsv.gz"), open="wb")
write.table(rownames(sample1_data), file=fh, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
close(fh)

mh = file.path("sample2", "filtered_feature_bc_matrix", "matrix.mtx")
Matrix::writeMM(sample2_data, file=mh)
R.utils::gzip(mh, overwrite=TRUE)

bh = gzfile(file.path("sample2", "filtered_feature_bc_matrix", "barcodes.tsv.gz"), open="wb")
write(colnames(sample2_data), file=bh)
close(bh)

fh = gzfile(file.path("sample2", "filtered_feature_bc_matrix", "features.tsv.gz"), open="wb")
write.table(rownames(sample2_data), file=fh, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
close(fh)

unlink("filtered_feature_bc_matrix", recursive=T)
