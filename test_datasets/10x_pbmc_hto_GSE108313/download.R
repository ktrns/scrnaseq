# Downloaded 10X data set, 'pbmc_hto_mtx.rds' and 'pbmc_umi_mtx.rds'
# https://www.dropbox.com/sh/ntc33ium7cg1za1/AAD_8XIDmu4F7lJ-5sp-rGFYa?dl=0
#
# Important: Run this script in its directory

# Clear old files
unlink("demultiplexed", recursive=TRUE)
unlink("counts", recursive=TRUE)
dir.create("counts", showWarnings=FALSE)
unlink("download", recursive=TRUE)
dir.create("download", showWarnings=FALSE)

# Download
url = 'https://www.dropbox.com/sh/ntc33ium7cg1za1/AAAp2nqrR0c3WOH6tjegOVm9a/pbmc_hto_mtx.rds?dl=1'
curl::curl_download(url=url, destfile="download/pbmc_hto_mtx.rds")

url = 'https://www.dropbox.com/sh/ntc33ium7cg1za1/AACJ_JzX7ZYblSn1-HANBM_7a/pbmc_umi_mtx.rds?dl=1'
curl::curl_download(url=url, destfile="download/pbmc_umi_mtx.rds")

# Read downloaded files
pbmc_umi = readRDS("download/pbmc_umi_mtx.rds")
pbmc_hto = readRDS("download/pbmc_hto_mtx.rds")

# Rename HTOs
rownames(pbmc_hto) = gsub(rownames(pbmc_hto), pattern="HTO", replacement="hto", fixed=TRUE)
rownames(pbmc_hto) = gsub(rownames(pbmc_hto), pattern="_", replacement="", fixed=TRUE)
hto_names = rownames(pbmc_hto)

# Map to Ensembl
annot_mart = biomaRt::useEnsembl("ensembl", dataset="hsapiens_gene_ensembl", mirror="www")
annot_ensembl = biomaRt::getBM(mart=annot_mart, attributes=c("ensembl_gene_id", "external_gene_name"), useCache=FALSE)
idx = match(rownames(pbmc_umi), annot_ensembl$external_gene_name) # if multiple symbols match, take the first one
genes_ensembl = setNames(annot_ensembl$ensembl_gene_id[idx], rownames(pbmc_umi))

# Remove rows in the counts table with no Ensembl ID
message("Number of gene symbols not found through biomaRt: ", sum(is.na(idx)))
message("Number of genes with no Ensembl: ", sum(is.na(genes_ensembl)))
pbmc_umi = pbmc_umi[!is.na(genes_ensembl), ]

# Merge count matrices
joint_barcodes = intersect(colnames(pbmc_umi), colnames(pbmc_hto))
pbmc_umi = pbmc_umi[, joint_barcodes]
pbmc_hto = pbmc_hto[, joint_barcodes]
pbmc_counts = rbind(pbmc_umi, pbmc_hto)

feat = data.frame(V1=genes_ensembl[rownames(pbmc_umi)], V2=rownames(pbmc_umi), V3="Gene Expression")
feat = rbind(feat, data.frame(V1=rownames(pbmc_hto), V2=rownames(pbmc_hto), V3="Antibody Capture"))

# Write matrix.mtx.gz, barcodes.tsv.gz and features.tsv.gz
mh = file.path("counts", "matrix.mtx")
Matrix::writeMM(pbmc_counts, file=mh)
R.utils::gzip(mh, overwrite=TRUE)

bh = gzfile(file.path("counts", "barcodes.tsv.gz"), open="wb")
write(colnames(pbmc_counts), file=bh)
close(bh)

fh = gzfile(file.path("counts", "features.tsv.gz"), open="wb")
write.table(feat, file=fh, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
close(fh)

unlink("download", recursive=TRUE)
