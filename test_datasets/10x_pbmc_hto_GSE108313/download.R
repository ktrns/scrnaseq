# PBMC dataset with hash tags from "Cell hashing enable sample multiplexing, multiplet identification and super-loading on droplet-based single cell RNA-sequencing platforms" (Stoeckius et al 2018)

# Downloaded 10X data set, 'pbmc_hto_mtx.rds' and 'pbmc_umi_mtx.rds'
# https://www.dropbox.com/sh/ntc33ium7cg1za1/AAD_8XIDmu4F7lJ-5sp-rGFYa?dl=0

unlink("counts", recursive=T)

# Read downloaded files
pbmc_umi = readRDS("pbmc_umi_mtx.rds")
pbmc_hto = readRDS("pbmc_hto_mtx.rds")

rownames(pbmc_hto) = gsub(x=rownames(pbmc_hto), pattern="HTO", replacement="hto", fixed=TRUE)
rownames(pbmc_hto) = gsub(x=rownames(pbmc_hto), pattern="_", replacement="", fixed=TRUE)

# Merge count matrices
joint_barcodes=intersect(colnames(pbmc_umi), colnames(pbmc_hto))
pbmc_umi = pbmc_umi[,joint_barcodes]
pbmc_hto = pbmc_hto[,joint_barcodes]
pbmc_counts = rbind(pbmc_umi, pbmc_hto)

# get ensembl ids from biomart
annot_mart = biomaRt::useEnsembl("ensembl", dataset="hsapiens_gene_ensembl", mirror="asia")
annot_ensembl = biomaRt::getBM(mart=annot_mart, attributes=c("ensembl_gene_id", "external_gene_name"))
annot_ensembl = annot_ensembl[!duplicated(annot_ensembl$external_gene_name),]

rows_with_ensembl = which(rownames(pbmc_counts) %in% annot_ensembl$external_gene_name)
pbmc_counts = pbmc_counts[rows_with_ensembl,]

# write matrix.mtx.gz, barcodes.tsv.gz and features.tsv.gz
dir.create("counts", showWarnings=FALSE)
mh = file.path("counts","matrix.mtx")
Matrix::writeMM(pbmc_counts,file=mh)
R.utils::gzip(mh,overwrite=T)

bh = gzfile(file.path("counts", "barcodes.tsv.gz"), open="wb")
write(colnames(pbmc_counts), file=bh)
close(bh)


feat = merge(annot_ensembl,data.frame(rowname = rownames(pbmc_counts), row_num = 1:nrow(pbmc_counts),stringsAsFactors=FALSE), by.x=2, by.y=1)
feat = feat[order(feat$row_num),]
fh = gzfile(file.path("counts", "features.tsv.gz"), open="wb")
write.table(feat[,c(2,1)], file=fh, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
close(fh)
