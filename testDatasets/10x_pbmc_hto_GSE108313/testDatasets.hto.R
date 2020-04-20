# Downloaded 10X data set, 'pbmc_hto_mtx.rds' and 'pbmc_umi_mtx.rds'
# https://www.dropbox.com/sh/ntc33ium7cg1za1/AAD_8XIDmu4F7lJ-5sp-rGFYa?dl=0

# Read downloaded files
library(Seurat)
pbmc.umi = readRDS("pbmc_umi_mtx.rds")
pbmc.hto = readRDS("pbmc_hto_mtx.rds")

rownames(pbmc.hto) = gsub(rownames(pbmc.hto), pattern="HTO", replacement="hto", fixed=TRUE)
rownames(pbmc.hto) = gsub(rownames(pbmc.hto), pattern="_", replacement="", fixed=TRUE)

# Merge count matrices
joint.barcodes=intersect(colnames(pbmc.umi), colnames(pbmc.hto))
pbmc.umi = pbmc.umi[,joint.barcodes]
pbmc.hto = pbmc.hto[,joint.barcodes]
pbmc.counts = rbind(pbmc.umi, pbmc.hto)

# Write to 3 CellRanger output files
library(DropletUtils)
write10xCounts(path="./cellranger", x=pbmc.counts)

# Add HTO to genes.tsv
genes = read.delim("./cellranger/genes.tsv", header=FALSE)
genes = cbind(genes, "Gene Expression")
colnames(genes) = c("Ensembl", "GeneSymbol", "Type")
genes$Type = as.character(genes$Type)
genes$Type[grep("^hto", genes$GeneSymbol, perl=TRUE)] = "Antibody Capture"

# Add Ensembl IDs to features.tsv
annot.ensembl = read.delim("~/ensembl/hsapiens_gene_ensembl.98.csv")
idx = match(genes$GeneSymbol, annot.ensembl$Gene_Symbol)
genes$Ensembl = annot.ensembl$Ensembl_ID[idx]

# Write genes to features.tsv
write.table(genes, file="cellranger/features.tsv", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

# Gzip all three files
system("gzip cellranger/barcodes.tsv")
system("gzip cellranger/features.tsv")
system("gzip cellranger/matrix.mtx")

# Remove genes.tsv
# This is necessary, so Seurat recognises the CellRanger version
system("rm cellranger/genes.tsv")

# For using this test dataset, you need to provide the path to the three generated files
#   barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz
