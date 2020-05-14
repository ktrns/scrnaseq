# Downloaded 10X data set, 'pbmc_hto_mtx.rds' and 'pbmc_umi_mtx.rds'
# https://www.dropbox.com/sh/ntc33ium7cg1za1/AAD_8XIDmu4F7lJ-5sp-rGFYa?dl=0

# Set working directory accordingly
setwd("~/scrnaseq/test_datasets/10x_pbmc_hto_GSE108313")
unlink("counts", recursive=T)

# Read downloaded files
pbmc_umi = readRDS("pbmc_umi_mtx.rds")
pbmc_hto = readRDS("pbmc_hto_mtx.rds")

# Rename HTOs
rownames(pbmc_hto) = gsub(rownames(pbmc_hto), pattern="HTO", replacement="hto", fixed=TRUE)
rownames(pbmc_hto) = gsub(rownames(pbmc_hto), pattern="_", replacement="", fixed=TRUE)

# Merge count matrices
joint_barcodes = intersect(colnames(pbmc_umi), colnames(pbmc_hto))
pbmc_umi = pbmc_umi[,joint_barcodes]
pbmc_hto = pbmc_hto[,joint_barcodes]
pbmc_counts = rbind(pbmc_umi, pbmc_hto)

# Map to Ensembl
annot_mart = biomaRt::useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", mirror="www")
annot_ensembl = biomaRt::getBM(mart = annot_mart, attributes = c("ensembl_gene_id", "external_gene_name"))
idx = match(rownames(pbmc_counts), annot_ensembl$external_gene_name) # if multiple symbols match, take the first one
genes_ensembl = setNames(annot_ensembl$ensembl_gene_id[idx], rownames(pbmc_counts))

# Remove rows in the counts table with no Ensembl ID
message("Number of gene symbols not found through biomaRt: ", sum(is.na(idx)))
message("Number of genes with no Ensembl: ", sum(is.na(genes_ensembl)))
pbmc_counts = pbmc_counts[!is.na(genes_ensembl),]

# Write to 3 CellRanger output files
DropletUtils::write10xCounts(path="./counts", x=pbmc_counts)

# Add HTO to genes.tsv
genes = read.delim("./counts/genes.tsv", header=FALSE)
genes = cbind(genes, "Gene Expression")
colnames(genes) = c("Ensembl", "GeneSymbol", "Type")
genes$Type = as.character(genes$Type)
genes$Type[grep("^hto", genes$GeneSymbol, perl=TRUE)] = "Antibody Capture"

# Add Ensembl IDs to features.tsv
genes$Ensembl = genes_ensembl[genes$GeneSymbol]

# Write genes to features.tsv
write.table(genes, file="counts/features.tsv", sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE)

# Gzip all three files
system("gzip counts/barcodes.tsv")
system("gzip counts/features.tsv")
system("gzip counts/matrix.mtx")

# Remove genes.tsv
# This is necessary, so Seurat recognises the CellRanger version
system("rm counts/genes.tsv")
system("rm pbmc_umi_mtx.rds pbmc_hto_mtx.rds")

# For using this test dataset, you need to provide the path to the three generated files
#   barcodes.tsv.gz, features.tsv.gz, matrix.mtx.gz
