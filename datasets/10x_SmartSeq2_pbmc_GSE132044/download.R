# 10x and Smartseq2 datasets of PBMC cells from "Systematic comparative analysis of single cell RNA-sequencing methods" (Ding et al 2018)

# CAREFULL: BIG DATASET

# make download directory
unlink("download", recursive=TRUE)
dir.create("download", showWarnings=FALSE)
unlink("counts", recursive=TRUE)
dir.create("counts", showWarnings=FALSE)

# download from GEO
counts_url = "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE132nnn/GSE132044/suppl/GSE132044_pbmc_hg38_count_matrix.mtx.gz"
curl::curl_download(url=counts_url, destfile=file.path("download", "matrix.mtx.gz"))

genes_url = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE132nnn/GSE132044/suppl/GSE132044_pbmc_hg38_gene.tsv.gz'
curl::curl_download(url=genes_url, destfile=file.path("download", "features.orig.tsv.gz"))

cells_url = 'https://ftp.ncbi.nlm.nih.gov/geo/series/GSE132nnn/GSE132044/suppl/GSE132044_pbmc_hg38_cell.tsv.gz'
cells_file = basename(path=cells_url)
curl::curl_download(url=cells_url, destfile=file.path("download", "barcodes.tsv.gz"))

# the features.tsv.gz contains only one column with Ensembl ID and gene symbol concatenated by an underscore - need to be separated
feat = read.table(file.path("download","features.orig.tsv.gz"), header=FALSE)
feat = as.data.frame(stringr::str_split_fixed(string=feat[, 1, drop=TRUE],'_', n=2), stringsAsFactors=FALSE)
feat$V2 = make.unique(feat$V2)
feat$V3 = "Gene Expression"
fh = gzfile(file.path("download", "features.tsv.gz"), open="wb")
write.table(feat, file=fh, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
close(fh)

# now read entire dataset, then split by tech
sc_data = Seurat::Read10X("download")
cell_tech_info = stringr::str_split_fixed(string=colnames(sc_data), pattern="\\.", n=3)
# samples
#
# sample
# "PBMC1" "PBMC2
# techs
# "Smart-seq2" "CEL-Seq2" "10x-Chromium-v2-A" "10x-Chromium-v2-B" "10x-Chromium-v3"   "Drop-seq" "Seq-Well"         "inDrops" "10x-Chromium-v2" 

# smartseq2
dir.create(file.path("counts","smartseq2"), showWarnings=FALSE)
smartseq2_data = as.data.frame(sc_data[,cell_tech_info[,1]=="PBMC1" & cell_tech_info[,2]=="Smart-seq2"])
colnames(smartseq2_data) = gsub(x=gsub(x=colnames(smartseq2_data), pattern="PBMC1.Smart-seq2", replacement="pbmc_smartseq2_sample1", fixed=TRUE), pattern=".p", replacement="_", fixed=TRUE)

col_nms = colnames(smartseq2_data)
smartseq2_data$GeneId = feat$V1
smartseq2_data$GeneSymbol = feat$V2
smartseq2_data = smartseq2_data[,c("GeneId","GeneSymbol", setdiff(colnames(smartseq2_data), c("GeneId","GeneSymbol")))]
smartseq2_data$GeneSymbol = NULL
fh = gzfile(file.path("counts","smartseq2","counts_table.tsv.gz"), open="wb")
readr::write_delim(smartseq2_data, file=fh, delim="\t", col_names=TRUE)
close(fh)

# 10x (10x-Chromium-v3)
dir.create(file.path("counts","10x"), showWarnings=FALSE)
tenx_data = sc_data[,cell_tech_info[,1]=="PBMC1" & cell_tech_info[,2]=="10x-Chromium-v3"]
colnames(tenx_data) = gsub(x=colnames(tenx_data), pattern="PBMC1.10x-Chromium-v3.", replacement="PBMC1_10x_", fixed=TRUE)

# write matrix.mtx.gz, barcodes.tsv.gz and features.tsv.gz
mh = file.path("counts", "10x", "matrix.mtx")
Matrix::writeMM(tenx_data, file=mh)
R.utils::gzip(mh, overwrite=TRUE)

bh = gzfile(file.path("counts", "10x", "barcodes.tsv.gz"), open="wb")
write(colnames(tenx_data), file=bh)
close(bh)

fh = gzfile(file.path("counts", "10x", "features.tsv.gz"), open="wb")
write.table(feat, file=fh, row.names=FALSE, col.names=FALSE, quote=FALSE, sep="\t")
close(fh)

unlink("download", recursive=TRUE)

# write marker file
openxlsx::write.xlsx(data.frame(bcell=c("ENSG00000105369", "ENSG00000156738"),
                                tcell=c("ENSG00000167286", NA),
                                tcell_cd8p=c("ENSG00000153563", "ENSG00000172116"),
                                nk=c("ENSG00000115523", "ENSG00000105374"),
                                myeloid=c("ENSG00000101439", "ENSG00000090382"),
                                monocytes=c("ENSG00000203747", NA),
                                dendritic=c("ENSG00000179639", NA)),
                     "known_markers.xlsx")

# write path_data.csv
write.csv(data.frame(name=c("pbmc_10x", "pbmc_smartseq2"),
                     type=c("10x", "smartseq2"),
                     path=c(normalizePath(file.path("counts", "10x")), normalizePath(file.path("counts", "smartseq2", "counts_table.tsv.gz"))),
                     stats=c(NA, NA)),
          file="path_data.csv",
          row.names=FALSE,
          quote=FALSE)

# write deg_contrasts.xlsx
openxlsx::write.xlsx(data.frame(condition_column=c("orig.ident", "orig.ident", "Phase"),
                                condition_group1=c("pbmc_10x", "pbmc_10x", "G1"),
                                condition_group2=c("pbmc_smartseq2_sample1", "pbmc_smartseq2_sample1", "G2M"),
                                subset_column=c(NA, "seurat_clusters", "seurat_clusters"),
                                subset_group=c(NA, " ", "1;2")),
                     "deg_contrasts.xlsx")
