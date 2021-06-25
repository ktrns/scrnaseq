paramsList = list()

paramsList$project_id = "pbmc"
paramsList$path_data = data.frame(name=c("pbmc_10x","pbmc_smartseq2"),
					type=c("10x","smartseq2"),
					path=c("test_datasets/10x_SmartSeq2_pbmc_GSE132044/counts/10x/", "test_datasets/10x_SmartSeq2_pbmc_GSE132044/counts/smartseq2/counts_table.tsv.gz"),  
					stats=c(NA, NA))
paramsList$downsample_cells_n = NULL
paramsList$path_out = "test_datasets/10x_SmartSeq2_pbmc_GSE132044/results/"
paramsList$file_known_markers = "test_datasets/10x_SmartSeq2_pbmc_GSE132044/known_markers.xlsx"
paramsList$mart_dataset = "hsapiens_gene_ensembl"
paramsList$annot_version = 98
paramsList$annot_main = c(ensembl="ensembl_gene_id", symbol="external_gene_name", entrez="entrezgene_accession")
paramsList$mart_attributes = c(paramsList$annot_main, 
			c("chromosome_name", "start_position", "end_position", 
			  "percentage_gene_gc_content", "gene_biotype", "strand", "description"))
paramsList$biomart_mirror = NULL
paramsList$file_annot = NULL
paramsList$mt = "^MT-"
paramsList$cell_filter = list(nFeature_RNA=c(200, NA), percent_mt=c(NA, 20))
paramsList$feature_filter = list(min_counts=1, min_cells=3) # feature has to be found by at least one count in one cell
paramsList$samples_to_drop = NULL
paramsList$samples_min_cells = 10
paramsList$norm = "RNA"
paramsList$cc_remove = FALSE
paramsList$cc_remove_all = FALSE
paramsList$cc_rescore_after_merge = TRUE
paramsList$vars_to_regress = NULL
paramsList$integrate_samples = list(method="integrate", dimensions=30, reference=NULL, use_reciprocal_pca=FALSE)
paramsList$pc_n = 10
paramsList$cluster_resolution = 0.5
paramsList$marker_padj = 0.05
paramsList$marker_log2FC = log2(2)
paramsList$marker_pct = 0.25
paramsList$latent_vars = NULL
paramsList$deg_contrasts = data.frame(condition_column=c("orig.ident", "orig.ident", "Phase"),
						condition_group1=c("pbmc_10x", "pbmc_10x", "G1"),
						condition_group2=c("pbmc_smartseq2_sample1", "pbmc_smartseq2_sample1", "G2M"),
						subset_column=c(NA, "seurat_clusters", "seurat_clusters"),
						subset_group=c(NA, "", "1;2"),
						downsample_cells_n=c(NA, 50, 30))
paramsList$enrichr_padj = 0.05
paramsList$enrichr_dbs = c("GO_Molecular_Function_2018", "GO_Biological_Process_2018", "GO_Cellular_Component_2018")
paramsList$col = "palevioletred"
paramsList$col_palette_samples = "ggsci::pal_jama"
paramsList$col_palette_clusters = "ggsci::pal_igv"
paramsList$path_to_git = "."
paramsList$debugging_mode = "default_debugging"



rmarkdown::render(
	"scrnaseq.Rmd",
	output_format = "html_document",
	output_dir = ".",
	output_file = "scrnaseq",
	params = paramsList)
	
