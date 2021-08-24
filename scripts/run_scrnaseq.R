# Rscript scripts/run_scrnaseq.R --project-id=pbmc --path-data-csv=test_datasets/10x_SmartSeq2_pbmc_GSE132044/path_data.csv --path-out=test_datasets/10x_SmartSeq2_pbmc_GSE132044/results --deg-contrasts=test_datasets/10x_SmartSeq2_pbmc_GSE132044/deg_contrasts.xlsx --path-to-git=.
#!/usr/bin/env Rscript

# Get arguments from commandline
arguments = commandArgs(trailingOnly=FALSE)

# Path to script
script_index = grep("--file", arguments)
script_dir = dirname(sub("--file=", "", arguments[script_index]))
script_dir = normalizePath(script_dir)
script_name = basename(sub("--file=", "", arguments[script_index]))

# Which packages are needed
library(knitr)
suppressPackageStartupMessages(library("argparse"))

# Prepare parser for arguments
parser = ArgumentParser(
  add_help=TRUE,
  prog=script_name,
  description="Bioinformatics analysis workflow for single-cell RNA-seq analysis. The workflow is based on Seurat, and contains additional visualisations, tables and documentation to better understand the analysis. The workflow supports RNA sequencing data from one or more samples processed with 10X Genomics and SmartSeq-2.",
  epilog="Upon successful completion, the output directory will contain an HTML report, ....")

# Set mutually exclusive groups
analysisGroup = parser$add_mutually_exclusive_group(required=FALSE)
analysisGroup$add_argument("--standard-RNA", action="store_true", help="Run with logNormalization and standard parameters", required=FALSE)
analysisGroup$add_argument("--standard-SCT", action="store_true", help="Run with SCTransform and standard parameters", required=FALSE)

speciesGroup = parser$add_mutually_exclusive_group(required=FALSE)
speciesGroup$add_argument("--human", action="store_true", help="Annotation for human", required=FALSE)
speciesGroup$add_argument("--mouse", action="store_true", help="Annotation for mouse", required=FALSE)
speciesGroup$add_argument("--zebrafish", action="store_true", help="Annotation for zebrafish", required=FALSE)

# Get project and sample description
parser$add_argument("--project-id", action="store", help="Project ID", dest="project_id", required=TRUE)

parser$add_argument("--path-data-csv", action="store", help="CSV file containing the path data of the samples", dest="path_data_csv",required=TRUE)

parser$add_argument("--downsample-cells-n", action="store", type="integer", help="Downsample data to at most N cells (default: #default)", dest="downsample_cells_n", default=0)

parser$add_argument("--path-out", action="store", help="Path to an output directory", dest="path_out", required=TRUE)

# Get individual parameters
parser$add_argument("--file-known-markers", action="store", help="Path to a xlsx file with marker genes based on literature, translated to Ensembl IDs (one list per column, first row as header and Ensembl IDs below); #default if no known marker genes should be plotted", dest="file_known_markers", default=NULL)

parser$add_argument("--mart-dataset", action="store", help="Annotation via biomaRt; Dataset", dest="mart_dataset", required=TRUE)
parser$add_argument("--annot-version", action="store", help="Annotation via biomaRt; Version", dest="annot_version",type="integer", required=TRUE)
parser$add_argument("--annot-main", action="store", help="Annotation via biomaRt; Main (default: #default)", dest="annot_main", default="ensembl=ensembl_gene_id,symbol=external_gene_name,entrez=entrezgene_accession")
parser$add_argument("--mart-attributes", action="store", help="Annotation via biomaRt; Attributes (default: #default)", dest="mart_attributes",default="chromosome_name,start_position,end_position,percentage_gene_gc_content,gene_biotype,strand,description")
parser$add_argument("--biomart-mirror", action="store", help="Annotation via biomaRt; Mirror (default: #default)", dest="biomart_mirror", default=NULL)

parser$add_argument("--file-annot", action="store", help="Alternatively, read a previously compiled annotation table from file (default: #default)", dest="file_annot", default=NULL)

parser$add_argument("--mt", action="store", help="Prefix of mitochondrial genes (default: #default)", dest="mt", default="^MT-")

parser$add_argument("--cell-filter-nFeature-RNA", action="store", help="Filter cells for nFeature_RNA (lower and upper threshold seperated by ',') (default: #default)", dest="cell_filter_nFeature_RNA", default="200, NA")
parser$add_argument("--cell-filter-percent-mt", action="store", help="Filter cells for percent_mt (lower and upper threshold seperated by ',') (default: #default)", dest="cell_filter_percent_mt", default="NA, 25")

parser$add_argument("--feature-filter-min-counts", action="store", help="Filter features for min_counts; Minimum count of one feature (default: #default)", type="integer", dest="feature_filter_min_counts", default=1)
parser$add_argument("--feature-filter-min-cells", action="store", help="Filter features for min_cells; Minimum number of cells in which a feature has to be found (default: #default)", type="integer", dest="feature_filter_min_cells", default=3)

parser$add_argument("--samples-to-drop", action="store", help="Cells from these samples will be dropped after initial QC (name of the datasets_subsamples (e.g. pbmc_NC) seperated by ',') (default: #default)", dest="samples_to_drop", default=NULL)

parser$add_argument("--samples-min-cells", action="store", type="integer", help="Drop samples with too few cells (default: #default)", dest="samples_min_cells", default=10)

parser$add_argument("--norm", action="store", help="Use 'RNA' or 'SCT' for normalisation (default: #default)", dest="norm", default="RNA", choices=c("RNA","SCT"))
parser$add_argument("--cc-remove", action="store_true", help="Remove cell cycle effects (default: #default)", dest="cc_remove", default=FALSE)
parser$add_argument("--cc-remove-all", action="store_true", help="Remove all cell cycle effects, or (if FALSE) only the difference between profilerating cells (G2M and S phase) (default: #default) (Read https://satijalab.org/seurat/v3.1/cell_cycle_vignette.html, for an explanation)", dest="cc_remove_all", default=FALSE)
parser$add_argument("--cc-rescore-after-merge", action="store_true", help="Re-score cell cycle effects after data from different samples have been merged/integrated (default: #default)", dest="cc_rescore_after_merge", default=TRUE)
parser$add_argument("--vars-to-regress", action="store", help="Additional (unwanted) variables that will be regressed out for visualisation and clustering (default: #default)", dest="vars_to_regress", default=NULL)

parser$add_argument("--integrate-samples-method", action="store", help="Single or multiple datasets and how to combine them (method='single','merge', or 'integrate'); 'single': Default when there is only one dataset after filtering, no integration is needed. 'merge': Merge (in other words, concatenate) data when no integration is needed, e.g. when samples were multiplexed on the same chip. 'integrate': Anchors are computed for all pairs of datasets which will give all datasets the same weight during dataset integration but can be computationally intensive. (default: #default)", dest="integrate_samples_method", default="integrate", choices=c("single","merge","integrate"))
parser$add_argument("--integrate-samples-integrate", action="store", help="Additional options for the 'integrate' method (default: #default); dimensions: Number of dimensions to consider for integration. reference: Use one or more (seperated by ',') datasets as reference and compute anchors for all other datasets (computationally faster but less accurate). use_reciprocal_pca: Compute anchors in PCA space (computationally faster but even less accurate). k.filter: How many neighbors to use when filtering anchors (default: min(200, minimum number of cells in a sample)). k.weight: Number of neighbors to consider when weighting anchors (default: min(100, minimum number of cells in a sample)). k.anchor: How many neighbors to use when picking anchors (default: min(5, minimum number of cells in a sample))", dest="integrate_samples_integrate", default="dimensions=30, reference=NULL, use_reciprocal_pca=FALSE")

parser$add_argument("--pc-n", action="store", type="integer", help="The number of PCs to use; adjust this parameter based on the Elbowplot (default: #default)", dest="pc_n", default=10)

parser$add_argument("--cluster-resolution", action="store", type="double", help="Resolution of clusters; low values will lead to fewer clusters of cells (default: #default)", dest="cluster_resolution", default=0.5)

parser$add_argument("--marker-padj", action="store", type="double", help="Thresholds to define marker genes; Adjusted p-value (default: #default)", dest="marker_padj", default=0.05)
parser$add_argument("--marker-log2FC", action="store", type="double", help="Thresholds to define marker genes; Fold change (default: #default)", dest="marker_log2FC", default=log2(2))
parser$add_argument("--marker-pct", action="store", type="double", help="Thresholds to define marker genes; Percentage of marker expressing cells (default: #default)", dest="marker_pct", default=0.25)

parser$add_argument("--latent-vars", action="store", help="Additional (unwanted) variables to account for in statistical tests (default: #default)", dest="latent_vars", default=NULL)

parser$add_argument("--deg-contrasts", action="store", help="Contrasts to find differentially expressed genes; xlsx file (Required columns: condition_column, condition_group1, condition_group2; Optional columns: subset_column, subset_group, assay, slot, padj, log2FC, min_pct, test, downsample_cells_n, latent_vars)", dest="deg_contrasts", required=TRUE)

parser$add_argument("--enrichr-padj", action="store", type="double", help="P-value threshold for functional enrichment tests (default: #default)", dest="enrichr_padj", default=0.05)

parser$add_argument("--enrichr-dbs", action="store", help="Enrichr databases of interest seperated by ',' (default: #default)", dest="enrichr_dbs", default="GO_Molecular_Function_2018,GO_Biological_Process_2018,GO_Cellular_Component_2018")

parser$add_argument("--col", action="store", help="Main colour to use for plots (default: #default)", dest="col", default="palevioletred")
parser$add_argument("--col-palette-samples", action="store", help="Colour palette and colours used for samples (default: #default)", dest="col_palette_samples", default="ggsci::pal_jama")
parser$add_argument("--col-palette-clusters", action="store", help="Colour palette and colours used for cluster (default: #default)", dest="col_palette_clusters", default="ggsci::pal_igv")

parser$add_argument("--path-to-git", action="store", help="Path to the git repository (default: #default)", dest="path_to_git", default=".")

parser$add_argument("--debugging-mode", action="store", help="Debugging mode: 'default_debugging' for default, 'terminal_debugger' for debugging without X11, 'print_traceback' for non-interactive sessions (default: #default)", dest="debugging_mode", default="print_traceback", choices=c("default_debugging","terminal_debugger","print_traceback"))

parser$add_argument("--report-name",  action="store", help="Name of the HTML report (default: #default)", default="scrnaseq.html", dest="report_name")
parser$add_argument("--verbose",  action="store_true", help="Be verbose (default: #default)", default=TRUE, dest="verbose")

opt = parser$parse_args()

#
# Check params
#
param = list()

# Choose sets of parameters
if (opt[["standard_RNA"]]) {
  param[["norm"]] = "RNA"
  param[["vars_to_regress"]] = c("nCount_RNA", "percent_mt")
  param[["pc_n"]] = 10
  param[["cluster_resolution"]] = 0.5
}

if (opt[["standard_SCT"]]) {
  param[["norm"]] = "SCT"
  param[["vars_to_regress"]] = c("percent_mt")
  param[["pc_n"]] = 30
  param[["cluster_resolution"]] = 0.5
}

if (opt[["human"]]) {
  param[["mart_dataset"]] = "hsapiens_gene_ensembl"
  param[["annot_version"]] = 98
  param[["annot_main"]] = c(ensembl="ensembl_gene_id", symbol="external_gene_name", entrez="entrezgene_accession")
  param[["mart_attributes"]] = c(param$annot_main, c("chromosome_name", "start_position", "end_position",
                                                     "percentage_gene_gc_content", "gene_biotype", "strand", "description"))
  param[["mt"]] = "^MT-"
  param[["enrichr_dbs"]] = c("GO_Biological_Process_2018", "KEGG_2019_Human", "WikiPathways_2019_Human", "GO_Molecular_Function_2018", "GO_Cellular_Component_2018")
}

if (opt[["mouse"]]) {
  param[["mart_dataset"]] = "mmusculus_gene_ensembl"
  param[["annot_version"]] = 103
  param[["annot_main"]] = c(ensembl="ensembl_gene_id", symbol="external_gene_name", entrez="entrezgene_accession")
  param[["mart_attributes"]] = c(param$annot_main, c("chromosome_name", "start_position", "end_position",
                                                     "percentage_gene_gc_content", "gene_biotype", "strand", "description"))
  param[["mt"]] = "^mt-"
  param[["enrichr_dbs"]] = c("GO_Biological_Process_2018", "KEGG_2019_Mouse", "WikiPathways_2019_Mouse", "GO_Molecular_Function_2018", "GO_Cellular_Component_2018")
}

if (opt[["zebrafish"]]) {
  param[["mart_dataset"]] = "drerio_gene_ensembl"
  param[["annot_version"]] = 103
  param[["annot_main"]] = c(ensembl="ensembl_gene_id", symbol="external_gene_name", entrez="external_gene_name")
  param[["mart_attributes"]] = c(param$annot_main, c("chromosome_name", "start_position", "end_position",
                                                     "percentage_gene_gc_content", "gene_biotype", "strand", "description"))
  param[["mt"]] = "^mt-"
  param[["enrichr_dbs"]] = c("GO_Biological_Process_2018", "GO_Molecular_Function_2018", "GO_Cellular_Component_2018")
}

# Write project and sample description
# Param project id
param[["project_id"]] = as.character(opt[["project_id"]])

# Param path_data
# Read path_data csv file and check values
if (!file.exists(opt[["path_data_csv"]]))
  stop("Path to a path data definition file specified via --path-data-csv does not exist or is not a file!")
param[["path_data"]] = read.csv(normalizePath(opt[["path_data_csv"]]), header = TRUE, sep = ",", quote = "\"",
                                dec = ".", fill = TRUE, comment.char = "")

if (sum(!names(param$path_data) %in% c("name", "type", "path", "stats")) > 0)
  stop("Need to provide all columns (name, type, path, stats) in path data definition file specified via --path-data-csv!")

if (sum(!param$path_data$type %in% c("10x","smartseq2")) > 0)
  stop("Need to provide valid types ('10x' or 'smartseq2') in path-data-csv!")

if (sum(!file.exists(param$path_data$path)) > 0) 
  stop("At least one path to input data specified in --path-data-csv does not exist!")

if (sum(!is.na(param$path_data$stats)) > 0) {
  if (sum(!file.exists(param$path_data$stats[!is.na(param$path_data$stats)])) > 0) 
    stop("At least one path to a metrics_summary.csv file specified in --path-data-csv does not exist!")
}

# Param downsample_cells_n
param[["downsample_cells_n"]] = opt[["downsample_cells_n"]]

# Param path_out
param[["path_out"]] = normalizePath(opt[["path_out"]])

# Change additional parameters individually
# Param file_known_markers
if (!is.null(opt[["file_known_markers"]])) {
  if (!file.exists(opt[["file_known_markers"]])) 
    stop("Path to a known markers file specified via --file-known-markers does not exist or is not a file!")
  param[["file_known_markers"]] = normalizePath(opt[["file_known_markers"]])
} else {
  param[["file_known_markers"]] = NULL
}

# Param mart_dataset
param[["mart_dataset"]] = opt[["mart_dataset"]]

# Param annot_version
param[["annot_version"]] = opt[["annot_version"]]

# Param annot_main
param[["annot_main"]] = trimws(unlist(strsplit(trimws(unlist(strsplit(opt[["annot_main"]], ","))), "=")))[c(F,T)]
names(param[["annot_main"]]) = trimws(unlist(strsplit(trimws(unlist(strsplit(opt[["annot_main"]], ","))), "=")))[c(T,F)]

# Param mart_attributes
param[["mart_attributes"]] = unique(c(param[["annot_main"]], trimws(unlist(strsplit(opt[["mart_attributes"]], ",")))))

# Param biomart_mirror
param[["biomart_mirror"]] = opt[["biomart_mirror"]]

# Param file_annot
if (!is.null(opt[["file_annot"]])) {
  if (!file.exists(opt[["file_annot"]])) 
    stop("Path to annotation file specified via --file-annot does not exists or is not a file!")
  param[["file_annot"]] = normalizePath(opt[["file_annot"]])
} else {
  param[["file_annot"]] = NULL
}

# Param mt
param[["mt"]] = opt[["mt"]]

# Param cell_filter
param[["cell_filter"]][["nFeature_RNA"]] = as.numeric(trimws(unlist(strsplit(opt[["cell_filter_nFeature_RNA"]], ","))))
param[["cell_filter"]][["percent_mt"]] = as.numeric(trimws(unlist(strsplit(opt[["cell_filter_percent_mt"]], ","))))

# Param feature_filter
param[["feature_filter"]][["min_counts"]] = opt[["feature_filter_min_counts"]]
param[["feature_filter"]][["min_cells"]] = opt[["feature_filter_min_cells"]]
 
# Param samples_to_drop
if (!is.null(opt[["samples_to_drop"]])) {
  param[["samples_to_drop"]] = trimws(unlist(strsplit(opt[["samples_to_drop"]], ",")))
} else {
  param[["samples_to_drop"]] = opt[["samples_to_drop"]]
}
# 
# Param samples_min_cells
param[["samples_min_cells"]] = opt[["samples_min_cells"]]

# Param norm
param[["norm"]] = opt[["norm"]]

# Param cc_remove
param[["cc_remove"]] = opt[["cc_remove"]] 

# Param cc_remove_all
param[["cc_remove_all"]] = opt[["cc_remove_all"]]

# Param cc_rescore_after_merge
if ("cc_rescore_after_merge" %in% names(opt)) 
	param[["cc_rescore_after_merge"]] = opt[["cc_rescore_after_merge"]]

# Param vars_to_regress
if (!is.null(opt[["vars_to_regress"]])) { 
  param[["vars_to_regress"]] = trimws(unlist(strsplit(opt[["vars_to_regress"]], ",")))
} else {
  param[["vars_to_regress"]] = opt[["vars_to_regress"]]
}

# Param integrate_samples
if (opt[["integrate_samples_method"]] == "integrate") {
  param[["integrate_samples"]] = list(method=opt[["integrate_samples_method"]])
  integrate_list = trimws(unlist(strsplit(opt[["integrate_samples_integrate"]], ",")))
  if (!is.na(grep("dimensions",integrate_list)[1])) {
    param[["integrate_samples"]]["dimensions"] = as.numeric(trimws(unlist(strsplit(integrate_list[grep("dimensions",integrate_list)[1]],"=")))[2])
  }
  if (!is.na(grep("reference",integrate_list)[1])) {
    param[["integrate_samples"]]["reference"] = trimws(unlist(strsplit(integrate_list[grep("reference",integrate_list)[1]],"=")))[2]
	if (param[["integrate_samples"]]["reference"] == "NULL") {
	  param[["integrate_samples"]]["reference"] = NULL
	} else {
	  param[["integrate_samples"]]["reference"] = trimws(unlist(strsplit(param[["integrate_samples"]]["reference"], ",")))
	}
  }
  if (!is.na(grep("use_reciprocal_pca",integrate_list)[1])) {
    param[["integrate_samples"]]["use_reciprocal_pca"] = as.logical(trimws(unlist(strsplit(integrate_list[grep("use_reciprocal_pca",integrate_list)[1]],"=")))[2])
  }
  if (!is.na(grep("k.filter",integrate_list)[1])) {
    param[["integrate_samples"]]["k.filter"] = as.numeric(trimws(unlist(strsplit(integrate_list[grep("k.filter",integrate_list)[1]],"=")))[2])
  }
  if (!is.na(grep("k.weight",integrate_list)[1])) {
    param[["integrate_samples"]]["k.weight"] = as.numeric(trimws(unlist(strsplit(integrate_list[grep("k.weight",integrate_list)[1]],"=")))[2])
  }
  if (!is.na(grep("k.anchor",integrate_list)[1])) {
    param[["integrate_samples"]]["k.anchor"] = as.numeric(trimws(unlist(strsplit(integrate_list[grep("k.anchor",integrate_list)[1]],"=")))[2])
  }
} else {
  param[["integrate_samples"]] = list(method=opt[["integrate_samples_method"]])
}

# Param pc_n
param[["pc_n"]] = opt[["pc_n"]]

# Param cluster_resolution
param[["cluster_resolution"]] = opt[["cluster_resolution"]]

# Param marker_padj
param[["marker_padj"]] = opt[["marker_padj"]]

# Param marker_log2FC
param[["marker_log2FC"]] = opt[["marker_log2FC"]]

# Param marker_pct
param[["marker_pct"]] = opt[["marker_pct"]]

# Param latent_vars
if (!is.null(opt[["latent_vars"]])) { 
  param[["latent_vars"]] = trimws(unlist(strsplit(opt[["latent_vars"]], ",")))
} else {
  param[["latent_vars"]] = opt[["latent_vars"]]
}

# deg_contrasts
param[["deg_contrasts"]] = opt[["deg_contrasts"]]

# Param enrichr_padj
param[["enrichr_padj"]] = opt[["enrichr_padj"]]

# Param enrichr_dbs
param[["enrichr_dbs"]] = trimws(unlist(strsplit(opt[["enrichr_dbs"]], ",")))

# Param col
if (!opt[["col"]] %in% colours()) {
  stop("Only R colours allowed for --col!")
} else {
  param[["col"]] = opt[["col"]]
}

# Param col_palette_samples
param[["col_palette_samples"]] = opt[["col_palette_samples"]]

# Param col_palette_clusters
param[["col_palette_clusters"]] = opt[["col_palette_clusters"]]

# Param path_to_git
param[["path_to_git"]] = opt[["path_to_git"]]

# Param debugging_mode
param[["debugging_mode"]] = opt[["debugging_mode"]]


#
# Run run knitr
#
# Note: a small R script will be generated so that knitr can be run again if needed

# The name of the script
script_name = paste0(opt[["report_name"]], ".rerun.r")

# This is the template for the script
script_template = '
#!/usr/bin/env Rscript

library(knitr)

rmarkdown::render(
    file.path("{{rmd_dir}}","scrnaseq.Rmd"),
    output_format="html_document",
    output_dir="{{output_dir}}",
    intermediates_dir="{{intermediates_dir}}",
    output_file="{{output_file}}",
    knit_root_dir="{{knit_root_dir}}",
    quiet={{quiet}},
    params={{params_lst}})
'

# Now fill the script with content
script_content = knitr::knit_expand(text=script_template,
                                    rmd_dir=param[["path_to_git"]],
                                    output_dir=param[["path_out"]],
                                    intermediates_dir=param[["path_out"]],
                                    output_file=opt[["report_name"]],
                                    knit_root_dir=param[["path_out"]],
                                    quiet=!opt[["verbose"]],
                                    params_lst=deparse(param, control="all"))


writeLines(script_content, file.path(param[["path_out"]], script_name))

# Now run actual knitr process
rmarkdown::render(
    file.path(param[["path_to_git"]],"scrnaseq.Rmd"),
    output_format="html_document",
    output_dir=param[["path_out"]],
    intermediates_dir=param[["path_out"]],
    output_file=opt[["report_name"]],
    knit_root_dir=param[["path_out"]],
    quiet=!opt[["verbose"]],
    params=param)
