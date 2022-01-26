#!/usr/bin/env Rscript
# Rscript scripts/run_scrnaseq.R --project-id=pbmc --path-data-csv=test_datasets/10x_SmartSeq2_pbmc_GSE132044/path_data.csv --path-out=test_datasets/10x_SmartSeq2_pbmc_GSE132044/results --mart-dataset=hsapiens_gene_ensembl --annot-version=98 --deg-contrasts=test_datasets/10x_SmartSeq2_pbmc_GSE132044/deg_contrasts.xlsx --path-to-git=.

# options(rlang_trace_top_env = rlang::caller_env())
# options(error = function() {
#   sink()
#   print(rlang::trace_back(bottom = sys.frame(-1)), simplify = "none")
# })

# Required arguments are: 
# --path-data-csv
# --path-out
# --mart-dataset
# --annot-version

# Get command line call, identify script name and directory and extract the actual arguments; if there are no, add a "help" argument
arguments = commandArgs(trailingOnly=FALSE)
script_index = grep("--file",arguments)
script_dir = dirname(sub("--file=","",arguments[script_index]))
script_dir = normalizePath(script_dir)
script_name = basename(sub("--file=", "", arguments[script_index]))

arguments_length = length(arguments)
arguments_index = grep("--args",arguments)[1]
if(is.na(arguments_index)){
  arguments = c("--args", "--help")
} else{
  arguments = arguments[(arguments_index+1):arguments_length]
}

# Which packages are needed
library(knitr)
library(magrittr) # %>% operator
suppressWarnings(suppressPackageStartupMessages(library("argparse")))

# Prepare parser for arguments (help: https://docs.python.org/3/library/argparse.html#the-add-argument-method)
parser = ArgumentParser(
  add_help=TRUE,
  prog=script_name,
  description="Bioinformatics analysis workflow for single-cell RNA-seq analysis. The workflow is based on Seurat, and contains additional visualisations, tables and documentation to better understand the analysis. The workflow supports RNA sequencing data from one or more samples processed with 10X Genomics and SmartSeq-2.",
  epilog="The workflow analyses single-cell RNA-seq data in R. An extensive HTML report is generated, and data is exported from R so you can further explore it by yourself."
)

# Set mutually exclusive groups
analysisGroup = parser$add_mutually_exclusive_group(required=FALSE)
analysisGroup$add_argument(
  "--standard-RNA",
  action="store_true",
  help="Run log normalisation and apply a set of standard parameters",
  required=FALSE
)

analysisGroup$add_argument(
  "--standard-SCT",
  action="store_true",
  help="Run SCTransform and apply a set of standard parameters",
  required=FALSE
)

speciesGroup = parser$add_mutually_exclusive_group(required=FALSE)
speciesGroup$add_argument(
  "--human",
  action="store_true",
  help="Annotation for human",
  required=FALSE
)

speciesGroup$add_argument(
  "--mouse",
  action="store_true",
  help="Annotation for mouse",
  required=FALSE
)

speciesGroup$add_argument(
  "--zebrafish",
  action="store_true",
  help="Annotation for zebrafish",
  required=FALSE
)

# Get project Id
parser$add_argument(
  "--project-id",
  action="store",
  help="Project ID",
  dest="project_id",
  default="Project"
)

parser$add_argument(
  "--path-data-csv",
  action="store",
  help="CSV file containing the paths to the samples",
  dest="path_data_csv",
  required=TRUE
)

parser$add_argument(
  "--assay-raw",
  action="store",
  help="Single-cell RNA-seq, or spatial RNA-seq data (default: %(default)s)",
  dest="assay_raw",
  default="RNA",
  choices=c("RNA", "Spatial")
)

parser$add_argument(
  "--downsample-cells-n",
  action="store",
  type="integer",
  help="Downsample data to at most N cells (default: %(default)s)",
  dest="downsample_cells_n",
  default=NULL
)

parser$add_argument(
  "--path-out",
  action="store",
  help="Path to the output directory",
  dest="path_out",
  required=TRUE
)

# Get individual parameters
parser$add_argument(
  "--file-known-markers",
  action="store",
  help="Path to an xlsx file with marker genes based on literature, translated to Ensembl IDs (one list per column, first row as header and Ensembl IDs below); #default if no known marker genes should be plotted",
  dest="file_known_markers",
  default=NULL
)

parser$add_argument(
  "--mart-dataset",
  action="store",
  help="Annotation via biomaRt; Dataset",
  dest="mart_dataset",
  required=TRUE
)

parser$add_argument(
  "--annot-version",
  action="store",
  help="Annotation via biomaRt; Version",
  dest="annot_version",
  type="integer",
  required=TRUE
)

parser$add_argument(
  "--annot-main",
  action="store",
  nargs="*",
  help="Annotation via biomaRt; Main (default: %(default)s)",
  dest="annot_main",
  default=c("ensembl=ensembl_gene_id", "symbol=external_gene_name", "entrez=entrezgene_accession")
)

parser$add_argument(
  "--mart-attributes",
  action="store",
  nargs="*",
  help="Annotation via biomaRt; Attributes (default: %(default)s)",
  dest="mart_attributes",
  default=c("chromosome_name", "start_position", "end_position", "percentage_gene_gc_content", "gene_biotype", "strand", "description")
)

parser$add_argument(
  "--biomart-mirror",
  action="store",
  help="Annotation via biomaRt; Mirror (default: %(default)s)",
  dest="biomart_mirror",
  default=NULL
)

parser$add_argument(
  "--mt",
  action="store",
  help="Prefix of mitochondrial genes (default: %(default)s)",
  dest="mt",
  default="^MT-"
)

parser$add_argument(
  "--cell-filter-nCount",
  action="store",
  nargs=1,
  help="Filter cells for number of counts (lower and upper threshold separated by ', ') (default: %(default)s)",
  dest="cell_filter_nCount",
  default=NULL
)

parser$add_argument(
  "--cell-filter-nFeature",
  action="store",
  nargs=1,
  help="Filter cells for number of features (lower and upper threshold separated by ', ') (default: %(default)s)",
  dest="cell_filter_nFeature",
  default=NULL
)

parser$add_argument(
  "--cell-filter-percent-mt",
  action="store",
  nargs=1,
  help="Filter cells for percent mitochondrial (lower and upper threshold seperated by ', ') (default: %(default)s)",
  dest="cell_filter_percent_mt",
  default=NULL
)

parser$add_argument(
  "--cell-filter",
  action="store",
  nargs="*",
  help="Filter cells using one or more cell metadata columns (column_name=min,max or column_name=category1,category2,category3,...) (default: %(default)s)",
  dest="cell_filter",
  default=NULL
)

parser$add_argument(
  "--feature-filter-min-counts",
  action="store",
  help="Filter features for min_counts; Minimum count of one feature (default: %(default)s)",
  type="integer",
  dest="feature_filter_min_counts",
  default=1
)

parser$add_argument(
  "--feature-filter-min-cells",
  action="store",
  help="Filter features for min_cells; Minimum number of cells in which a feature has to be found (default: %(default)s)",
  type="integer",
  dest="feature_filter_min_cells",
  default=3
)

parser$add_argument(
  "--samples-to-drop",
  action="store",
  nargs="*",
  help="Cells from these samples will be dropped after initial QC (name of the datasets_subsamples) (default: %(default)s)",
  dest="samples_to_drop",
  default=NULL
)

parser$add_argument(
  "--samples-min-cells",
  action="store",
  type="integer",
  help="Drop samples with less than N cells (default: %(default)s)",
  dest="samples_min_cells",
  default=10
)

parser$add_argument(
  "--norm",
  action="store",
  help="Use 'RNA' or 'SCT' for normalisation (default: %(default)s)",
  dest="norm",
  default="RNA",
  choices=c("RNA", "SCT")
)

parser$add_argument(
  "--cc-remove",
  action="store_true",
  help="Remove cell cycle effects (default: %(default)s)",
  dest="cc_remove",
  default=FALSE
)

parser$add_argument(
  "--cc-remove-all",
  action="store_true",
  help="Remove all cell cycle effects, or (if FALSE) only the difference between profilerating cells (G2M and S phase) (default: #default) (Read https://satijalab.org/seurat/v3.1/cell_cycle_vignette.html, for an explanation)",
  dest="cc_remove_all",
  default=FALSE
)

parser$add_argument(
  "--cc-rescore-after-merge",
  action="store_true",
  help="Re-score cell cycle effects after data from different samples have been merged/integrated (default: %(default)s)",
  dest="cc_rescore_after_merge",
  default=TRUE
)

parser$add_argument(
  "--vars-to-regress",
  action="store",
  nargs="*",
  help="Additional (unwanted) variables that will be regressed out (default: %(default)s)",
  dest="vars_to_regress",
  default=NULL
)

parser$add_argument(
  "--combine-samples-method",
  action="store",
  help="Single or multiple datasets and how to combine them (method='single', 'merge', or 'integrate'); 'single': Default when there is only one dataset after filtering, no integration is needed. 'merge': Merge (in other words, concatenate) data when no integration is needed, e.g. when samples were multiplexed on the same chip. 'integrate': Anchors are computed for all pairs of datasets which will give all datasets the same weight during dataset integration but can be computationally intensive. (default: %(default)s)",
  dest="combine_samples_method",
  default="integrate",
  choices=c("single", "merge", "integrate")
)

parser$add_argument(
  "--integrate-dimensions",
  action="store",
  type="integer",
  help="Number of dimensions to consider for integration (default: %(default)s)",
  dest="integrate_samples_dimensions",
  default=30
)

parser$add_argument(
  "--integrate-reference",
  action="store",
  type="character",
  help="Use one or more (seperated by ', ') datasets as reference and compute anchors for all other datasets; computationally faster but less accurate (default: %(default)s)",
  dest="integrate_samples_reference",
  default=NULL
)

parser$add_argument(
  "--integrate-use-reciprocal-pca",
  action="store_true",
  help="Set if anchors should be computed in PCA space; computationally faster but even less accurate (default: %(default)s)",
  dest="integrate_samples_use_reciprocal_pca",
  default=FALSE
)

parser$add_argument(
  "--integrate-k-filter",
  action="store",
  type="integer",
  help="How many neighbors to use when filtering anchors (default: min(%(default)s, minimum number of cells in a sample))",
  dest="integrate_samples_k_filter",
  default=200
)

parser$add_argument(
  "--integrate-k-weight",
  action="store",
  type="integer",
  help="Number of neighbors to consider when weighting anchors (default: min(%(default)s, minimum number of cells in a sample))",
  dest="integrate_samples_k_weight",
  default=100
)

parser$add_argument(
  "--integrate-k-anchor",
  action="store",
  type="integer",
  help="How many neighbors to use when picking anchors (default: min(%(default)s, minimum number of cells in a sample))",
  dest="integrate_samples_k_anchor",
  default=5
)

parser$add_argument(
  "--integrate-k-score",
  action="store",
  type="integer",
  help="How many neighbors to use when scoring anchors (default: min(%(default)s, minimum number of cells in a sample))",
  dest="integrate_samples_k_score",
  default=30
)

parser$add_argument(
  "--pc-n",
  action="store",
  type="integer",
  help="The number of PCs to use; adjust this parameter based on the Elbow plot (default: %(default)s)",
  dest="pc_n",
  default=10
)

parser$add_argument(
  "--cluster-k",
  action="store",
  type="integer",
  help="The number of k nearest neighbors to find clusters (default: %(default)s)",
  dest="cluster_k",
  default=20
)

parser$add_argument(
  "--umap-k",
  action="store",
  type="integer",
  help="The number of k nearest neighbors to construct the UMAP (default: %(default)s)",
  dest="umap_k",
  default=30
)

parser$add_argument(
  "--cluster-resolution",
  action="store",
  type="double",
  help="Resolution of clusters; low values will lead to fewer clusters of cells (default: %(default)s)",
  dest="cluster_resolution",
  default=1
)

parser$add_argument(
  "--cluster-resolution-test",
  nargs="*",
  action="store",
  type="double",
  help="Other clusters resolutions to test; will not be used for further analysis (default: %(default)s)",
  dest="cluster_resolution_test",
  default=NULL
)

parser$add_argument(
  "--marker-padj",
  action="store",
  type="double",
  help="Thresholds to define marker genes; adjusted p-value (default: %(default)s)",
  dest="marker_padj",
  default=0.05
)

parser$add_argument(
  "--marker-log2FC",
  action="store",
  type="double",
  help="Thresholds to define marker genes; fold change (default: %(default)s)",
  dest="marker_log2FC",
  default=log2(2)
)

parser$add_argument(
  "--marker-pct",
  action="store",
  type="double",
  help="Thresholds to define marker genes; only genes expressed in at least N%% of cells in each group are tested (default: %(default)s)",
  dest="marker_pct",
  default=0.25
)

parser$add_argument(
  "--latent-vars",
  action="store",
  nargs="*",
  help="Additional (unwanted) variables to account for in statistical tests (default: %(default)s)",
  dest="latent_vars",
  default=NULL
)

parser$add_argument(
  "--deg-contrasts",
  action="store",
  help="Contrasts to find differentially expressed genes; xlsx file (Required columns: condition_column, condition_group1, condition_group2; Optional columns: subset_column, subset_group, assay, slot, padj, log2FC, min_pct, test, downsample_cells_n, latent_vars)",
  dest="deg_contrasts"
)

parser$add_argument(
  "--enrichr-padj",
  action="store",
  type="double",
  help="P-value threshold for functional enrichment tests (default: %(default)s)",
  dest="enrichr_padj",
  default=0.05
)

parser$add_argument(
  "--enrichr-dbs",
  action="store",
  nargs="*",
  help="Enrichr databases of interest seperated by ', ' (default: %(default)s)",
  dest="enrichr_dbs",
  default=c("GO_Molecular_Function_2018", "GO_Biological_Process_2018", "GO_Cellular_Component_2018", "Azimuth_Cell_Types_2021", "CellMarker_Augmented_2021", "Descartes_Cell_Types_and_Tissue_2021")
)

parser$add_argument(
  "--enrichr-site",
  action="store",
  help="Enrichr site (default: %(default)s)",
  dest="enrichr_site",
  default="Enrichr",
  choices=c("Enrichr", "FlyEnrichr", "WormEnrichr", "YeastEnrichr", "FishEnrichr")
)

parser$add_argument(
  "--col",
  action="store",
  help="Main colour to use for plots (default: %(default)s)",
  dest="col",
  default="palevioletred"
)

parser$add_argument(
  "--col-palette-samples",
  action="store",
  help="Colour palette and colours used for samples (default: %(default)s)",
  dest="col_palette_samples",
  default="ggsci::pal_jama"
)

parser$add_argument(
  "--col-palette-clusters",
  action="store",
  help="Colour palette and colours used for cluster (default: %(default)s)",
  dest="col_palette_clusters",
  default="ggsci::pal_igv"
)

parser$add_argument(
  "--path-to-git",
  action="store",
  help="Path to the git repository (default: %(default)s)",
  dest="path_to_git",
  default=normalizePath(file.path(script_dir, "/.."))
)

parser$add_argument(
  "--debugging-mode",
  action="store",
  help="Debugging mode: 'default_debugging' for default, 'terminal_debugger' for debugging without X11, 'print_traceback' for non-interactive sessions (default: %(default)s)",
  dest="debugging_mode",
  default="default_debugging",
  choices=c("default_debugging", "terminal_debugger", "print_traceback")
)

parser$add_argument(
  "--cores",
  action="store",
  type="integer",
  help="The number of cores to use for parallel computations (default: %(default)s)",
  dest="cores",
  default=4
)

parser$add_argument(
  "--report-name",
  action="store",
  help="Name of the HTML report (default: %(default)s)",
  default="scrnaseq.html",
  dest="report_name"
)

parser$add_argument(
  "--verbose",
  action="store_true",
  help="Be verbose (default: %(default)s)",
  default=TRUE,
  dest="verbose"
)

# Parse arguments
opt = parser$parse_args(arguments)

#
# Check and process parameters
#
param = list()

# Choose sets of parameters (RCUG standards)
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
  param[["enrichr_dbs"]] = c("GO_Biological_Process_2018", "KEGG_2019_Human", "WikiPathways_2019_Human", "GO_Molecular_Function_2018", "GO_Cellular_Component_2018", "Azimuth_Cell_Types_2021", "CellMarker_Augmented_2021", "Descartes_Cell_Types_and_Tissue_2021")
}

if (opt[["mouse"]]) {
  param[["mart_dataset"]] = "mmusculus_gene_ensembl"
  param[["annot_version"]] = 103
  param[["annot_main"]] = c(ensembl="ensembl_gene_id", symbol="external_gene_name", entrez="entrezgene_accession")
  param[["mart_attributes"]] = c(param$annot_main, c("chromosome_name", "start_position", "end_position", 
                                                     "percentage_gene_gc_content", "gene_biotype", "strand", "description"))
  param[["mt"]] = "^mt-"
  param[["enrichr_dbs"]] = c("GO_Biological_Process_2018", "KEGG_2019_Mouse", "WikiPathways_2019_Mouse", "GO_Molecular_Function_2018", "GO_Cellular_Component_2018", "Azimuth_Cell_Types_2021", "CellMarker_Augmented_2021", "Descartes_Cell_Types_and_Tissue_2021")
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

### ****************************************************
### Walk through all parameters and process if necessary
### ****************************************************
param[["project_id"]] = as.character(opt[["project_id"]])

# Read path_data csv file and check values; make paths absolute
if (!file.exists(opt[["path_data_csv"]]))
  stop("Path to a path data definition file specified via --path-data-csv does not exist or is not a file!")
param[["path_data"]] = read.csv(normalizePath(opt[["path_data_csv"]]), 
                                header=TRUE, sep=",", quote="\"", dec=".", fill=TRUE, comment.char="")

if (sum(!names(param$path_data) %in% c("name", "type", "path", "stats")) > 0)
  stop("Need to provide all columns (name, type, path, stats) in path data definition file specified via --path-data-csv!")

if (sum(!param$path_data$type %in% c("10x", "smartseq2")) > 0)
  stop("Need to provide valid types ('10x' or 'smartseq2') in path-data-csv!")

if (sum(!file.exists(param$path_data$path)) > 0) 
  stop("At least one path to input data specified in --path-data-csv does not exist!")

if (sum(!is.na(param$path_data$stats)) > 0) {
  if (sum(!file.exists(param$path_data$stats[!is.na(param$path_data$stats)])) > 0) 
    stop("At least one path to a metrics_summary.csv file specified in --path-data-csv does not exist!")
}

param$path_data$path = normalizePath(param$path_data$path)
not_na = which(!is.na(param$path_data$stats))
if (length(not_na) > 0) param$path_data$stats[not_na] = normalizePath(param$path_data$stats[not_na])

# Basic parameters
param[["assay_raw"]] = opt[["assay_raw"]]
if (!is.null(opt[["downsample_cells_n"]])) {
    param[["downsample_cells_n"]] = opt[["downsample_cells_n"]]
} else {
    param["downsample_cells_n"] = list(NULL)
}

param[["path_out"]] = opt[["path_out"]]

# Known marker genes 
if (!is.null(opt[["file_known_markers"]])) {
  if (!file.exists(opt[["file_known_markers"]])) 
    stop("Path to a known markers file specified via --file-known-markers does not exist or is not a file!")
  param[["file_known_markers"]] = normalizePath(opt[["file_known_markers"]])
} else {
  param["file_known_markers"] = list(NULL)
}

# Annotation parameters
param[["mart_dataset"]] = opt[["mart_dataset"]]
param[["annot_version"]] = opt[["annot_version"]]
annot_main_tbl = mapply(function(x) strsplit(x, split="=", fixed=TRUE)[[1]], opt[["annot_main"]])
param[["annot_main"]] = annot_main_tbl[2, ]
names(param[["annot_main"]]) = annot_main_tbl[1, ]
param[["mart_attributes"]] = c(param[["annot_main"]], opt[["mart_attributes"]]) %>% trimws() %>% unique()
if (!is.null(opt[["biomart_mirror"]])) {
    param[["biomart_mirror"]] = opt[["biomart_mirror"]]
} else {
    param["biomart_mirror"] = list(NULL)
}
param[["mt"]] = opt[["mt"]]

# Filter
cell_filter = list()
if (!is.null(opt[["cell_filter_nCount"]])) cell_filter[[paste0("nCount_", param[["assay_raw"]])]] = opt[["cell_filter_nCount"]] %>% strsplit(",") %>% unlist() %>% trimws() %>% gsub(pattern="NA", replacement=NA) %>% as.numeric()
if (!is.null(opt[["cell_filter_nFeature"]])) cell_filter[[paste0("nFeature_", param[["assay_raw"]])]] = opt[["cell_filter_nFeature"]] %>% strsplit(",") %>% unlist() %>% trimws() %>% gsub(pattern="NA", replacement=NA) %>% as.numeric()
if (!is.null(opt[["cell_filter_percent_mt"]])) cell_filter[["percent_mt"]] = opt[["cell_filter_percent_mt"]] %>% strsplit(",") %>% unlist() %>% trimws() %>% gsub(pattern="NA", replacement=NA) %>% as.numeric()

if (!is.null(opt[["cell_filter"]])) {
  filter_tbl = mapply(function(x) strsplit(x, split="=", fixed=TRUE)[[1]], opt[["cell_filter"]])
  filter_names = filter_tbl[1, ] %>% unname() %>% trimws()
  filter_values = lapply(filter_tbl[2, ], function(f) {
      f = f %>% strsplit(",", fixed=TRUE) %>% unlist() %>% trimws() %>% gsub(pattern="NA", replacement=NA)
      f = ifelse (length(f)==2 & all(!is.na(as.numeric(f))), as.numeric(f), as.character(f))
      return(f)
  })
  cell_filter[filter_names] = filter_values
}
param[["cell_filter"]] = cell_filter

param[["feature_filter"]][["min_counts"]] = opt[["feature_filter_min_counts"]]
param[["feature_filter"]][["min_cells"]] = opt[["feature_filter_min_cells"]]

if (!is.null(opt[["samples_to_drop"]])) {
    param[["samples_to_drop"]] = opt[["samples_to_drop"]]
} else {
    param["samples_to_drop"] = list(NULL)
}
param[["samples_min_cells"]] = opt[["samples_min_cells"]]

# Pre-processing
param[["norm"]] = opt[["norm"]]
param[["cc_remove"]] = opt[["cc_remove"]] 
param[["cc_remove_all"]] = opt[["cc_remove_all"]]
if ("cc_rescore_after_merge" %in% names(opt)) param[["cc_rescore_after_merge"]] = opt[["cc_rescore_after_merge"]]

if (!is.null(opt[["vars_to_regress"]])) {
    param[["vars_to_regress"]] = opt[["vars_to_regress"]]
} else {
    param["vars_to_regress"] = list(NULL)
}

# Integration of multiple samples 
param[["integrate_samples"]] = list(method=opt[["combine_samples_method"]])
if (param[["integrate_samples"]] == "integrate") {
  param[["integrate_samples"]][["dimensions"]] = opt[["integrate_samples_dimensions"]]
  param[["integrate_samples"]][["reference"]] = opt[["integrate_samples_reference"]]
  param[["integrate_samples"]][["use_reciprocal_pca"]] = ifelse(opt[["integrate_samples_use_reciprocal_pca"]], TRUE, FALSE)
  param[["integrate_samples"]][["k_filter"]] = opt[["integrate_samples_k_filter"]]
  param[["integrate_samples"]][["k_weight"]] = opt[["integrate_samples_k_weight"]]
  param[["integrate_samples"]][["k_anchor"]] = opt[["integrate_samples_k_anchor"]]
  param[["integrate_samples"]][["k_score"]] = opt[["integrate_samples_k_score"]]
}

# Analysis
param[["pc_n"]] = opt[["pc_n"]]
param[["cluster_k"]] = opt[["cluster_k"]]
param[["umap_k"]] = opt[["umap_k"]]
param[["cluster_resolution"]] = opt[["cluster_resolution"]]
if (!is.null(opt[["cluster_resolution_test"]])) {
  param[["cluster_resolution_test"]] = opt[["cluster_resolution_test"]]
} else {
  param["cluster_resolution_test"] = list(NULL)
}

# Marker genes
param[["marker_padj"]] = opt[["marker_padj"]]
param[["marker_log2FC"]] = opt[["marker_log2FC"]]
param[["marker_pct"]] = opt[["marker_pct"]]

if (!is.null(opt[["latent_vars"]])) {
    param[["latent_vars"]] = opt[["latent_vars"]]
} else {
    param["latent_vars"] = list(NULL)
}

# Differentially expressed genes
if (!is.null(opt[["deg_contrasts"]])) {
    param[["deg_contrasts"]] = normalizePath(opt[["deg_contrasts"]])
} else {
    param["deg_contrasts"] = list(NULL)
}

# Functional enrichment
param[["enrichr_padj"]] = opt[["enrichr_padj"]]
param[["enrichr_dbs"]] = opt[["enrichr_dbs"]]
param[["enrichr_site"]] = opt[["enrichr_site"]]

# Colors
if (!opt[["col"]] %in% colours()) {
  stop("Only R colours allowed for --col!")
} else {
  param[["col"]] = opt[["col"]]
}

param[["col_palette_samples"]] = opt[["col_palette_samples"]]
param[["col_palette_clusters"]] = opt[["col_palette_clusters"]]

# Other parameters
param[["path_to_git"]] = normalizePath(opt[["path_to_git"]])
param[["debugging_mode"]] = opt[["debugging_mode"]]
param[["cores"]] = opt[["cores"]]

# Print parameter for log
print(param)

#
# Run run knitr
#
# Note: a small R script will be generated so that knitr can be run again if needed

# We already need to path_out here
if (dir.exists(param[["path_out"]])) {
  stop("The output directory already exsists. Please specify another path.")
}
dir.create(param[["path_out"]], recursive=TRUE, showWarnings=FALSE)
param[["path_out"]] = normalizePath(param[["path_out"]])

# This is the template for the script
script_template = '
#!/usr/bin/env Rscript

library(knitr)

rmarkdown::render(
    "{{rmd_dir}}/scrnaseq.Rmd", 
    output_format = "html_document", 
    output_dir = "{{output_dir}}", 
    intermediates_dir = "{{intermediates_dir}}", 
    output_file = "{{output_file}}", 
    knit_root_dir = "{{knit_root_dir}}", 
    quiet = {{quiet}}, 
    params = {{params_lst}})
'

# Now fill the script with content
script_content = knitr::knit_expand(text = script_template, 
                                    rmd_dir = param[["path_to_git"]], 
                                    output_dir = param[["path_out"]], 
                                    intermediates_dir = param[["path_out"]], 
                                    output_file = opt[["report_name"]], 
                                    knit_root_dir = param[["path_out"]], 
                                    quiet = !opt[["verbose"]], 
                                    params_lst = deparse(param, control="all"))
writeLines(script_content, file.path(param[["path_out"]], opt[["report_name"]] %>% paste0(".rerun.r")))

# Now run actual knitr process
rmarkdown::render(
    file.path(param[["path_to_git"]], "scrnaseq.Rmd"), 
    output_format = "html_document", 
    output_dir = param[["path_out"]], 
    intermediates_dir = param[["path_out"]], 
    output_file = opt[["report_name"]], 
    knit_root_dir = param[["path_out"]], 
    quiet = !opt[["verbose"]], 
    params = param)
