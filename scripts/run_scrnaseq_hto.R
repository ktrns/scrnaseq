#!/usr/bin/env Rscript

arguments = commandArgs(trailingOnly = F)
# get path to script
script_index = grep("--file",arguments)
script_dir = dirname(sub("--file=","",arguments[script_index]))
script_name = basename(sub("--file=","",arguments[script_index]))

# Which packages are needed
library(knitr)
suppressPackageStartupMessages(library("optparse"))

# Prepare parser for arguments
option_list = list( 
    make_option(c("--project-id"), action="store", type="character", help="Project ID", default="bfx", dest="project_id"),
    make_option(c("--path-data"), action="store", type="character", help="Path to a filtered counts directory created by the 10X cellranger pipeline", dest="path_data"),
    make_option(c("--stats"), action="store", type="character", help="Optional: path to a metrics_summary.csv file created by the 10X cellranger pipeline (default: %default)", default=NULL, dest="stats"),
    make_option(c("--path-out"), action="store", type="character", help="Path to an output directory (default: %default)", default="demultiplexed", dest="path_out"),
    make_option(c("--hto-features"), action="store", type="character", help="One or more HTO features for demultiplexing samples, separated by commas (e.g. HTO1,HTO2 ...). To assign a name to a sample, add a '=' to the HTO followed by the name (e.g. HTO1=SampleA,HTO2=SampleB). When the HTO features are not present, the script will exit with an error.", dest="hto_names"),
    make_option(c("--hto-features-regex"), action="store", type="character", help="Alternatively, specify a regular expression to identify HTO features. When the HTO features cannot be found, the script will exit without an error.", dest="hto_features_regex"),
    make_option(c("--norm"), action="store", type="character", help="Use LogNormalize or CLR for HTO normalisation (default: %default)", default="LogNormalize", dest="norm"),
    make_option(c("--col"),  action="store", type="character", help="R colour for quantitative plots (default: %default)", default="palevioletred", dest="col"),
    make_option(c("--mt-names"),  action="store", type="character", help="Prefix of mitochondrial genes (default: %default)", default="^MT-", dest="mt_names"),
    make_option(c("--downsample-cells-n"),  action="store", type="integer", help="Downsample data to at most N cells (default: %default)", default=NULL, dest="downsample_cells_n"),
    make_option(c("--path-to-git"),  action="store", type="character", help="Path to the scrnaseq git repository (default: %default)", default=dirname(script_dir), dest="path_to_git"),
    make_option(c("--report-name"),  action="store", type="character", help="Name of the HTML report (default: %default)", default="scrnaseq_hto.html", dest="report_name"),
    make_option(c("--verbose"),  action="store_true", type="logical", help="Be verbose (default: %default)", default=TRUE, dest="verbose"))

opt = parse_args(OptionParser(usage="Usage: %prog [options]", 
    option_list=option_list, 
    add_help_option=TRUE,
    prog=script_name,
    description="Runs the scrnaseq_hto.Rmd workflow.\n\nAdditional description still needed.", 
    epilogue="Additional description still needed."))

#
# Check params
#
paramsList = list()

# project id
paramsList[["project_id"]] = as.character(opt[["project_id"]])

# path_data
if (!"path_data" %in% names(opt)) stop("Need to provide the path to a filtered counts directory created by the 10X cellranger pipeline via --path-data!")
if (!dir.exists(opt[["path_data"]])) stop("Path to a filtered counts directory specified via --path-data does not exists or is not a directory!")
paramsList[["path_data"]] = opt[["path_data"]]

# stats
if ("stats" %in% names(opt)) {
    if (!file.exists(opt[["stats"]])) stop("Path to a metrics_summary.csv file specified via --stats does not exists or is not a file!")
    paramsList[["stats"]] = opt[["stats"]]
} else {
    paramsList["stats"] = list(NULL)
}

# path_out
paramsList[["path_out"]] = opt[["path_out"]]

# hto_names
if ("hto_names" %in% names(opt)) {
    hto_names_split = strsplit(trimws(unlist(strsplit(opt[["hto_names"]], ","))),"=")
    hto_names_n = sapply(hto_names_split, function(h) {return(h[1])})
    hto_names_s = sapply(hto_names_split, function(h) { if(length(h)==2) return(h[2]) else return(h[1])})
    
    paramsList[["hto_names"]] = setNames(hto_names_s, hto_names_n)
} else {
    paramsList["hto_names"] = list(NULL)
}

# hto_regex
if ("hto_features_regex" %in% names(opt)) {
    if (!is.null(paramsList[["hto_names"]])) stop("Cannot use --hto-features-regex together with --hto-features!")
    paramsList[["hto_regex"]] = trimws(opt[["hto_features_regex"]])
} else {
    paramsList["hto_regex"] = list(NULL)
}

if (is.null(paramsList[["hto_regex"]]) & is.null(paramsList[["hto_names"]])) stop("Need to specify HTO features for sample demultiplexing either via --hto-features or via --hto-features-regex!")

# norm
if ("norm" %in% names(opt)) {
    if (!opt[["norm"]] %in% c("CLR","LogNormalize")) stop("Only 'CLR' or 'LogNormalize' allowed for --norm!")
    paramsList[["norm"]] = opt[["norm"]]
}

# mt_names
if ("mt_names" %in% names(opt)) {
    paramsList[["mt_names"]] = opt[["mt_names"]]
}

# col
if ("col" %in% names(opt)) {
    if (!opt[["col"]] %in% colours()) stop("Only R colours allowed for --col!")
    paramsList[["col"]] = opt[["col"]]
}

# downsample_cells_n
if ("downsample_cells_n" %in% names(opt)) {
    paramsList[["downsample_cells_n"]] = opt[["downsample_cells_n"]]
}

# path_to_git
if ("path_to_git" %in% names(opt)) {
    paramsList[["path_to_git"]] = opt[["path_to_git"]]
}

# debugging_mode: 'default_debugging' for RStudio, 'terminal_debugger' for debugging on a terminal without X11, 'print_traceback' for non-interactive sessions
paramsList[["debugging_mode"]] = 'print_traceback'

#
# Also check whether the dataset actually has the HTOs specified with --hto-features: if not, exit with error. 
# For HTOs specified with --hto-features-regex: if not, just exit (without error).
#
features_tbl = read.table(file.path(paramsList[["path_data"]] ,"features.tsv.gz"), sep="\t", header=FALSE)
run_knitr = TRUE

if (!is.null(paramsList[["hto_names"]])) {
    # At least one HTO feature name is not present in the dataset, stop with error
    if (any(!names(paramsList[["hto_names"]]) %in% features_tbl[, 1, drop=TRUE])) stop("Some HTO features are not present in the dataset!") else run_knitr = TRUE
} else if (!is.null(paramsList[["hto_regex"]])) {
    # Regular expression could not find HTO feature, just set hto_features_found = FALSE
    if (sum(grepl(paramsList[["hto_regex"]], features_tbl[ , 1, drop=TRUE])) == 0) run_knitr = FALSE else run_knitr = TRUE
}

#
# Run run knitr (if needed)
#

if (run_knitr) {
    rmarkdown::render(
        file.path(paramsList[["path_to_git"]],"scrnaseq_hto.Rmd"),
        output_format="html_document",
        output_dir=paramsList[["path_out"]],
        intermediates_dir=paramsList[["path_out"]],
        output_file=opt[["report_name"]],
        knit_root_dir=getwd(),
        quiet=!opt[["verbose"]],
        params=paramsList)
}

