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
    description="Demultiplexes a single cell dataset containing samples tagged by hashtags with the scrnaseq_hto workflow.",
    epilog="Upon successfull completion, the output directory will contain an HTML report, the counts per sample and a CSV file with the cell identities for import into the Loupe Cell Browser.")

parser$add_argument(c("--project-id"), action="store", help="Project ID", default="bfx", dest="project_id")
parser$add_argument(c("--path-data"), action="store", help="Path to a filtered counts directory created by the 10X cellranger pipeline", dest="path_data")
parser$add_argument(c("--stats"), action="store", help="Optional: path to a metrics_summary.csv file created by the 10X cellranger pipeline (default: %default)", default=NULL, dest="stats")
parser$add_argument(c("--path-out"), action="store", help="Path to an output directory (default: %default)", default="demultiplexed", dest="path_out")
parser$add_argument(c("--hto-features"), action="store", help="One or more HTO features for demultiplexing samples, separated by commas (e.g. HTO1,HTO2 ...). To assign a name to a sample, add a '=' to the HTO followed by the name (e.g. HTO1=SampleA,HTO2=SampleB). When the HTO features are not present, the script will exit with an error.", dest="hto_names")
parser$add_argument(c("--hto-features-regex"), action="store", help="Alternatively, specify a regular expression to identify HTO features. When the HTO features cannot be found, the script will exit without an error.", dest="hto_regex")
parser$add_argument(c("--norm"), action="store", help="Use LogNormalize or CLR for HTO normalisation (default: %default)", default="CLR", dest="norm")
parser$add_argument(c("--col"),  action="store", help="R colour for quantitative plots (default: %default)", default="palevioletred", dest="col")
parser$add_argument(c("--mt-names"),  action="store", help="Prefix of mitochondrial genes (default: %default)", default="^MT-", dest="mt_names")
parser$add_argument(c("--downsample-cells-n"),  action="store", type="integer", help="Downsample data to at most N cells (default: %default)", default=NULL, dest="downsample_cells_n")
parser$add_argument(c("--path-to-git"),  action="store", help="Path to the scrnaseq git repository (default: %default)", default=dirname(script_dir), dest="path_to_git")
parser$add_argument(c("--report-name"),  action="store", help="Name of the HTML report (default: %default)", default="scrnaseq_hto.html", dest="report_name")
parser$add_argument(c("--verbose"),  action="store_true", help="Be verbose (default: %default)", default=TRUE, dest="verbose")
opt = parser$parse_args()

#
# Check params
#
paramsList = list()

# Param project id
paramsList[["project_id"]] = as.character(opt[["project_id"]])

# Param path_data
if (!"path_data" %in% names(opt) || is.null(opt[["path_data"]])) stop("Need to provide the path to a filtered counts directory created by the 10X cellranger pipeline via --path-data!")
if (!file.exists(opt[["path_data"]])) stop("Path to a filtered counts directory specified via --path-data does not exists or is not a directory!")
paramsList[["path_data"]] = normalizePath(opt[["path_data"]])

# Param stats
if ("stats" %in% names(opt) && !is.null(opt[["stats"]])) {
    if (!file.exists(opt[["stats"]])) stop("Path to a metrics_summary.csv file specified via --stats does not exists or is not a file!")
    paramsList[["stats"]] = normalizePath(opt[["stats"]])
} else {
    paramsList["stats"] = list(NULL)
}

# Param path_out
if (!file.exists(opt[["path_out"]])) dir.create(opt[["path_out"]])
paramsList[["path_out"]] = normalizePath(opt[["path_out"]])

# Param hto-features
if ("hto_names" %in% names(opt) && !is.null(opt[["hto_names"]])) {
    hto_names_split = strsplit(trimws(unlist(strsplit(opt[["hto_names"]], ","))),"=")
    hto_names_n = sapply(hto_names_split, function(h) {return(h[1])})
    hto_names_s = sapply(hto_names_split, function(h) { if(length(h)==2) return(h[2]) else return(h[1])})
    
    paramsList[["hto_names"]] = setNames(hto_names_s, hto_names_n)
} else {
    paramsList["hto_names"] = list(NULL)
}

# Param hto-features-regex
if ("hto_regex" %in% names(opt) && !is.null(opt[["hto_regex"]])) {
    if (!is.null(paramsList[["hto_names"]])) stop("Cannot use --hto-features-regex together with --hto-features!")
    paramsList[["hto_regex"]] = trimws(opt[["hto_regex"]])
} else {
    paramsList["hto_regex"] = list(NULL)
}

if (is.null(paramsList[["hto_regex"]]) & is.null(paramsList[["hto_names"]])) stop("Need to specify HTO features for sample demultiplexing either via --hto-features or via --hto-features-regex!")

# Param norm
if ("norm" %in% names(opt)) {
    if (!opt[["norm"]] %in% c("CLR","LogNormalize")) stop("Only 'CLR' or 'LogNormalize' allowed for --norm!")
    paramsList[["norm"]] = opt[["norm"]]
}

# Param mt_names
if ("mt_names" %in% names(opt)) {
    paramsList[["mt_names"]] = opt[["mt_names"]]
}

# Param col
if ("col" %in% names(opt)) {
    if (!opt[["col"]] %in% colours()) stop("Only R colours allowed for --col!")
    paramsList[["col"]] = opt[["col"]]
}

# Param downsample_cells_n
if ("downsample_cells_n" %in% names(opt) && !is.null(opt[["downsample_cells_n"]])) {
    paramsList[["downsample_cells_n"]] = opt[["downsample_cells_n"]]
}

# Param path_to_git
if ("path_to_git" %in% names(opt)) {
    paramsList[["path_to_git"]] = opt[["path_to_git"]]
}

# Param debugging_mode: 'default_debugging' for RStudio, 'terminal_debugger' for debugging on a terminal without X11, 'print_traceback' for non-interactive sessions
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
    if (sum(grepl(paramsList[["hto_regex"]], features_tbl[ , 1, drop=TRUE])) == 0) {
        run_knitr = FALSE
        message("No HTOs identified by regular expression. Script will exit.")
    }
}

#
# Run run knitr
#

if (run_knitr) {

    # Note: a small R script will be generated so that knitr can be run again if needed

    # The name of the script
    script_name = paste0(opt[["report_name"]], ".rerun.r")

    # This is the template for the script
    script_template = '
#!/usr/bin/env Rscript

library(knitr)

rmarkdown::render(
    file.path("{{rmd_dir}}","scrnaseq_hto.Rmd"),
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
                                        rmd_dir=paramsList[["path_to_git"]],
                                        output_dir=paramsList[["path_out"]],
                                        intermediates_dir=paramsList[["path_out"]],
                                        output_file=opt[["report_name"]],
                                        knit_root_dir=paramsList[["path_out"]],
                                        quiet=!opt[["verbose"]],
                                        params_lst=deparse(paramsList, control="all"))
    
    
    writeLines(script_content, file.path(paramsList[["path_out"]], script_name))

    # Now run actual knitr process
    rmarkdown::render(
        file.path(paramsList[["path_to_git"]],"scrnaseq_hto.Rmd"),
        output_format="html_document",
        output_dir=paramsList[["path_out"]],
        intermediates_dir=paramsList[["path_out"]],
        output_file=opt[["report_name"]],
        knit_root_dir=paramsList[["path_out"]],
        quiet=!opt[["verbose"]],
        params=paramsList)

}