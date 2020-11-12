# Introduction
**scrnaseq** is a bioinformatics analysis workflow for single-cell RNA-seq analysis. The workflow is based on Seurat, and contains additional visualisations, tables and documentation to better understand the analysis. The workflow supports RNA sequencing data from one or more samples processed with 10X Genomics and SmartSeq-2. 

# News
(2020-11-12)  
We just updated our master code branch, see [#59](https://github.com/ktrns/scrnaseq/pull/59). We updated our system to Seurat 4, python3, and all latest related R packages. We simplified the structure of our workflow, and focus on one type of normalisation per run. Cell cycle effects can be scored per sample, as well as after samples are combined. We added the visualisation of the highest expressed genes, and the number of cells per sample per cluster. We also added the export of average expression values alongside differentially expressed genes.  

We will next incorporate the analysis of differentially expressed genes between samples of origin, and will additionally allow for specific comparisons between individual clusters and samples. Once this is done, we will create our first code release! 

# Workflow summary
## Pre-Workflow: Demultiplexing with hashtag oligos
* Read input data
* Demutliplexing with hashtag oligos (HTOs)
   * Normalisation of HTO counts
   * Classification of cells based on normalised HTO data
* Visualisation of raw and normalised HTO data
* Remove cells classified as doublet or negative
* Preliminary visualisation of demultiplexed RNA data
   * Visualisation with UMAP
* Write out demultiplexed data
* Parameter table
* Software versions

## Workflow: Single-cell RNA-seq analysis 
* Dataset description
   * Project-specific parameters
* Read data
   * Read and print mapping statistics
   * Read gene annotation
   * Read scRNA-seq data
* Pre-processing
   * Quality control
   * Genes with highest expression
   * Filtering
   * Normalisation, scaling, variable genes, and cell cycle scoring
      * Variable genes
   * Combining multiple samples
   * Relative log expression
   * Dimensionality reduction
   * Dimensionality of the dataset
* Downstream analysis
   * Clustering
   * Visualisation with UMAP
   * Distribution of cells in clusters
   * Cell Cycle Effect
   * Cluster QC
   * Known marker genes
   * Differentially expressed genes, comparing one cluster against the rest (marker genes)
      * Table of top marker genes
      * Visualisation of top marker genes
      * Heatmaps of all differentially expressed genes
      * Functional enrichment analysis
* Further analysis with other tools
   * Export to Loupe Cell Browser
   * Export to the Cerebro Browser
* Output files
* Parameter table
* Software versions
* References

# Quick start
The workflow is inialised for test data in `test_datasets`. First, navigate to the respective test dataset folder(s), and download the test dataset(s) by running the `download.R` script(s). Once all test data is downloaded, you can knit the workflow to HTML. 

The repository provides several other useful test data that you can use to get to know the functionality of the workflow. To run the workflow for another than the initial dataset, you need to adapt the `project_parameters` code chunk and provide all relevant paths and parameters. 

# Documentation 

## Pre-Workflow: Demultiplexing with hashtag oligos

### Running the script
The Pre-Workflow can be run from outside the actual Rmarkdown script. When you are using the render function of the `rmarkdown` package you can run the script as follows:
```
rmarkdown::render(
    "scrnaseq_hto.Rmd",
    output_format = "html_document",
    output_dir = ".",
    output_file = "test",
    params = paramsList)
```
with paramsList set as:
```
paramsList = list()

paramsList$project = "HTO_testDataset"
paramsList$path_data = "test_datasets/10x_pbmc_hto_GSE108313/counts"
paramsList$path_out = "test_datasets/10x_pbmc_hto_GSE108313/demultiplexed"
paramsList$hto_names = setNames(c("htoA","htoB","htoC","htoD","htoE","htoF","htoG","htoH"), c("htoA","htoB","htoC","htoD","htoE","htoF","htoG","htoH"))
paramsList$mt = "^MT-"
paramsList$col = "palevioletred"
paramsList$sample_cells = NULL
```

### Arguments
#### `project`
ID of the project (Default: "HTO_testDataset").

#### `path_data`
Input directory where data are located (Default: "test_datasets/10x_pbmc_hto_GSE108313/counts").

#### `path_out`
Output directory where the results will be saved (Default: "test_datasets/10x_pbmc_hto_GSE108313/demultiplexed").

#### `hto_names`
HTOs have an ID that is included in the 'features.tsv' input file. We additionally ask for readable names that are used throughout the report. Names could look as follows, where `HTO1-3` are the IDs included in raw dataset: 
```param$hto.names = setNames(c("NameA", "NameB", "NameC"), c("HTO1", "HTO2", "HTO3"))```
(Default: `c("htoA", "htoB", "htoC", "htoD", "htoE", "htoF", "htoG", "htoH"), c("htoA", "htoB", "htoC", "htoD", "htoE", "htoF", "htoG", "htoH")`)

#### `mt`
Prefix of mitochondrial genes (Default: "^MT-").

#### `col`
Main colour(s) to use for plots (Defaults: "palevioletred").

#### `sample_cells`
Sample data to at most n cells (mainly for tests); set to NULL to deactivate (Default: NULL).

## Workflow: Single-cell RNA-seq analysis 

### Running the script
The main workflow is currently run from within Rstudio. Project-specific parameters are adapted in the `project_parameters` code chunk. 

# Credits
The [Seurat Vignettes](https://satijalab.org/seurat/vignettes.html) were initially used as templates for this workflow. 

The workflow was developed by [Katrin Sameith](https://github.com/ktrns) and [Andreas Petzold](https://github.com/andpet0101) at the [Dresden-concept Genome Center (Dresden, Germany)](https://genomecenter.tu-dresden.de/about-us). Through collaboration with the [Research Core Unit Genomics (Hannover, Germany)](https://www.mhh.de/genomics) the workflow has grown substantially and has been standardised. Many thanks to all who have contributed along the way, including (but not limited to): [Dimitra Alexopoulou](https://github.com/dimialex), [Mathias Lesche](https://github.com/mlesche), [Oliver Dittrich](https://github.com/Oliver-D-B), [Fabian Friedrich](https://github.com/Colorstorm), [Colin Davenport](https://github.com/colindaven), [Torsten Glomb](https://github.com/tglomb), and [Marius Rueve](https://github.com/mariusrueve).

# Citation
If you use the scrnaseq workflow to analyse your data, please cite it by mentioning the Dresden-concept Genome Center URL 'https://genomecenter.tu-dresden.de'. 
