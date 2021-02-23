# Table of contents
* [Introduction](#introduction) 
* [News](#news) 
* [Workflow summary](#workflow_summary)
* [Quick start](#quick_start)
* [Documentation](#documentation)
  * [Demultiplexing hashtag oligos](#documentation_hto)
    * [Running the script](#documentation_hto_script)
    * [Arguments](#documentation_hto_arguments)
  * [Single-cell RNA-seq analysis](#documentation_scrnaseq)
    * [Running the script](#documentation_scrnaseq_script)
    * [Arguments](#documentation_scrnaseq_arguments)
      * [General](#documentation_scrnaseq_arguments_general)
      * [Input](#documentation_scrnaseq_arguments_input)
      * [Filtering](#documentation_scrnaseq_arguments_filtering)
      * [Normalisation and integration, clustering and dimensionality reduction](#documentation_scrnaseq_arguments_normalisation)
      * [Marker genes and genes with differential expression](#documentation_scrnaseq_arguments_marker_degs)
* [Credits](#credits)
* [Citation](#citation)

# Introduction
<a name="introduction"/>

**scrnaseq** is a bioinformatics analysis workflow for single-cell RNA-seq analysis. The workflow is based on Seurat, and contains additional visualisations, tables and documentation to better understand the analysis. The workflow supports RNA sequencing data from one or more samples processed with 10X Genomics and SmartSeq-2. 


# News
<a name="news"/>

(2021-02-23)  
We again updated our master code branch (see [#67](https://github.com/ktrns/scrnaseq/pull/67) and [#68](https://github.com/ktrns/scrnaseq/pull/68)).

The biggest changes are:

* Introduced possibility to identify differentially expressed genes between any cell sets of interest, followed by functional enrichment analysis ([#45](https://github.com/ktrns/scrnaseq/issues/45))
* Updated README and improved documentation ([#60](https://github.com/ktrns/scrnaseq/issues/60))
* Figures are exported to png and pdf ([#62](https://github.com/ktrns/scrnaseq/issues/62))
* Fixed bugs ([#63](https://github.com/ktrns/scrnaseq/issues/63) and [#65](https://github.com/ktrns/scrnaseq/issues/65))
* Updated pre-workflow

Stay tuned, we will soon create our first code release! 

# Workflow summary
<a name="workflow_summary"/>

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
* Further analysis with other tools
   * Export to Loupe Cell Browser
   * Export to the Cerebro Browser
* Output files
* Parameter table
* Software versions
* References


# Quick start
<a name="quick_start"/>

The workflow is inialised for test data in `test_datasets`. First, navigate to the respective test dataset folder(s), and download the test dataset(s) by running the `download.R` script(s). Once all test data is downloaded, you can knit the workflow to HTML. 

The repository provides several other useful test data that you can use to get to know the functionality of the workflow. To run the workflow for another than the initial dataset, you need to adapt the `project_parameters` code chunk and provide all relevant paths and parameters. 


# Documentation 
<a name="documentation"/>

## Pre-Workflow: Demultiplexing with hashtag oligos
<a name="documentation_hto"/>

### Running the script
<a name="documentation_hto_script"/>

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
<a name="documentation_hto_arguments"/>

#### `project`
ID of the project (Default: "HTO_testDataset")

#### `path_data`
Input directory where data are located (Default: "test_datasets/10x_pbmc_hto_GSE108313/counts")

#### `path_out`
Output directory where the results will be saved (Default: "test_datasets/10x_pbmc_hto_GSE108313/demultiplexed")

#### `hto_names`
HTOs have an ID that is included in the "features.tsv" input file. We additionally ask for readable names that are used throughout the report. Names could look as follows, where `HTO1-3` are the IDs included in raw dataset: 
```param$hto.names = setNames(c("NameA", "NameB", "NameC"), c("HTO1", "HTO2", "HTO3"))```
(Default: `c("htoA", "htoB", "htoC", "htoD", "htoE", "htoF", "htoG", "htoH"), c("htoA", "htoB", "htoC", "htoD", "htoE", "htoF", "htoG", "htoH")`)

#### `mt`
Prefix of mitochondrial genes (Default: "^MT-")

#### `col`
Main colour(s) to use for plots (Defaults: "palevioletred")

#### `sample_cells`
Sample data to at most `n` cells (mainly for tests); set to NULL to deactivate (Default: NULL)


## Workflow: Single-cell RNA-seq analysis
<a name="documentation_scrnaseq"/>

### Running the script
<a name="documentation_scrnaseq_script"/>

The main workflow is currently run from within Rstudio. Project-specific parameters are adapted in the `project_parameters` code chunk. 

### Arguments
<a name="documentation_scrnaseq_arguments"/>

#### General
<a name="documentation_scrnaseq_arguments_general"/>

##### __`project_id`__
ID or name of the project

##### __`col`__

Main colour used for continuous data

##### __`col_palette_samples`__

Colour palette used for samples. Ideally, this should be one of ggsci palettes (https://cran.r-project.org/web/packages/ggsci/vignettes/ggsci.html), but it can also be another palette, as long as it is a function that provides a colour for a sample number. Make sure that the palette provides enough colours. 

##### __`col_palette_clusters`__

Colour palette used for clusters. Ideally, this should be one of ggsci palettes (https://cran.r-project.org/web/packages/ggsci/vignettes/ggsci.html), but it can also be another palette, as long as it is a function that provides a colour for a cluster number. Make sure that the palette provides enough colours.

##### __`path_to_git`__

Path to the local git repository

##### __`default_debugging`__

Debugging mode: "default_debugging" for default, "terminal_debugger" for debugging without X11, "print_traceback" for non-interactive sessions 


#### Input
<a name="documentation_scrnaseq_arguments_input"/>

##### __`path_data`__

R `data.frame` with information on the counts datasets. Each row corresponds to one dataset, and each column to one of the following information:

* `name`: Uunique identifier for a dataset
* `type`: `10x` for sparse matrix-like datasets as generated by the 10x pipeline Cell Ranger, or `smartseq2` for non-sparse counts tables
* `path`: Path to the respective counts directory or file
* `stats`: Path to a mapping statistics file, e.g. "metrics_summary.csv" generated by Cell Ranger, but can be any table or `NA` if not available

For dataset type `10x`, there should be a path to a directory with the following files:

* _features.tsv.gz_: Tab-separated file with feature id, feature name and type
* _barcodes.tsv.gz_: File with the cell names
* _matrix.mtx.gz_: Sparse matrix file
* _metadata.tsv.gz_: Tab-separated cell metadata file with the cell barcode in the first column (optional)

For dataset type `smartseq2`, there should be:

* Gzipped tab-separated counts file with the first column containing the gene id and the remaining columns containing the counts for the cells
* Additional cell metadata in a tab-separated _metadata.tsv.gz_ file that contains the cell barcode in the first column and that must be in the same directory (optional)

For dataset type `smartseq2`, cell names can contain plate information in the following format `<sample>_<plate_number>_<row><column>` where `sample` is a string, `plate_number` a two digit number, `row` a capital letter and `column` a number. In this case, plate information will be parsed and cells will be grouped into multiple samples based on the `sample` string in their cell name.

When providing multiple datasets, cell names are expected to be unique. If not, they will be made unique by adding a numerical suffix. Furthermore, it is highly recommended to use Ensembl ids as gene ids. If not, then it should be made sure that the gene ids are unique and only contain only letters, digits, hyphens and dots.

##### __`downsample_cells_n`__

Downsample data  to `n` cells for each sample (mainly for tests); set to NULL to deactivate

##### __`path_out`__

Path to an output directory (will be created if it does not exist)

##### __`file_known_markers`__

Path to an Excel file with marker genes based on literature. There should be one list per column with the first row as header and the gene ids listed below. The expression of known marker genes will be visualised in the report. 

##### __`mart_dataset`__

When using Ensembl ids as gene ids, the name of the Ensembl dataset for further information on the genes . For human, this will be `hsapiens_gene_ensembl`, for mouse `mmusculus_gene_ensembl` and for zebrafish `drerio_gene_ensembl`. For other species, there will be an analogous name.

##### __`annot_version`__

When using Ensembl ids as gene ids, the Ensembl version to use. This should match the version that was used for obtaining the raw counts.

##### __`annot_main`__

When using Ensembl ids as gene ids, the argument is an R named list that specifies which Ensembl attributes will be used as gene id, as gene symbol and as entrez accession, respectively.

##### __`file_annot`__

Instead of fetching the gene annotation from Ensembl servers, read it from an external file. Note that by default the gene annotation will be fetched only once and is then saved in a separate file. Use this argument only if there is no access to Ensembl at all.

##### __`mart_attributes`__

Gene annotation attributes to include when fetching gene annotation from Ensembl

##### __`biomart_mirror`__

`biomaRt` mirror to query Ensembl (`www`, `useast`, `uswest` or `asia`; if `NULL`, the mirror will be automatically selected)

#### Filtering
<a name="documentation_scrnaseq_arguments_filtering"/>

##### __`mt`__

Prefix of mitochondrial genes, usually `^MT-` for human, `^Mt-` for mouse and `^mt-` for zebrafish. The caret indicates the prefix.

##### __`cell_filter`__

R named list with filters for cells. Filters can be specified based on numerical and categorical column names in the cell metadata. Numerical filters must have two values indicating the minimum and maximum. If there is no minimum or maximum, the respective value can be set to `NA`. Categorial filters consist of a character vector of allowed values. The following cell metadata columns are computed by default:

* `nCount_RNA`: Total number of counts calculated based on the gene expression data
* `nFeature_RNA`: Total number of detected genes calculated based on the gene expression data
* `percent_mt`: Percent mitochondrial content calculated based on the gene expression data
* `nCount_ERCC`: Total number of counts calculated based on the ERCC spike-in data (if available)
* `nFeature_ERCC`: Total number of ERCCs calculated based on the ERCC spike-in data (if available)
* `percent_ercc`: Percent ERCC content calculated based on the ERCC and gene expression data (if available)
* `orig.ident`: Name of the sample

By default, top-level filters will be applied to all samples. To apply sample-specific filters, the R named list can contain sublists with the sample names and the respective sample-specific filters, which will overwrite the top-level filters. The names of the samples correspond to the dataset names or alternatively, if a dataset contains multiple samples (SmartSeq-2) to the dataset and sample name separated by a dot.

Here are some examples. The filters will keep:

* `list(nFeature_RNA=c(200, NA), percent_mt=c(NA, 20))`: cells with at least 200 genes expressed and at most a mitochondrial content of 20%
* `list(nFeature_RNA=c(200, NA), orig.ident=c("sampleA", "sampleB""))`: cells of sample A and B with at least 200 genes expressed
* `list(sampleA=list(nFeature_RNA=c(200, NA)), sampleB=list(nFeature_RNA=c(2000, NA)))`: cells of sample A with at least 200 genes and cells of sample B with at least 2000 genes
* `list(nFeature_RNA=c(200, NA), sampleB=list(nFeature_RNA=c(2000, NA)))`: cells of sample B with at least 2000 genes and cells of all other samples with at least 200 genes

Filtering can also be done based on user-provided cell metadata. This can be useful for example for subclustering. In a first pass, the big dataset is clustered into general clusters. In a second pass, the cluster information is provided as additional metadata to analyse only cells of a specific cluster.

##### __`feature_filter`__

R named list to filters for genes in the gene expression assay. The filter `min_counts` specifies the minimum number of counts a gene must have to be considered expressed. The filter `min_cells` specifies the minimum number of cells that must express a gene. Similarly, to the `cell_filter`, there can be sample-specific filtering, i.e. the R named list can contain sublists with the sample names and the respective sample-specific filters.

Here are some examples. The filters will keep:

* `list(min_counts=1, min_cells=3)`: genes with at least one count expressed in at least three cells
* `list(sampleA=list(min_counts=1, min_cells=3), sampleB=list(min_counts=1, min_cells=10))`: genes with at least one count and expressed in at least three cells for sampleA and at least 10 cells for sampleB
* `list(sampleB=list(min_cells=10), min_counts=1, min_cells=3)`: genes with at least one count and expressed in at least 10 cells for sampleB and at least three cells for all other samples

##### __`samples_to_drop`__

Names of samples to drop after initial QC. This can be useful to remove technical controls (e.g. NC/negative control). The names of the samples correspond to the dataset names or, if a dataset contains multiple samples (SmartSeq-2) to the dataset and sample name separated by a dot.

##### __`samples_min_cells`__

Minimum number of cells a sample must have. This can be useful to remove samples which are part of a dataset but to few cells for analysis.

#### Normalisation and integration, clustering and dimensionality reduction
<a name="documentation_scrnaseq_arguments_normalisation"/>

##### __`norm`__

Normalisation method to use for the analysis. Can be:

* `RNA`: Counts are scaled to 10,000 and then natural log-transformed (default method by Seurat). Works well for all kinds of single cell data.
* `SCT`: SCTransform (Hafemeister 2019). Uses regularized negative binomial regression to account for varying gene expression levels. Expects molecule counts.

The `SCT` method is recommended for analyses where all datasets are derived from UMI-based technologies that generate molecule counts. If at least one dataset is derived from read-based methods such as SmartSeq-2, it is better to use the `RNA` method. Also, at the moment, it is not possible to use different, sample-specific normalisation methods.

##### __`cc_remove`__
Whether to remove the effect of cell cycle genes prior to dimensionality reduction, visualisation and clustering. Gene expression in cycling cells can be dominated by a set of cell cycle-related genes which will mask the true relationship between cells. By default, a cell cycle score will be calculated for each cell based on a defined list of known cell cycle marker genes and the cell will be classified into G1, G2M or S phase. If `TRUE`, this information will then be used to regress out (that is remove) potential cell cycle-related contributions for each gene prior to further analysis. Note that this will likely work only for human and mouse. 

##### __`cc_remove_all`__

If `TRUE`, cell cycle-related contributions for a gene will be assessed on the assigned cell cycle phases of the cells and all differences that likely originate from being in different cell cycle phases will be removed. However, when proliferation is an import aspect of the dataset, for example when analysing a developmental trajectory, this removal may mask important structures in the dataset. In this case, a more sensitive approach is to remove only the difference between proliferation cells in G2M and S phase (set to `FALSE`). Please read https://satijalab.org/seurat/v3.1/cell_cycle_vignette.html for more explanations.

##### __`cc_rescore_after_merge`__

When there are multiple datasets, should cell cycle-effects be scored on the individual datasets (`FALSE`) or on the combined dataset (`TRUE`). In general, cell cycle-scores are more reliable when calculated with more cells. However, there have been cases where calculation based on individual datasets worked better.

##### __`vars_to_regress`__

Additional variables that need to be regressed out prior to dimensionality reduction, visualisation and clustering. This can be useful when different individuals introduce an unwanted batch effect. Continuous variables such as the sequencing depth can be used as well. Variables must be names of the cell metadata.  

##### __`integrate_samples`__

When there are multiple datasets, how to integrate and combine them. R named list with:

* `method`:
  * `single`: Default, if there is only one dataset after filtering and no integration is needed 
  * `merge`: Datasets are merged (concatenated), and then processed as a whole. Since this is the least "invasive" method, it is recommended to start with this method and then to check whether there are unwanted sample-specific effects.
  * `standard`: This is the standard method for _integrating_ samples as implemented by Seurat. Anchors are computed for all pairs of datasets based on a set of variable (hence informative) genes. These anchors are then used to sensibly group cells of different samples and to harmonize the gene expression between the samples based on the correct grouping.
  * `reference`: One dataset is used as reference and anchors are computed for all other datasets. This method is less accurate but computationally faster (not yet implemented).
  * `reciprocal`: Anchors are computed in PCA space instead based on the actual data. This method is even less accurate but computationally faster especially for very big datasets (not yet implemented).
* `reference_dataset`: When using the method `reference`, which dataset is the reference? Can be numeric or name of the dataset.
* `dimensions`: Number of dimensions to consider for integration. More dimensions will cause integration to be more pronounced on details.

Note that the integration methods `standard`, `reference` and `reciprocal` are only done on a set of variable (hence informative) genes, which will later be used for dimensionality reduction, visualisation and clustering. Finding marker genes or identifying genes with differential expression will still be done on the full non-integrated dataset but will include the sample information as additonal factor to consider. The integrated data will be stored in a separate assay. 

##### __`pc_n`__

The number of principle components (PCs) used for dimensionaly reduction. After normalisation, scaling and integration (if requested), a principle component analysis (PCA) is done. To denoise the dataset and to reduce the complexity of the analysis, only the first `n` components are kept for further analysis such as clustering and visualisation. This parameter can be adjusted based on the Elbow plot in the report, which shows the percentage of variation explained by the PCs. A good number of PCs can be estimated based on where the plot reaches a plateau phase such that further PCs contribute only marginally to the variation. 

##### __`cluster_resolution`__

The resolution of the clustering algorithm. This controls the granularity of the clustering, i.e. lower values will lead to fewer clusters and higher values will lead to more clusters. Usually between 0.2 and 2 (but can be higher).

#### Marker genes and genes with differential expression
<a name="documentation_scrnaseq_arguments_degs"/>

##### __`marker_padj`__

Adjusted p-value for defining a marker gene

##### __`marker_log2FC`__

Minimum absolute log2 fold change for defining a marker gene

##### __`marker_pct`__

Minimum fraction of cells expressing a marker gene

##### __`latent_vars`__

Confounding variables to adjust for in finding markers and differentially expressed genes

##### __`deg_contrasts`__

Contrasts for testing differential expression. Can be a table in form of an R `data.frame` or an Excel file. The table needs to have the following columns:

* `condition_column`: The categorial column in the cell metadata to test. Special columns are `orig.ident` for sample, `seurat_clusters` for cluster and `Phase` cell cycle phase.
* `condition_group1`: The condition level(s) in group 1 (nominator). The group can contain multiple levels concatenated by the `+` character. An empty string means that all levels not in group 2 will be used.
* `condition_group2`: The condition level(s) in group 2 (denominator). The group can contain multiple levels concatenated by the `+` character. An empty string means that all levels not in group 1 will be used.

These columns are optional:

* `subset_column`: Prior to testing, subset cells based on this cell metadata column. Special columns are `orig.ident` for sample, `seurat_clusters` for cluster and `Phase` cell cycle phase. Must be used together with `subset_group` (default: none).
* `subset_group`: Levels to subset before testing. Multiple levels to be analysed can be given separate by semicolons. An empty string means that all levels will be analysed. Levels can be concatenated with the `+` character. Must be used together with `subset_column` (default: none).
* `assay`: Assay or dimensionality reduction to test on (default: `RNA`).
* `slot`: When an assay is selected, which slot to use:
  * raw `counts`: Default when `test` is one of `negbinom`, `poisson` or `DESeq2`
  * normalised `data`: Default for all other tests 
  * normalised and scaled `scale.data`: Should not be used
* `padj`: Maximum adjusted p-value (default: 0.05)
* `log2FC`: Minimum absolute log2 fold change (default: 0)
* `min_pct`: Minimum percentage of cells that need to express a gene in each group (default: 0.1)
* `test`: Test to use. Can be:
  * `wilcox`: Wilcoxon Rank Sum test (default)
  * `bimod`: Likelihood-ratio test for single cell gene expression, (McDavid et al., Bioinformatics, 2013)
  * `roc`: Identifies markers of gene expression using ROC analysis.
  * `t`: t-test
  * `negbinom`: tests based on a negative binomial generalized linear model
  * `poisson`: tests based on a poisson generalized linear model (only UMI datasets)
  * `LR`: tests based on a logistic regression model
  * `MAST`:  uses a hurdle model tailored to scRNA-seq data implemented in the the MAST package
  * `DESeq2`: tests using the DESeq2 package
  * For more information on the tests, please see the documentation on the `FindMarkers` function of the Seurat R package 
* `latent_vars`: Additional variables to account for when testing (e.g. batches). Can be one or more cell metadata columns, can contain multiple column names concatenated by semicolons.

Here are some examples for a better understanding:

* Compare cluster 1 versus cluster 2:
  * `condition_column`: `seurat_clusters`
  * `condition_group1`: `1`
  * `condition_group2`: `2`
* Compare cluster 1 versus cluster 2 and cluster 3 combined :
  * `condition_column`: `seurat_clusters`
  * `condition_group1`: `1`
  * `condition_group2`: `2+3`
* Compare cluster 1 versus rest:
  * `condition_column`: `seurat_clusters`
  * `condition_group1`: `1`
  * `condition_group2`: ` `
* Compare all G1 cells versus all S cells:
  * `condition_column`: `Phase`
  * `condition_group1`: `G1`
  * `condition_group2`: `S`
* Compare G1 cells versus S cells for cluster 1:
  * `condition_column`: `Phase`
  * `condition_group1`: `G1`
  * `condition_group2`: `S`
  * `subset_column`: `seurat_clusters`
  * `subset_group`: `1`
* Compare G1 cells versus S cells, separately for cluster 1, 2, and 3:
  * `condition_column`: `Phase`
  * `condition_group1`: `G1`
  * `condition_group2`: `S`
  * `subset_column`: `seurat_clusters`
  * `subset_group`: `1;2;3`
  * this will result in 3 tests
* Compare G1 cells versus S cells, combined for cluster 1+2+3:
  * `condition_column`: `Phase`
  * `condition_group1`: `G1`
  * `condition_group2`: `S`
  * `subset_column`: `seurat_clusters`
  * `subset_group`: `1+2+3`
  * this will result in 1 test
* Compare G1 cells versus S cells, combined for cluster 1+2+3, and separately for cluster 4:
  * `condition_column`: `Phase`
  * `condition_group1`: `G1`
  * `condition_group2`: `S`
  * `subset_column`: `seurat_clusters`
  * `subset_group`: `1+2+3;4`
  * this will result in 2 tests

##### __`enrichr_padj`__

P-value threshold for functional enrichment tests by Enrichr

##### __`enrichr_dbs`__

Enrichr databases for functional enrichment tests. Please see the corresponding table in the HTML output.


# Credits
<a name="credits"/>

The [Seurat Vignettes](https://satijalab.org/seurat/vignettes.html) were initially used as templates for this workflow. 

The workflow was developed by [Katrin Sameith](https://github.com/ktrns) and [Andreas Petzold](https://github.com/andpet0101) at the [Dresden-concept Genome Center (Dresden, Germany)](https://genomecenter.tu-dresden.de/about-us). Through collaboration with the [Research Core Unit Genomics (Hannover, Germany)](https://www.mhh.de/genomics) the workflow has grown substantially and has been standardised. Many thanks to all who have contributed along the way, including (but not limited to): [Dimitra Alexopoulou](https://github.com/dimialex), [Mathias Lesche](https://github.com/mlesche), [Oliver Dittrich](https://github.com/Oliver-D-B), [Fabian Friedrich](https://github.com/Colorstorm), [Colin Davenport](https://github.com/colindaven), [Torsten Glomb](https://github.com/tglomb), and [Marius Rueve](https://github.com/mariusrueve).

# Citation
<a name="citation"/>

If you use the scrnaseq workflow to analyse your data, please cite it by mentioning the Dresden-concept Genome Center URL "https://genomecenter.tu-dresden.de". 
