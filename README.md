# Table of contents
* [Introduction](#introduction) 
* [News](#news) 
* [Workflow summary](#workflow_summary)
* [Quick start](#quick_start)
* [Documentation](#documentation)
* [Credits](#credits)
* [Contributions and support](#contributions_and_support)
* [Citation](#citation)

# Introduction
<a name="introduction"/>

**scrnaseq** is a bioinformatics analysis workflow for single-cell RNA-seq analysis. The workflow is based on Seurat, and contains additional visualisations, tables and documentation to better understand the analysis. The workflow supports RNA sequencing data from one or more samples processed with 10X Genomics and SmartSeq-2. 

The workflow generates an extensive HTML report. Are you curious about what that looks like and if it would be useful for your own data? If so, you can  download the GitHub repository and open `scrnaseq.html`. This report has been generated for test data as mentioned [below](#quick_start). 

If you are a researcher, and you would like to start analysing your own data, the workflow can be your starting point. If you work in a *bioinformatics core facility* and frequently support other researchers with bioinformatics analyses, the workflow can be run in a standardised fashion both interactively in `Rstudio` and on command line. We typically first run the workflow with default parameters, and communicate with our collaborators. We then optimise the parameters in further rounds to improve the results. 


# Workflow summary
<a name="workflow_summary"/>


## Workflow: Single-cell RNA-seq analysis 
* Read data
   * Read and print mapping statistics
   * Read gene annotation
   * Read scRNA-seq data
* Pre-processing
   * Quality control
   * Genes with highest expression
   * Filtering
   * Quality control post filtering
   * Normalisation
   * Combining multiple samples
   * Relative log expression
   * Dimensionality reduction
   * Dimensionality of the dataset
* Clustering
* Visualisation with UMAP
* Distribution of samples in clusters
* Cell Cycle Effect
* Cluster QC
* Known marker genes
* Marker genes
    * Table to top marker genes
    * Visualisation of top marker genes
    * Heatmaps
    * Functional enrichment analysis
* Differentially expressed genes
* Output   
   * Export to Loupe Cell Browser
   * Export to cellxgene browser
   * Export to Cerebro browser
   * Other output files
* Parameter table
* Software versions
* Credits and References


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


# Quick start
<a name="quick_start"/>

The workflow is inialised for test data in `test_datasets`. First, navigate to the respective test dataset folder(s), and download the test dataset(s) by running the `download.R` script(s). Once all test data is downloaded, you can knit the workflow to HTML. 

The repository provides several other useful test data that you can use to get to know the functionality of the workflow. To run the workflow for another than the initial dataset, you need to adapt the `project_parameters` code chunk and provide all relevant paths and parameters. 


# Documentation 
<a name="documentation"/>

The scrnaseq workflow comes with a good amount of documentation, found in the `docs/` directory:
 
[Installation](docs/installation.md)   
[Running the workflow](docs/usage_workflow.md)   
[Running the pre-workflow](docs/usage_preworkflow.md)   

# Credits
<a name="credits"/>

The [Seurat Vignettes](https://satijalab.org/seurat/vignettes.html) were initially used as templates for this workflow. 

The workflow was developed by [Katrin Sameith](https://github.com/ktrns) and [Andreas Petzold](https://github.com/andpet0101) at the [Dresden-concept Genome Center (Dresden, Germany)](https://genomecenter.tu-dresden.de/about-us). Through collaboration with the [Research Core Unit Genomics (Hannover, Germany)](https://www.mhh.de/genomics) the workflow has grown substantially and has been standardised. Many thanks to all who have contributed along the way, including (but not limited to): [Oliver Dittrich](https://github.com/Oliver-D-B), [Torsten Glomb](https://github.com/tglomb), [Maike Kosanke](https://github.com/kosankem), [Dimitra Alexopoulou](https://github.com/dimialex), [Mathias Lesche](https://github.com/mlesche), [Fabian Friedrich](https://github.com/Colorstorm), [Colin Davenport](https://github.com/colindaven), and [Marius Rueve](https://github.com/mariusrueve).

# Contributions and Support
<a name="contributions_and_support"/>

If you would like to contribute to this workflow, please first create your own fork of the GitHub repository. You can then work on your master branch, or create feature branches for developement that you later merge into your master branch. Once your code is finalised and working, you can create a pull request. 

# Citation
<a name="citation"/>

If you use the scrnaseq workflow to analyse your data, please cite it by mentioning the Dresden-concept Genome Center URL "https://genomecenter.tu-dresden.de". 
