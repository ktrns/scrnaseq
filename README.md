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

The workflow generates an extensive HTML report. Are you curious about what that looks like and if it would be useful for your own data? If so, you can  download the GitHub repository and open `scrnaseq.html`. This report generated for test data as mentioned [below](#quick_start). 

# News
<a name="news"/>

(2021-12-02)  
We updated our master code branch, and merged two recent pull requests (see [#92](https://github.com/ktrns/scrnaseq/pull/92) and [#104](https://github.com/ktrns/scrnaseq/pull/104)). 

The most important changes are: 

* Introduced the possiblity to run the workflow directly from command-line ([#86](https://github.com/ktrns/scrnaseq/issues/86))
* Added UMAPs that display samples separately ([#98](https://github.com/ktrns/scrnaseq/issues/98))
* Added possibility to adapt dot size in UMAPs ([#99](https://github.com/ktrns/scrnaseq/issues/99))
* Added UMAPs that show the effect of different cluster resolution values ([#97](https://github.com/ktrns/scrnaseq/issues/97))
* Updated tables that show QC filtering results ([#95](https://github.com/ktrns/scrnaseq/issues/95))
* Added extensive explanations into the report ([#43](https://github.com/ktrns/scrnaseq/issues/42))

Stay tuned, we will soon create our first code release! 

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
* Distribution of cells in clusters
* Cell Cycle Effect
* Cluster QC
* Known marker genes
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
[Output and how to interpret the results](docs/output.md)

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
