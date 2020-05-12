# Introduction
**scrnaseq** is a bioinformatics analysis workflow for single-cell RNA-seq analysis using Seurat. The workflow currently supports RNA sequencing data derived for single samples processed with 10X Genomics and SmartSeq-2. 

# Workflow summary
## Pre-Workflow: Demultiplexing with hashtag oligos
* Dataset description 
  + Project-specific parameters  
* Read input data  
* Demutliplexing with hashtag oligos (HTOs)  
  + Normalisation of HTO counts  
  + Classification of cells based on normalised HTO data 
* Visualisation of raw and normalised HTO data 
* Remove cells classified as doublet or negative 
* Preliminary pre-processing of RNA data 
  + Visualisation of demultiplexed RNA data 
* Write out demultiplexed data 

## Workflow: Single-cell RNA-seq analysis 
* Dataset description
  + Project-specific parameters  
* Read input data  
  + Read and print mapping statistics  
  + Setup the Seurat object  
  + Read gene annotation  
* Pre-processing  
  + Quality control  
  + Normalisation, feature selection and scaling  
  + Dimensional reduction  
  + Dimensionality of the dataset  
* Downstream analysis  
  + Clustering  
  + Visualisation with UMAP  
  + Feature plots QC  
  + Feature plots for known marker genes  
  + Differentially expressed genes  
  + Visualisation of differentially expressed genes  
  + Functional enrichment analysis  
* Cell Cycle Effect  
* Loupe Cell Browser integration  
* Output files  

# Quick start
TODO: add this dataset into the repository and change text here
The Rmd script is initialised and tested for the following published 10X dataset:  
https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.0/pbmc_1k_v3

You need to download:  
* Feature / cell matrix (filtered)
* Summary CSV

... and modify the section "Project-specific paths" in the code chunk "projectParameters". 

# Documentation 
TODO

# Credits
The [Seurat Vignettes]() were initially used as templates for this workflow. 

The workflow was developed by [Katrin Sameith](https://github.com/ktrns) and [Andreas Petzold](https://github.com/andpet0101) at the [Dresden-concept Genome Center (Dresden, Germany)](https://genomecenter.tu-dresden.de/about-us). Through collaboration with the [Research Core Unit Genomics (Hannover, Germany)](https://www.mhh.de/genomics) the workflow has grown substantially and has been standardised. Many thanks to all who have contributed along the way, including (but not limited to): [Dimitra Alexopoulou](https://github.com/dimialex), [Mathias Lesche](https://github.com/mlesche), [Oliver Dittrich](https://github.com/Oliver-D-B), [Fabian Friedrich](https://github.com/Colorstorm), [Colin Davenport](https://github.com/colindaven), [Torsten Glomb](https://github.com/tglomb), and [Marius Rueve](https://github.com/mariusrueve).

# Citation
If you use the scrnaseq workflow to analyse your data, please cite it by mentioning the Dresden-concept Genome Center URL 'https://genomecenter.tu-dresden.de'. 
