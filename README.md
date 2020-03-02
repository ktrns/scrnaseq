# Introduction
Rmarkdown workflow for single-cell RNA-seq analysis using Seurat 

# Workflow summary 
* Current project  
  + Project-specific parameters  
* Read data  
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
The Rmd script is initialised and tested for the following published 10X dataset:  
https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.0/pbmc_1k_v3

You need to download:  
* Feature / cell matrix (filtered)
* Summary CSV

... and modify the section "Project-specific paths" in the code chunk "projectParameters". 

# Documentation 

# Credits
The workflow was originally written by [Katrin Sameith](https://github.com/ktrns) and her colleagues [Dimitra Alexopoulou](https://github.com/dimialex), [Andreas Petzold](https://github.com/andpet0101), and Mathias Lesche at the [Dresden-concept Genome Center, in Dresden, Germany](https://genomecenter.tu-dresden.de/about-us). 

The [Seurat Vignette](https://satijalab.org/seurat/v3.1/pbmc3k_tutorial.html) was initially used as a template for this workflow. 

# Citation
If you use the scrnaseq workflow to analyse your data, please cite it by mentioning the Dresden-concept Genome Center URL 'https://genomecenter.tu-dresden.de'. 
