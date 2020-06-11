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
* Preliminary visualisation of demultiplexed RNA data
  + Visualisation with UMAP
* Write out demultiplexed data

**Usage**  
rmarkdown::render with parameter list

Example:  
```
paramsList = list()

paramsList$project = "HTO_testDataset"
paramsList$path_data = "test_datasets/10x_pbmc_hto_GSE108313/counts"
paramsList$path_out = "test_datasets/10x_pbmc_hto_GSE108313/demultiplexed"
paramsList$hto_names = setNames(c("htoA","htoB","htoC","htoD","htoE","htoF","htoG","htoH"), c("htoA","htoB","htoC","htoD","htoE","htoF","htoG","htoH"))
paramsList$mt = "^MT-"
paramsList$col = "palevioletred"
paramsList$sample_cells = NULL

rmarkdown::render(
    "scrnaseq_hto.Rmd",
    output_format = "html_document",
    output_dir = ".",
    output_file = "test",
    params = paramsList)
```

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
You can have a quick start using a 10X Genomics PMBC dataset, downloaded from [here](https://support.10xgenomics.com/single-cell-gene-expression/datasets/3.0.0/pbmc_1k_v3) and included in <code>test_datasets</code>. The workflow <code>scrnaseq.Rmd</code> is already initialised for this dataset, so all you have to do is to knit to html. 

The repository provides several other test datasets that you can use for a quick start. All you have to do is to modify the code chunk <code>Project-specific parameters</code>. 

# Documentation 
TODO

# Credits
The [Seurat Vignettes](https://satijalab.org/seurat/vignettes.html) were initially used as templates for this workflow. 

The workflow was developed by [Katrin Sameith](https://github.com/ktrns) and [Andreas Petzold](https://github.com/andpet0101) at the [Dresden-concept Genome Center (Dresden, Germany)](https://genomecenter.tu-dresden.de/about-us). Through collaboration with the [Research Core Unit Genomics (Hannover, Germany)](https://www.mhh.de/genomics) the workflow has grown substantially and has been standardised. Many thanks to all who have contributed along the way, including (but not limited to): [Dimitra Alexopoulou](https://github.com/dimialex), [Mathias Lesche](https://github.com/mlesche), [Oliver Dittrich](https://github.com/Oliver-D-B), [Fabian Friedrich](https://github.com/Colorstorm), [Colin Davenport](https://github.com/colindaven), [Torsten Glomb](https://github.com/tglomb), and [Marius Rueve](https://github.com/mariusrueve).

# Citation
If you use the scrnaseq workflow to analyse your data, please cite it by mentioning the Dresden-concept Genome Center URL 'https://genomecenter.tu-dresden.de'. 
