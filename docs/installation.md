<img src="https://user-images.githubusercontent.com/9417927/155297862-5d3e1dc0-68f4-4bfe-933d-22f2936c9581.gif" alt="drawing" width="200"/>

To make our analyses as portable and reproducible as possible, we work with [singularity containers](https://sylabs.io/guides/3.0/user-guide/index.html). 

Our `scrnaseq` singularity container includes all software required to run the workflows. The container is generated once, and all team members can use it on their system. 

The container version is written into the workflow output. This way, we can later go back to the same container and rely on software versions to reproduce the analysis. 

As containers are big in space, we decided to share our [DcGC singularity recipes](https://github.com/dcgc-bfx/singularity-single-cell). You can download the recipe and build the container on your own.

