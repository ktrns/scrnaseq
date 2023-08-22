<img src="https://user-images.githubusercontent.com/9417927/155297862-5d3e1dc0-68f4-4bfe-933d-22f2936c9581.gif" alt="drawing" width="200"/>

To make our analyses as portable and reproducible as possible, we work with [singularity containers](https://sylabs.io/guides/3.0/user-guide/index.html). 

Our `scrnaseq` singularity container includes all software required to run the workflows. The container is generated once, and all team members can use it on their system. 

The container version is written into the workflow output. This way, we can later go back to the same container and rely on software versions to reproduce the analysis. 

You can download the container as follows: 

1. Install Singularity
2. Pull container: `singularity pull oras://gcr.hrz.tu-chemnitz.de/dcgc-bfx/singularity/singularity-single-cell:v1.4.0`

You can then run Rstudio-Server:

3. Generate `<my_dir>/rstudio-server/run` and `<my_dir>/rstudio-server/lib` directories
4. `singularity run -B <my_dir>/rstudio-server/run:/var/run/rstudio-server -B <my_dir>/rstudio-server/lib:/var/lib/rstudio-server --writable-tmpfs --app rserver singularity-single-cell_v1.4.0.sif 8789`
5. Open localhost:8789 in your web browser
