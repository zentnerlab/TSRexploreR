# TSRexploreR

TSRexploreR allows for the processing and exploration of data generated from global transcription start site (TSS) mapping technologies such as CAGE and STRIPE-seq.
TSSs can be supplied as mapped reads in BAM/SAM format, or as aggregated TSSs in GRanges, bigwig, bedgraph, or CTSS format.
TSSs can then be normalized and clustered into transcription start regions (TSRs), or optionally TSRs generated from other software such as CAGEr can be imported.

One strength of TSRexploreR is the flexibility and ease in which data can be visualized.
Numerous plots can be generated for TSSs and TSRs such as density plots, heatmaps, sequence logos, gene tracks, correlation heatmaps, and much more.
For most plots we also provide the ability to filter, order, split, and group TSSs or TSRs by various calculated metrics such as IQR, quantiles, peak shape, etc.

In additional to differential TSSs and TSRs, we introduce a new algorithm for discovery of TSS cluster shifts, termed earth movers score (EMS).
Previous literature has established global TSS shifts in zebrafish development, as well as a phenomenon of mutations in certain transcription factors.
EMS is a sensible and intuitive approach to quantify the shift in position of discrete distributions,
and lies on a scale from -1 to 1, indicating upstream and downstream shifts respectively.

## Installing TSRexploreR

TSRexploreR can be installed wholly from github (with or without a conda environment), or entirely contained within a singularity container.

### Github

```
Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS"="true")
devtools::install_github("zentnerlab/TSRexploreR")
```

### Github and Conda

Conda can be optionally used before the above github installation step to download most of the package dependencies.
This also has the benefit of creating an R environment separate from your main software environment,
which helps avoid dependency version conflicts, and makes removing associated sofware easy if desired.
Installation instructions for miniconda3 can be found [here](https://docs.conda.io/projects/conda/en/latest/user-guide/install/).

First, create the conda environment.

```
conda create -n TSRexploreR -y -c conda-forge -c bioconda \
qpdf r-base>=4.0.0 r-tidyverse r-data.table r-devtools r-ggseqlogo r-cowplot r-rcpp \
r-assertthat r-testthat r-cairo r-ggrastr r-pkgdown \
bioconductor-apeglm bioconductor-genomicranges bioconductor-genomicfeatures \
bioconductor-genomicalignments bioconductor-biostrings bioconductor-rsamtools \
bioconductor-chipseeker bioconductor-edger bioconductor-deseq2 bioconductor-clusterProfiler \
bioconductor-complexheatmap bioconductor-cager bioconductor-tsrchitect bioconductor-gviz \
bioconductor-rtracklayer bioconductor-biocgenerics bioconductor-plyranges \
bioconductor-pcatools bioconductor-genomeinfodb bioconductor-bsgenome bioconductor-bioccheck \
bioconductor-bsgenome.scerevisiae.ucsc.saccer3 bioconductor-txdb.scerevisiae.ucsc.saccer3.sgdgene
```

After the environment has been created, it must be activated to be used.

```
conda activate TSRexploreR
```

TSRexploreR can now be installed into this environment,
where it now will be available when the environment is activated.

```
Sys.setenv("R_REMOTES_NO_ERRORS_FROM_WARNINGS"="true")
devtools::install_github("zentnerlab/TSRexploreR")
```

### Singularity.

The TSRexploreR singularity container is akin to a box that contains TSRexploreR and all required depndencies to run.
The advantage of containerized software is that they require no installation (beyond the main container software),
are self contained, and provide stable and reproducible results over time.
Singularity installation instructions can be found [here](https://sylabs.io/docs/).

After singularity is installed, the container can be pulled and run in a few commands.
Create a working directory with all required files, and run the following commands.

First, download the singularity container.

```
singularity pull --arch amd64 library://zentlab/default/tsrexplorer:main
```

You can then Start the R version from within the container that has TSRexploreR installed.

```
singularity exec -eCB `pwd` -H `pwd` tsrexplorer_main.sif R
```

## Using TSRexploreR

We prevoide numerous vignettes to showcase the various functionalities of TSRexplorer.

- [Standard Analysis](documentation/STANDARD_ANALYSIS.pdf) - Common workflow steps for TSS and TSR analysis and visualization.
- [BAM Processing](documentation/BAM_PROCESSING.pdf) - Importing and processing aligned reads from BAM/SAM format.
- [Differential TSSs/TSRs](documentation/DIFF_FEATURES.pdf) - Using edgeR or DESeq2 to find differential TSSs or TSRs between conditions.
- [TSS Cluster Shifts](documentation/FEATURE_SHIFT.pdf) - TSS cluster shifts between conditions using earth movers score (EMS).
- [Data Conditionals](documentation/DATA_CONDITIONING.pdf) - Flexible filtering, quantiling, ordering, and grouping TSSs or TSRs during plotting.
