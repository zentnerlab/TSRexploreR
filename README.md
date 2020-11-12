# TSRexploreR

## Installing TSRexploreR

**Create conda environment**
```
conda create -n tsrexplorer -y -c conda-forge -c bioconda \
r-tidyverse \
r-data.table \
r-devtools \
r-ggseqlogo \
r-ggally \
r-cowplot \
r-uwot \
r-rcpp \
bioconductor-genomicranges \
bioconductor-genomicfeatures \
bioconductor-biostrings \
bioconductor-rsamtools \
bioconductor-chipseeker \
bioconductor-edger \
bioconductor-clusterProfiler \
bioconductor-complexheatmap \
bioconductor-cager \
bioconductor-tsrchitect \
bioconductor-gviz \
bioconductor-rtracklayer
bioconductor-plyranges
bioconductor-deseq2
```

**Install latest version of tsrexplorer**
```
devtools::install_github("zentnerlab/TSRexploreR", ref = "main")
```

## Using TSRexploreR

### Importing Data

TSRexplereR was designed to allow high flexibility and interoperability when importing data. TSSs and TSRs can be imported from common genomic formats such as bedgraphs and beds. Furthermore, data can be imported directly from other R libraries such as TSRchitect and CAGEr.

[Data Import Vignette](./github/DATA_IMPORT.md)

### Standard Analysis

This vignette goes over the standard workflow for analysis of TSSs and TSRs.

[Standard Analysis Vignette](./github/STANDARD_ANALYSIS.md)

### Differential Features

This vignette describes detection of TSSs and TSRs that differ between samples.

[Differential Features Vignette](./github/DIFF_FEATURES.md)
