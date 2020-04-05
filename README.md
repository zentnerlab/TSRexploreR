# TSSexploreR

## Installing TSSexploreR

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
```

**Install latest version of tsrexplorer**
```
devtools::install_github("rpolicastro/tsrexplorer", ref = "clean")
```

## Using TSRexplorer

### Importing Data

TSSexplereR was designed to allow high flexibility and interoperability when import data.
TSSs and TSRs can be imported from common genomic formats such as bedgraphs and beds.
Furthermore, data can be imported directly from other R libraries such as TSRchitect and CAGEr.

[Data Import Vignette](./github/DATA_IMPORT.md)

### Standard Analysis

This vignette goes over the standard analysis of TSSs and TSRs.

[Standard Analysis Vignette](./github/STANDARD_ANALYSIS.md)

### Differential Features

Differential TSSs or TSRs can be calculated based on different conditions.

[Differential Features Vignette](./github/DIFF_FEATURES.md)
