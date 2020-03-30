# TSSexploreR

## Installing TSSexploreR

**Create conda environment**
```
conda create -n tsrexplorer -y -c conda-forge -c bioconda \
r-tidyverse \
r-devtools \
r-ggseqlogo \
r-ggally \
bioconductor-genomicranges \
bioconductor-genomicfeatures \
bioconductor-biostrings \
bioconductor-rsamtools \
bioconductor-chipseeker \
bioconductor-edger \
bioconductor-clusterProfiler \
bioconductor-complexheatmap
```

**Install latest version of tsrexplorer**
```
devtools::install_github("rpolicastro/tsrexplorer")
```

## Using TSRexplorer

### Importing Data

TSSexplereR was designed to allow high flexibility and interoperability when import data.
TSSs and TSRs can be imported from common genomic formats such as bedgraphs and beds.
Furthermore, data can be imported directly from other R libraries such as TSRchitect and CAGEr.

[Data Import Vignette](./github/DATA_IMPORT.md)

### TSS Analysis


