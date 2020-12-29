# Installation

## Method 1: via conda

```bash
#Download the Anaconda3 installer from https://www.anaconda.com/, then e.g.:
#
bash ~/Downloads/Anaconda3-2020.07-Linux-x86_64.sh
\rm  ~/Downloads/Anaconda3-2020.07-Linux-x86_64.sh

#In a new shell:
#
conda config --set auto_activate_base False
conda deactivate
conda --version
conda info --envs
```

Now create and activate the TSRexploreR environment:
```bash
conda create -n TSRexploreR -y -c conda-forge -c bioconda \
r-tidyverse \
r-data.table \
r-devtools \
r-ggseqlogo \
r-ggally \
r-cowplot \
r-rcpp \
r-assertthat \
r-testthat \
bioconductor-apeglm \
bioconductor-genomicranges \
bioconductor-genomicfeatures \
bioconductor-genomicalignments \
bioconductor-biostrings \
bioconductor-rsamtools \
bioconductor-chipseeker \
bioconductor-edger \
bioconductor-deseq2 \
bioconductor-clusterProfiler \
bioconductor-complexheatmap \
bioconductor-cager \
bioconductor-tsrchitect \
bioconductor-gviz \
bioconductor-rtracklayer \
bioconductor-biocgenerics \
bioconductor-plyranges \
bioconductor-pcatools \
bioconductor-GenomeInfoDb
bioconductor-bsgenome.scerevisiae.ucsc.saccer3 \
bioconductor-txdb.scerevisiae.ucsc.saccer3.sgdgene

conda activate TSRexploreR
```

Now in R
```bash
devtools::install_github("zentnerlab/TSRexploreR", ref="clean")
```
will load the TSRexploreR library.


## Method 2: create a local TSRexploreR librarry

In an R console:
```bash
install.packages("BiocManager")

BiocManager::install(c("annotate", "AnnotationDbi", "AnnotationFilter", "AnnotationHub", "askpass", "assertthat", "backports", "base64enc", "beanplot", "Biobase", "BiocFileCache", "BiocGenerics", "BiocParallel", "BiocVersion", "biomaRt", "Biostrings", "biovizBase", "bit", "bit64", "bitops", "blob", "boot", "BSgenome", "CAGEr", "Cairo", "caTools", "cellranger", "checkmate", "ChIPseeker", "circlize", "clue", "cluster", "colorspace", "compiler", "ComplexHeatmap", "cowplot", "crayon", "curl", "data.table", "DBI", "dbplyr", "DelayedArray", "DESeq2", "dichromat", "digest", "DO.db", "DOSE", "dplyr", "edgeR", "ellipsis", "enrichplot", "ensembldb", "farver", "fastmap", "fastmatch", "fgsea", "forcats", "foreign", "Formula", "formula.tools", "genefilter", "geneplotter", "generics", "GenomeInfoDb", "GenomeInfoDbData", "GenomicAlignments", "GenomicFeatures", "GenomicRanges", "GetoptLong", "ggforce", "ggplot2", "ggraph", "ggrepel", "ggseqlogo", "GlobalOptions", "glue", "GO.db", "GOSemSim", "gplots", "graphlayouts", "grid", "gridExtra", "gtable", "gtools", "Gviz", "Hmisc", "hms", "htmlTable", "htmltools", "htmlwidgets", "httpuv", "httr", "igraph", "interactiveDisplayBase", "IRanges", "jpeg", "KernSmooth", "knitr", "later", "lattice", "latticeExtra", "lazyeval", "lifecycle", "limma", "locfit", "magrittr", "MASS", "Matrix", "MatrixGenerics", "matrixStats", "memoise", "mgcv", "mime", "MultiAssayExperiment", "munsell", "nlme", "nnet", "openssl", "operator.tools", "parallel", "permute", "pillar", "pkgconfig", "plotrix", "plyr", "plyranges", "png", "polyclip", "prettyunits", "progress", "promises", "ProtGenerics", "purrr", "qvalue", "R6", "rappdirs", "RColorBrewer", "Rcpp", "RCurl", "readxl", "reshape", "reshape2", "rjson", "rlang", "rpart", "Rsamtools", "RSQLite", "rstudioapi", "rtracklayer", "rvcheck", "S4Vectors", "scales", "scatterpie", "shadowtext", "shape", "shiny", "som", "splines", "stats4", "stringdist", "stringi", "stringr", "SummarizedExperiment", "survival", "tibble", "tidygraph", "tidyr", "tidyselect", "tools", "TSRchitect", "tweenr", "TxDb.Hsapiens.UCSC.hg19.knownGene", "uwot", "VariantAnnotation", "vctrs", "vegan", "VGAM", "viridis", "viridisLite", "xfun", "XML", "xml2", "xtable", "XVector", "yaml", "zlibbioc"))
```
Note: You will have to run the BiocManager::install command several times until all dependencies are taken care of in proper order (the above listing is simply alphabetical).

Then create the TSRexploreR library in directory LOCALR in your home directory (or wherever else you desire; the command below are easily modified):

```bash
mkdir ~/LOCALR
git clone https://zentnerlab/TSRexploreR
R CMD build TSRexploreR
R CMD INSTALL -l ${HOME}/LOCALR/ TSRexploreR_0.1.tar.gz 
```

Now in an R console:
```
homeloc <- system2("echo","$HOME",stdout=TRUE)
libloc <- paste(libloc,"LOCALR",sep="/")
library("TSRexploreR",lib.loc=libloc)
```

will load the TSRexploreR library.
