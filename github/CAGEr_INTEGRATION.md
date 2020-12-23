# CAGEr Integration

The [Bioconductor](https://bioconductor.org/) library [CAGEr](https://bioconductor.org/packages/release/bioc/html/CAGEr.html) has been extensively developed, and contains various normalization and clustering scheme options that may be desired.
In order to facilitate interoperability, convenience functions have been
provided with [TSRexploreR](https://github.com/zentnerlab/TSRexploreR) to import various bits of data.

## Preparing Data

This vignette will follow the steps in the [CAGEr](https://bioconductor.org/packages/release/bioc/html/CAGEr.html) vignette using the supplied *D. rerio* CAGE dervived TSSs and TSRs.

### Preparing CAGEexp Object

The first step is to create a [CAGEr](https://bioconductor.org/packages/release/bioc/html/CAGEr.html) object and import the data supplied by their package.

```
library("CAGEr")

inputFiles <- list.files(
	system.file("extdata", package = "CAGEr"),
	"ctss$", full.names = TRUE
)

ce <- CAGEexp(
	genomeName = "BSgenome.Drerio.UCSC.danRer7",
	inputFiles = inputFiles,
	inputFilesType = "ctss",
        sampleLabels = sub( ".chr17.ctss", "", basename(inputFiles))
)

getCTSS(ce)
```

### Normalizing TSS Counts

[CAGEr](https://bioconductor.org/packages/release/bioc/html/CAGEr.html) provides the option to normalize TSS counts using a power-law based normalization.

```
normalizeTagCount(
	ce, method = "powerLaw", fitInRange = c(5, 1000),
	alpha = 1.2, T = 5*10^4)
)
```

## TSS Clustering

TSSs can be clustered into TSRs, and various CAGE based TSR metrics can be calculated.

### Clustering TSSs.

[CAGEr](https://bioconductor.org/packages/release/bioc/html/CAGEr.html) provides various clustering algorithms, such as nieve distance clustering and clustering algorithms such as paraclu. This example uses nieve distance clustering.

```
clusterCTSS(
	object = ce, threshold = 1, thresholdIsTpm = TRUE, nrPassThreshold = 1,
        method = "distclu", maxDist = 20, removeSingletons = TRUE,
        keepSingletonsAbove = 5
)
```

## Consensus Clusters

All TSRs can be merged into consensus clusters.

```
cumulativeCTSSdistribution(ce, clusters = "tagClusters", useMulticore = T)
quantilePositions(ce, clusters = "tagClusters", qLow = 0.1, qUp = 0.9)
aggregateTagClusters(ce, tpmThreshold = 5, qLow = 0.1, qUp = 0.9, maxDist = 100)
```


