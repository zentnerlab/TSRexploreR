# Standard Analysis

Transcription Start Sites (TSSs) represent the first base to be transcribed by RNA Polymerase.
Most genes tend to have a heterogenous collection of TSSs as opposed to a single position.
This is interesting because choice of TSS can affect gene isoform, stability, and translational efficiency.
Furthermore, TSS choice has been implicated in develoment, homseostasis, and disease.

tsrexploreR has a series of analysis and plotting functions to allow deep exploration of TSSs.

## Preparing Data

This example will use a set of *S. cerevisiae* TSSs collected using the STRIPE-seq method.
There are many ways to import TSSs into tsrexplorer.
This example data uses a named list of GRanges as imput into the function that creates the tsrexplorer object.

```
library("tsrexplorer")
library("magrittr")

TSSs <- system.file("extdata", "S288C_TSSs.RDS", package = "tsrexplorer")
TSSs <- readRDS(TSSs)

# Keep only the 3 WT samples for now.
TSSs <- names(TSSs) %>%
	stringr::str_detect("WT") %>%
	purrr::keep(TSSs, .)

exp <- tsr_explorer(TSSs)
```

## Initial TSS Processing

After the TSSs are loaded into the tsrexplorer object,
there are a few processing steps to go through to get the data ready for analysis.

### Format TSSs

The first step is to convert the TSSs into a table that will facilitate downstream analysis.

```
exp <- format_counts(exp, data_type = "tss")
```

### Normalize TSSs

The next step is to Counts Per Million (CPM) normalize the TSSs.
This step is optional, and if the counts you inputed were normalized already this step can safely be skipped.

```
exp <- cpm_normalize(exp, data_type = "tss")
```

## TSS Annotation

After formatting the counts and optionally CPM normalizing them, the TSSs will be annotated relative to known features.
This function takes either the path and file name of a 'GTF' or 'GFF' file, or a 'TxDb' package from bioconductor.
The example below uses a 'GTF' file from Ensembl (R64-1-1 Release 99), and will annotate each TSS to the closest transcript.

```
annotation <- system.file("extdata", "S288C_Annotation.gtf", package = "tsrexplorer")

exp <- annotate_features(
        exp, annotation_data = annotation,
        data_type = "tss", feature_type = "transcript"
)
```

### Naive Threshold Exploration

In most TSS mapping methods there are low count TSSs that are difficult to distinguish between signal and noise.
A good amount of this background can be removed by requiring a certain number of reads per TSS.
If a genome annotation is available, a thresholding plot may help in picking this read threshold.

```
threshold_data <- explore_thresholds(exp, max_threshold = 25)

p <- plot_threshold_exploration(threshold_data, ncol = 3, point_size = 0.5) +
	ggplot2::theme(text = element_text(size = 4), legend.key.size = unit(0.3, "cm"))

ggsave("tss_thresholding.png", plot = p, device = "png", type = "cairo", height = 1.25, width = 5)
```

![tss_thresholds](../inst/images/tss_thresholding.png)
