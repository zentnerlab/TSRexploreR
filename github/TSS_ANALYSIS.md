# TSS Analysis

Transcription Start Sites (TSSs) represent the first base to be transcribed by RNA Polymerase.
Most genes tend to have a heterogenous collection of TSSs as opposed to a single position.
This is interesting because choice of TSS can affect gene isoform, stability, and translational efficiency.
Furthermore, TSS choice has been implicated in develoment, homseostasis, and disease.

tsrexploreR has a series of analysis and plotting functions to allow deep exploration of TSSs.

## Preparing Data

This example will use a set of *S. cerevisiae* TSSs collected using the STRIPE-seq method.
There are many ways to import TSSs into tsrexplorer.
This example data uses a named list of GRanges as imput into the function that creates the *tsrexplorer object*.

```
library("tsrexplorer")

TSSs <- system.file("extdata", "S288C_TSSs.RDS", package = "tsrexplorer")
TSSs <- readRDS(TSSs)

exp <- tsr_explorer(TSSs)
```

## Processing of TSSs

After the TSSs are loaded into the *tsrexplorer object*,
there are a few processing steps to go through to get the data ready for analysis.
These steps include converting the TSSs into a *RangedSummarizedExperiment*,
CPM normalizing the TSSs, and then marking the dominant TSS.

The first step is to convert the TSSs into a *RangedSummarizedExperiment*.
This object type is a convenient container that stores range, count, and sample information.

```
exp <- format_counts(exp, data_type = "tss")
```

The next step is to Counts Per Million (CPM) normalize the TSSs and store them as an additional assay with the original counts.
This step is optional, and if the counts you inputed were normalized already this step can safely be skipped.

```
exp <- cpm_normalize(exp, data_type = "tss")
```

After formatting the counts and optionally CPM normalizing them, the TSSs will be annotated relative to known features.
This function takes either the path and file name of a *.GTF* or *.GFF* file, or a *BSgenome* package from bioconductor.
The annotation information will be added onto the range of the *RangedSummarizedExperiment*.
In the example below each TSS is annotated relative to the closest transcript.

```
library("TxDb.Scerevisiae.UCSC.sacCer3.sgdGene")

exp <- annotate_features(
        exp, annotation_data = TxDb.Scerevisiae.UCSC.sacCer3.sgdGene,
        data_type = "tss", feature_type = "transcript"
)
```
