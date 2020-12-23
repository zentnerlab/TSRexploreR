# Export Data

[TSRexploreR](https://github.com/zentnerlab/TSRexploreR) provides flexible options to export TSSs and TSRs.

The currently supported export formats are supported"

* TSSs: Bedgraphs or tables
* TSRs: Beds or tables.

## Export TSSs.

TSSs are often saved as bedgraph files split by the positive and minus strands.

We will first prepare an example [TSRexploreR](https://github.com/zentnerlab/TSRexploreR) object with TSSs we wish to export.

```
# Load two example TSS files.
TSSs <- system.file("extdata", "S288C_TSSs.RDS", package = "TSRexploreR")
TSSs <- readRDS(TSSs)[c("S288C_WT_1", "S288C_WT_2")]

# Keep only the TSSs from chromosome I (for example purposes).
TSSs <- map(TSSs, function(x) x[seqnames(x) == "I"])

# Prepare the TSSs.

exp <- tsr_explorer(TSSs)
exp <- format_counts(exp, data_type = "tss")
```

Now that we have a [TSRexploreR](https://github.com/zentnerlab/TSRexploreR) object with TSSs, we can export them as bedgraph files.
A reminder that they can also be exported as a tabl delimited table.

```
tss_export(exp)
```

## Export TSRs.

TSRs are often saved as bed files.

The example TSSs from above must first be clustered.

```
exp <- tss_clustering(exp, threshold = 3, max_distance = 25)
```

You can now save the TSRs as bed files.

```
tsr_export(exp)
```
