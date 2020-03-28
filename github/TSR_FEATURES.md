# TSR Features

## Loading TSSs

```
library("tsrexplorer")

TSSs <- system.file("extdata", "S288C_TSSs.RDS", package = "tsrexplorer")
TSSs <- readRDS(TSSs)

# Keep only the 3 WT samples for now.
TSSs <- names(TSSs) %>%
        stringr::str_detect("WT") %>%
        purrr::keep(TSSs, .)

exp <- tsr_explorer(TSSs)
```

## Formatting TSSs

### Processing TSSs

```
exp <- format_counts(exp, data_type = "tss")
```

### Clustering TSSs

```
exp <- tss_clustering(exp, threshold = 3, max_distance = 25)
```

## Associate TSRs with TSSs

```
exp <- associate_with_tsr(
	exp, use_sample_sheet = FALSE,
	sample_list <- list(
		"S288C_WT_1" = "S288C_WT_1",
		"S288C_WT_2" = "S288C_WT_2",
		"S288C_WT_3" = "S288C_WT_3"
	)
)
```

## TSR Metrics

```
exp <- tsr_metrics(exp)
```

## Annotate TSSs

```
annotation <- system.file("extdata", "S288C_Annotation.gtf", package = "tsrexplorer")

exp <- annotate_features(
        exp, annotation_data = annotation,
        data_type = "tss", feature_type = "transcript"
)
```

## Annotate TSRs

```
annotation <- system.file("extdata", "S288C_Annotation.gtf", package = "tsrexplorer")

exp <- annotate_features(
        exp, annotation_data = annotation,
        data_type = "tsr", feature_type = "transcript"
)
```

## TSS Plots

### Genomic Distribution

```
conditions <- list(grouping = "shape_class")
tss_distribution <- genomic_distribution(exp, data_type = "tss", threshold = 3, data_conditions = conditions)
```

### Average Plots

```
NULL
```

### Heatmaps

```
conditions <- list(order_by = "score", order_group = "transcriptId", grouping = "shape_class")
mat <- tss_heatmap_matrix(
	exp, upstream = 250, downstream = 250,
	threshold = 3, data_conditions = conditions
)

p <- plot_heatmap(mat)
```
