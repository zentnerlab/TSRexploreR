
# TSS Shifting

```
library("TSRexploreR")
library("magrittr")

TSSs <- system.file("extdata", "S288C_TSSs.RDS", package = "TSRexploreR")
TSSs <- readRDS(TSSs)

annotation <- system.file("extdata", "S288C_Annotation.gtf", package = "TSRexploreR")
```

```
exp <- tsr_explorer(TSSs) %>%
	## Merge the TSSs so they can be later linked to the merged TSRs.
	format_counts(data_type = "tss") %>%
	merge_samples(data_type = "tss", threshold = 3, sample_list = list(
		Untreated = c("S288C_WT_1", "S288C_WT_2", "S288C_WT_3"),
		Diamide = c("S288C_D_1", "S288C_D_2", "S288C_D_3")
	)) %>%
	format_counts(data_type = "tss", samples = c("Untreated", "Diamide")) %>%
	## Cluster the unmerged TSSs into TSRs.
        tss_clustering(threshold = 3, max_distance = 25, samples = c(
		"S288C_WT_1", "S288C_WT_2", "S288C_WT_3",
		"S288C_D_1", "S288C_D_2", "S288C_D_3"
	)) %>%
	## Merge the TSRs called from the unclustered TSSs.
	merge_samples(data_type = "tsr", sample_list = list(
		Untreated = c("S288C_WT_1", "S288C_WT_2", "S288C_WT_3"),
		Diamide = c("S288C_D_1", "S288C_D_2", "S288C_D_3")
	)) %>%
	format_counts(data_type = "tsr", samples = c("Untreated", "Diamide")) %>%
	## Associate the merged TSSs from the beginning to the merged TSRs.
	associate_with_tsr(use_sample_sheet = FALSE, sample_list = list(
		Untreated = "Untreated", Diamide = "Diamide"
	))

```
## Calculate the shifting scores and return as table.
shift_results <- tss_shift(
	exp, compare_samples = c("Untreated", "Diamide"), min_distance = 100,
	min_threshold = 10, n_resamples = 1000L
)
```
