
# TSS Shifting

```
library("tsrexplorer")
library("magrittr")

TSSs <- system.file("extdata", "S288C_TSSs.RDS", package = "tsrexplorer")
TSSs <- readRDS(TSSs)
```

```
exp <- tsr_explorer(TSSs) %>%
	format_counts(data_type = "tss") %>%
	merge_samples(data_type = "tss", threshold = 3, sample_list = list(
		Untreated = c("S288C_WT_1", "S288C_WT_2", "S288C_WT_3"),
		Diamide = c("S288C_D_1", "S288C_D_2", "S288C_D_3")
	)) %>%
	format_counts(data_type = "tss", samples = c("Untreated", "Diamide")) %>%
        tss_clustering(threshold = 3, max_distance = 25, samples = c(
		"S288C_WT_1", "S288C_WT_2", "S288C_WT_3",
		"S288C_D_1", "S288C_D_2", "S288C_D_3"
	)) %>%
	merge_samples(data_type = "tsr", sample_list = list(
		Untreated = c("S288C_WT_1", "S288C_WT_2", "S288C_WT_3"),
		Diamide = c("S288C_D_1", "S288C_D_2", "S288C_D_3")
	)) %>%
	format_counts(data_type = "tsr", samples = c("Untreated", "Diamide")) %>%
	associate_with_tsr(use_sample_sheet = FALSE, sample_list = list(
		Untreated = "Untreated", Diamide = "Diamide"
	))
