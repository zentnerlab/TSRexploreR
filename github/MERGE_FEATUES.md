# Merge Features


```
library("tsrexplorer")
library("magrittr")

TSSs <- system.file("extdata", "S288C_TSSs.RDS", package = "tsrexplorer")
TSSs <- readRDS(TSSs)

exp <- tsr_explorer(TSSs)
```

```
exp <- format_counts(exp, data_type = "tss") %>%
        tss_clustering(threshold = 3, max_distance = 25)
```

```
exp <- merge_samples(
	exp, data_type = "tsr",
	sample_list = list(
		Untreated = c("S288C_WT_1", "S288C_WT_2", "S288C_WT_3"),
		Diamide = c("S288C_D_1", "S288C_D_2", "S288C_D_3")
	)
)
```
