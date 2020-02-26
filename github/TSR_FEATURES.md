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
