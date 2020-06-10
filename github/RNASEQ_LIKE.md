# RNA-seq Like Analaysis

Prepare the TSSs for feature counting.

```
library("tsrexplorer")
library("magrittr")

TSSs <- system.file("extdata", "S288C_TSSs.RDS", package = "tsrexplorer")
TSSs <- readRDS(TSSs)

exp <- tsr_explorer(TSSs)
exp <- format_counts(exp, data_type = "tss")
```

Annotate the TSSs to the nearest transcript.

```
annotation <- system.file("extdata", "S288C_Annotation.gtf", package = "tsrexplorer")
exp <- annotate_features(
        exp, annotation_data = annotation,
        data_type = "tss", feature_type = "transcript"
)
```

Summarize the reads per gene.

```
exp <- count_features(exp, "tss")
```

TMM normalize for correlation.

```
exp <- count_matrix(exp, "tss_features")
exp <- tmm_normalize(exp, "tss_features", threshold = 3, n_samples = 3)
```

Correlation matrix plot.

```
p <- plot_correlation(exp, data_type = "tss_features", font_size = 2, pt_size = 0.4) +
        ggplot2::theme_bw() +
        ggplot2::theme(text = element_text(size = 3), panel.grid = element_blank())

ggsave("tss_feature_correlation.png", plot = p, device = "png", type = "cairo", height = 3, width = 3)
```

![tss_feature_corr_plot](../inst/images/tss_feature_correlation.png)

Get number of detected features.

```
features <- detect_features(exp, data_type = "tss_features", threshold = 3)

p <- plot_detected_features(features) +
        ggplot2::theme(text = element_text(size = 3))

ggsave("tss_features_feature_plot.png", plot = p, device = "png", type = "cairo", height = 1, width = 1.5)
```

![tss_features_feature_plot](../inst/images/tss_features_feature_plot.png)
