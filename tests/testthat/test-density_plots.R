
source("setup.R")

test_that("Density plots.", {

  tsre <- TSSs[1] %>%
    tsr_explorer(genome_annotation=annotation) %>%
    format_counts(data_type="tss") %>%
    normalize_counts(data_type="tss", method="CPM") %>%
    annotate_features(data_type="tss") %>%
    tss_clustering %>%
    annotate_features(data_type="tsr")

  ## TSS density plots.
  plot_density(tsre, data_type="tss") %>%
    expect_s3_class("ggplot")
  plot_density(tsre, data_type="tss", consider_score=TRUE) %>%
    expect_s3_class("ggplot")
  plot_density(tsre, data_type="tss", threshold=3) %>%
    expect_s3_class("ggplot")
  plot_density(tsre, data_type="tss", use_normalized=TRUE) %>%
    expect_s3_class("ggplot")
  plot_density(tsre, data_type="tss", exclude_antisense=FALSE) %>%
    expect_s3_class("ggplot")

  ## TSR density plots.
  plot_density(tsre, data_type="tsr") %>%
    expect_s3_class("ggplot")
  plot_density(tsre, data_type="tsr", consider_score=TRUE) %>%
    expect_s3_class("ggplot")
  plot_density(tsre, data_type="tsr", threshold=3) %>%
    expect_s3_class("ggplot")
  plot_density(tsre, data_type="tsr", use_normalized=TRUE) %>%
    expect_s3_class("ggplot")
  plot_density(tsre, data_type="tsr", exclude_antisense=FALSE) %>%
    expect_s3_class("ggplot")

})
