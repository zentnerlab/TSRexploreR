
source("setup.R")

test_that("Check TSS and TSR correlations.", {

  tsre <- TSSs[1:2] %>%
    tsr_explorer %>%
    format_counts(data_type="tss") %>%
    normalize_counts(method="CPM") %>%
    tss_clustering

  ## TSS correlation.
  plot_correlation(tsre, data_type="tss") %>%
    expect_s4_class("Heatmap")
  plot_correlation(tsre, data_type="tss", use_normalized=FALSE) %>%
    expect_s4_class("Heatmap")
  plot_correlation(tsre, data_type="tss", cluster_samples=TRUE) %>%
    expect_s4_class("Heatmap")
  plot_correlation(tsre, data_type="tss", show_values=FALSE) %>%
    expect_s4_class("Heatmap")
  plot_correlation(tsre, data_type="tss", return_matrix=TRUE) %>%
    {expect_true(is(., "matrix"))}

  ## TSR correlation.
  plot_correlation(tsre, data_type="tsr") %>%
    expect_s4_class("Heatmap")
  plot_correlation(tsre, data_type="tsr", use_normalized=FALSE) %>%
    expect_s4_class("Heatmap")
  plot_correlation(tsre, data_type="tsr", return_matrix=TRUE) %>%
    {expect_true(is(., "matrix"))}

})
