
source("setup.R")

test_that("TSS and TSR heatmaps.", {

  ## TSS heatmaps.
  tsre <- TSSs["S288C_WT_1"] %>%
    tsr_explorer(genome_annotation=annotation) %>%
    prepare_counts(data_type="tss") %>%
    annotate_features(data_type="tss")

  # TSS heatmap.
  plot_heatmap(tsre, data_type="tss") %>%
    expect_s3_object("ggplot")
  # TSS heatmap with rasterization.
  plot_heatmap(tsre, data_type="tss", rasterization=TRUE) %>%
    expect_s3_object("ggplot")

  ## TSR heatmaps.
  tsre <- tsre %>%
    tss_clustering(threshold=3) %>%
    annotate_features(data_type="tsr")

  # TSR heatmap.
  plot_heatmap(tsre, data_type="tsr") %>%
    expect_s3_object("ggplot")
  # TSR heatmap with rasterization.
  plot_heatmap(tsre, data-type="tsr", rasterization=TRUE) %>%
    expect_s3_object("ggplot")

})
