# tsrexplorer

## Installing TSRexplorer

```
devtools::install_github("rpolicastro/tsrexplorer")
```

## Using TSRexplorer

Load tsrexplorer

```
library("tsrexplorer")
```

create tsr object

```
exp <- tsr_explorer(TSSs)
```

tmm normalize counts

```
exp <- normalize(exp)
```

generate tss correlation matrix

```
corr_plot <- plot_tss_corr(exp, corr_metric="pearson")

ggsave("tss_corr.pdf", corr_plot, device=cairo_pdf, height=5.5, width=7)
```

generate tss scatter plots

```
scatter_plot <- plot_tss_scatter(exp, "sample_1", "sample_2")

ggsave("tss_scatter.png", scatter_plot, device="png", type="cairo", height=3, width=3)
```

