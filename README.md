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

generate correlation matrix

```
corr_plot <- plot_corr_matrix(exp, corr_metric="pearson")

ggsave("corr.pdf", corr_plot, device=cairo_pdf, height=5.5, width=7)
```
