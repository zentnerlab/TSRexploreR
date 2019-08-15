# tsrexplorer

## Installing TSRexplorer

**Create conda environment**
```
conda create -n tsrexplorer -y -c conda-forge -c bioconda \
r-tidyverse \
r-devtools \
r-ggseqlogo \
bioconductor-genomicranges \
bioconductor-genomicfeatures \
bioconductor-biostrings \
bioconductor-rsamtools \
bioconductor-chipseeker \
bioconductor-edger
```

**Install latest version of tsrexplorer**
```
devtools::install_github("rpolicastro/tsrexplorer")
```

## Using TSRexplorer

### Preparing TSRexplorer

**Load tsrexplorer**

```
library("tsrexplorer")
```

**Load example data**

```
TSSs <- system.file("extdata", "yeast_TSSs.RDS", package = "tsrexplorer")
TSSs <- readRDS(TSSs)

TSRs <- system.file("extdata", "yeast_TSRs.RDS", package = "tsrexplorer")
TSRs <- readRDS(TSRs)

annotation <- system.file("extdata", "yeast_annotation.gtf", package="tsrexplorer")
assembly <- system.file("extdata", "yeast_assembly.fasta", package="tsrexplorer")
```

**create tsr object**

```
exp <- tsr_explorer(TSSs, TSRs)
```

## TSS Analysis

### Count Normalization and Correlation

**tmm normalize counts**

```
exp <- tss_normalization(exp)
```

**tss correlation matrix**

```
p <- plot_tss_corr(exp, corr_metric="pearson")

ggsave("tss_corr.png", plot = p, device = "png", type = "cairo", height = 3.5, width = 5)
```
![tss_corr_plot](./inst/images/tss_corr.png)

**generate tss scatter plots**

```
p <- plot_tss_scatter(exp, sample_1 = "S288C_WT_100ng_1", sample_2 = "S288C_WT_100ng_2")

ggsave("tss_scatter.png", plot = p, device = "png", type = "cairo", height = 4, width = 4)
```
![tss_scatter_plot](./inst/images/tss_scatter.png)

### TSS Annotation

```
exp <- tss_annotation(exp, annotation_file = annotation, feature_type = "transcript")
```

### TSS Average Plot and Heatmap

**TSS average plot**

```
p <- plot_tss_average(exp, sample = "S288C_WT_100ng_1", threshold = 3)

ggsave("tss_average_plot.png", plot = p, device = "png", type = "cairo", height = 4, width = 4)
```

![tss_average_plot](./inst/images/tss_average_plot.png)

**TSS heatmap**

```
count_matrix <- tss_heatmap_matrix(exp, sample = "S288C_WT_100ng_1", threshold = 3, anno_type = "geneId")

p <- plot_tss_heatmap(count_matrix)

ggsave("tss_heatmap.png", plot = p, device = "png", type = "cairo", height = 4, width = 3)
```

![tss_heatmap](./inst/images/tss_heatmap.png)

### TSS Motif and Base Composition

**TSS sequence logo**

```
seqs <- tss_sequences(exp, sample = "S288C_WT_100ng_1", genome_assembly = assembly, threshold = 3)

p <- plot_sequence_logo(seqs)

ggsave("tss_seq_logo.png", plot = p, device = "png", type = "cairo", height = 1.5, width = 4)

```

![tss_sequence_logo](./inst/images/tss_seq_logo.png)

**TSS base color map**

```
p <- plot_sequence_colormap(seqs)

ggsave("tss_seq_colormap.png", plot = p, device = "png", type = "cairo", height = 4, width = 4)
```
![tss_sequence_colormap](./inst/images/tss_seq_colormap.png)

**TSS dinucleotide frequencies**

```
frequencies <- dinucleotide_frequencies(exp, sample = "S288C_WT_100ng_1", genome_assembly = assembly, threshold = 3)

p <- plot_dinucleotide_frequencies(frequencies)

ggsave("tss_dinucleotide_frequencies.png", plot = p, device = "png", type = "cairo", height = 3, width = 4)
```

![tss_dinucelotide_frequencies](./inst/images/tss_dinucleotide_frequencies.png)

### Misc TSS Plots

**Average distance of dominant TSS**

```
dominant <- dominant_tss(exp, sample = "S288C_WT_100ng_1", threshold = 3, feature_type = "geneId")

p <- plot_dominant_tss(dominant, upstream = 500, downstream = 500)

ggsave("dominant_tss.png", plot = p, device = "png", type = "cairo", height = 4, width = 4)
```

![dominant_tss](./inst/images/dominant_tss.png)

**Max UTR Length**
``` 
max <- max_utr(exp, sample = "S288C_WT_100ng_1", threshold = 3, feature_type = "geneId")

p <- plot_max_utr(max)

ggsave("max_utr.png", plot = p, device = "png", type = "cairo", height = 4, width = 4)
```

![max_utr](./inst/images/max_utr.png)

## TSR Analysis

### Count Normalization and Correlation

**TMM normalize counts**

```
exp <- tsr_normalization(exp)
```

**TSR correlation matrix**

```
p <- plot_tsr_corr(experiment, corr_metric = "pearson")

ggsave("tsr_corr.png", plot = p, device = "png", type = "cairo", height = 3.5, width = 5)
```

![tsr_correlation](./inst/images/tsr_corr.png)

**TSR scatter plot**

```
p <- plot_tsr_scatter(exp, sample_1 = "S288C_WT_100ng_1", sample_2 = "S288C_WT_100ng_2")

ggsave("tsr_scatter.png", plot = p, device = "png", type = "cairo", height = 4, width = 4)
```

![tsr_scatter](./inst/images/tsr_scatter.png)

### TSR Annotation and Genomic Distribution

**TSR Annotation**

```
exp <- tsr_annotation(exp, annotation_file = annotation, feature_type = "transcript")
```

**TSR Genomic Distribution**

```
tsr_distribution <- tsr_genomic_distribution(exp, sample = "S288C_WT_100ng_1")

p <- plot_tsr_genomic_distribution(tsr_distribution)

ggsave("tsr_genomic_distribution.png", plot = p, device = "png", type = "cairo", height = 2, width = 5)
```

![tsr_genomic_distribution](./inst/images/tsr_genomic_distribution.png)

### TSR Average Plot and Heatmap

**TSR Average Plot**

```
p <- plot_tsr_average(exp, sample = "S288C_WT_100ng_1")

ggsave("tsr_average_plot.png", plot = p, device = "png", type = "cairo", height = 4, width = 4)
```

![tsr_average_plot](./inst/images/tsr_average_plot.png)

**TSR Heatmap**

```
counts <- tsr_heatmap_matrix(exp, sample = "S288C_WT_100ng_1", feature_type = "transcriptId")

p <- plot_tsr_heatmap(counts)

ggsave("tsr_heatmap.png", plot = p, device = "png", type = "cairo", height = 4, width = 3)
```

![tsr_heatmap](./inst/images/tsr_heatmap.png)

## Differential TSRs

(Work in Progress)

**Find Differential TSRs**

```
edger_model <- fit_edger_model(
	exp, 
	samples = c(
		"S288C_WT_100ng_1",
		"S288C_WT_100ng_2",
		"S288C_WT_100ng_3",
		"S288C_Diamide_100ng_1",
		"S288C_Diamide_100ng_2",
		"S288C_Diamide_100ng_3"
	),
	groups = c(1, 1, 1, 2, 2, 2)
)

diff_tsrs <- differential_tsrs(edger_model, comparisons = c(1, 2))

**Annotate Differential TSRs**

```
annotated_diff_tsrs <- annotate_differential_tsrs(diff_tsrs, annotation_file = annotation, feature_type = "trasncript")
```
