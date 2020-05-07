
# Import Data

While tsrexplorer can accept a named list of GRanges as input for both TSSs and TSRs, 
A data import function is provided for convenience to allow compatability between various sources and packages.

The currently supported sources of data input include:

* Bedgraph and bigwig files for TSSs, and bed files for TSRs.
* TSRchitect tssObject containing TSSs or TSRs.
* CAGEr CAGEexp object containing CTSSs or tag clusters.

## Import TSSs

### Bedgraph or Bigwig Files

TSS data is most often presented as bedgraphs or bigwigs split into plus and minus strand files.
In order to import this data you must first create a sample sheet with three columns 'sample_name', 'file_1', and 'file_2'.
The first column 'sample_name' is the name you would like associated to that sample.
the second column 'file_1' is the path to the positive strand file.
Finally, the third column 'file_2' is the path to the negatively stranded file.
The sample sheet can be a tab delimited file, or a data.frame.

First, we will create some sample bedgraph files to import.

```
library("GRanges")
library("rtracklayer")
library("purrr")

# Load two example TSS files.
TSSs <- system.file("extdata", "S288C_TSSs.RDS", package = "tsrexplorer")
TSSs <- readRDS(TSSs)[c("S288C_WT_1", "S288C_WT_2")]

# Keep only the TSSs from chromosome I (for example purposes).
TSSs <- map(TSSs, function(x) x[seqnames(x) == "I"])

# Save the file as bedgraphs.
iwalk(TSSs, function(x, y) {
	pos <- x[strand(x) == "+"]
	export(pos, str_c(y, "_pos.bedgraph"))
	neg <- x[strand(x) == "-"]
	export(neg, str_c(y, "_neg.bedgraph"))
})

```

Now that there are example bedgraphs to import, we can make a sample sheet.

```
sample_sheet <- data.frame(
	"sample_name" = c("S288C_WT_1", "S288C_WT_2"),
	"file_1" = c("S288C_WT_1_pos.bedgraph", "S288C_WT_2_pos.bedgraph"),
	"file_2" = c("S288C_WT_1_neg.bedgraph", "S288C_WT_2_neg.bedgraph"),
	stringsAsFactors = FALSE
)
```

Finally, we can import the bedgraphs and add them to the tsrexplorer object.
A reminder that the sample sheet can also be a tab delimited file of the same format as the
data.frame created above.

```
# Initialize empty tsrexplorer object.
exp <- tsr_explorer()

# Import the bedgraphs.
exp <- tss_import(exp, sample_sheet)
```

### TSRchitect Derived TSSs

TSSs can be extracted from TSRchitect objects (tssObject).
In order to import the TSSs, they must first be derived from `inputToTSS`,
followed by merging and summing overlapping TSSs with `processTSS`.


```
exp <- tss_import(exp, tssObject)
```

## Import TSRs

### Bed Files

TSR data is most often represented as a bed file.
In order to import this data, you must first create a sample sheet with two columns, 'sample_name', and 'bed'.
The first column, 'sample_name' is the name you would like associated to that sample.
The second and last column 'bed' is the path to the bed files to import.
The sample sheet can be a tab delimited file, or a data.frame.

```
# Sample sheet as tab delimited file
exp <- tsr_import(exp, "example_path/sample_sheet.tsv")

or

# Sample sheet as data.frame 'df'.
exp <- tsr_import(exp, sample_sheet = df)
```

### TSRchitect Derived TSRs

TSRs can be extracted from TSRchitect objects (tssObject).
In order to import the TSRs, they must first be derived from `determineTSR`.

```
exp <- tsr_import(exp, tssObject)
```
