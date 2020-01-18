
# Import Data

While tsrexplorer can accept a named list of GRanges as input for both TSSs and TSRs, 
A data import function is provided for convenience to allow compatability between various sources and packages.

The currently supported sources of data input include:

* Bedgraph and bigwig files for TSSs, and bed files for TSRs.
* TSRchitect tssObject containing TSSs or TSRs.
* CAGEr CAGEexp object containing CTSSs or tag clusters.

## Prepare tsrexplorer Object

Before importing any data, you must initialize an empty tsrexplorer object to store the imported data.

```
exp <- tsr_explorer()
```

## Import TSSs

### Bedgraph or Bigwig files

TSS data is most often presented in bedgraph or bigwig files split into plus and minus strand files.
In order to import this data you must first create a sample sheet with three columns 'sample_name', 'pos', and 'neg'.
The first column 'sample_name' is the name you would like associated to that sample.
the second column 'pos' is the path to the positive strand file.
Finally, the third column 'neg' is the path to the negatively stranded file.
The sample sheet can be a tab delimited file, or a data.frame.

```
exp <- tss_import(exp, sample_sheet = "example_path/sample_sheet.tsv")

or

exp <- tss_import(exp, sample_sheet = df)
```

### TSRchitect derived TSSs

TSSs can be extracted from TSRchitect objects (tssObject).
In order to import the TSSs, they must first be derived from 'inputToTSS',
followed by merging and summing overlapping TSSs with 'processTSS'.


```
exp <- tss_import(exp, tssObject)
``` 
