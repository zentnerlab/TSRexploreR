
  ## Create example data.
  TSSs <- readRDS(system.file("extdata", "S288C_TSSs.RDS", package="TSRexploreR"))
  TSSs <- map(TSSs, function(x) {
    x <- x[seqnames(x) == "I"]
    return(x)
  })

  ## Genome assembly an annotation files.
  assembly <- system.file("extdata", "S288C_Assembly.fasta", package="TSRexploreR")
  annotation <- system.file("extdata", "S288C_Annotation.gtf", package="TSRexploreR")
