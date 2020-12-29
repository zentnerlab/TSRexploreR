
  ## Create example data.
  TSSs <- readRDS(system.file("extdata", "S288C_TSSs.RDS", package="TSRexploreR"))
  TSSs <- map(TSSs, function(x) {
    x <- x[seqnames(x) == "I"]
    return(x)
  }) 
