NWW <-
  function(covariate, x, kernel="gaussian", bw) {
    spa <- NULL
    len <- length(covariate)
    listg <- .C( "NWWeightsKernel", as.double(covariate), as.integer(len),
                 as.double(x), as.double(bw), as.character(kernel),
                 weight = double(len) )
    return(listg$weight)}
