LLW <-
  function(x, kernel = "gaussian", bw, t1){
    nobs <- length(x)
    return(.C("LLWeightsKernel",as.double(x), as.integer(nobs),
              as.double(t1), as.double(bw), as.character(kernel),
              w = double(nobs) )$w)}
