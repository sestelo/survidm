Beran <-
  function(time, status, covariate, delta, x, y, kernel = "gaussian", bw,
           lower.tail = FALSE) {


    if (missing(time)|is.null(time)) {
      stop("Argument \"time\" is missing, with no default")
    }
    if (missing(status)|is.null(time)) {
      stop("Argument \"status\" is missing, with no default")
    }


    spa <- NULL
    len <- length(time)
    if (missing(delta)) delta <- rep(1, len)


    res <- .C( "SurvBeranKernel", as.double(time), as.integer(status),
               as.double(covariate), as.integer(delta), as.integer(len),
               as.double(y), as.double(x), as.double(bw),
               as.character(kernel), p = as.double(1))$p
    if (lower.tail == TRUE) res <- 1 - res
    return(res)}
