KM <-
  function(time, status, t)
  {
    if (missing(time)) stop("Argument 'time' is missing with no default")
    if (missing(status)) stop("Argument 'status' is missing with no default")
    if (missing(t)) stop("Argument 't' is missing with no default")
    return( .C( "KaplanMeierValueSort", as.double(time), as.integer(status),
                as.integer( length(time) ), as.double(t), as.double(1) )[[5]] )
  }
