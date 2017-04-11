PKMW <-
  function(time, status){
    t1 <- max(time)
    len <- length(time)
    status1 <- fitted(glm(status ~ time, family = binomial))
    res <- .C( "WeightsKaplanMeierSortEx", time=as.double(time),
               status <- as.double(status1), as.integer(len), as.double(t1),
               weights = double(len) )
    return(res$weights)
  }
