tpLIDA <-
  function(object, s, conf = FALSE, n.boot = 199, conf.level = 0.95,
           cluster = FALSE, ncores = NULL)
  {
    if (missing(object))
      stop("Argument 'object' is missing, with no default")
   # if (!inherits(object, "survIDM")) stop("'object' must be of class 'survIDM'")
    if (missing(s))
      s <- 0

    ptimes <- which(object[[1]]$event1 == 1 | object[[1]]$event == 1)

    if (s > max(object[[1]]$time1)) stop("The value of 's' is too large")
    if (s < 0) stop("'s' must be nonnegative")
    if (length(s) > 1) stop("Length of 's' must be 1")

    t1 <- object[[1]]$time1[object[[1]]$event1 == 1]
    t2 <- object[[1]]$Stime[object[[1]]$event == 1]
    t <- c(s, t1, t2)
    t <- t[t >= s]
    t <- sort(unique(t))

    n <- length(object[[1]]$Stime)

    p0 <- which(object[[1]]$time1 <= s)
    p1 <- which(object[[1]]$time1 <= s & object[[1]]$Stime > s)
    p1new <- which(object[[1]]$time1 <= s & object[[1]]$Stime > s & object[[1]]$event == 1)
    state1 <- length(p1new)

    p00 <- rep(NA, length(t))
    p01 <- rep(NA, length(t))
    p02 <- rep(NA, length(t))
    p11 <- rep(NA, length(t))
    p12 <- rep(NA, length(t))

    den <- KM(object[[1]]$time1, object[[1]]$event1, s)
    den2 <- sum(KMW(object[[1]]$Stime, object[[1]]$event)[p1])
    kmw1 <- KMW(object[[1]]$Stime, object[[1]]$event)

    for (k in 1: length(t)) {
      p2 <- which(object[[1]]$time1 > s & object[[1]]$time1 <= t[k] & object[[1]]$Stime > t[k])
      p3 <- which(object[[1]]$time1 <= s & object[[1]]$Stime > t[k])
      n.p2 <- length(p2)
      n.p3 <- length(p3)
      p00[k] <- KM(object[[1]]$time1, object[[1]]$event1, t[k]) / den
      ifelse(n.p2 == 0, p01[k] <- 0, p01[k] <- sum(kmw1[p2]) / den)
      p02[k] <- 1 - p00[k] - p01[k]
      if(p01[k] >  1) p01[k] <- 1
      if(p02[k] <  0) p02[k] <- 0
      if (state1 > 0) {
        ifelse(n.p3 == 0, p11[k] <- 0, p11[k]	<- sum(kmw1[p3]) / den2)
        p12[k] <- 1 - p11[k]
      }
    }

    resu <- data.frame(cbind(t, p00, p01, p02, p11, p12))
    names(resu) <- c("t", "p00", "p01", "p02", "p11", "p12")


    p00.ci <- matrix(NA, length(t), 2)
    p01.ci <- matrix(NA, length(t), 2)
    p02.ci <- matrix(NA, length(t), 2)
    p11.ci <- matrix(NA, length(t), 2)
    p12.ci <- matrix(NA, length(t), 2)

    p00.std <- matrix(NA, length(t), 2)
    p01.std <- matrix(NA, length(t), 2)
    p02.std <- matrix(NA, length(t), 2)
    p11.std <- matrix(NA, length(t), 2)
    p12.std <- matrix(NA, length(t), 2)


    if (conf == TRUE) {
      res.ci <- array(NA, dim=c(length(t), n.boot, 5))

      for (j in 1:n.boot){
        xx <- sample.int(n, size = n, replace = TRUE)
        ndata <- object[[1]][xx,]
        p0 <- which(ndata$time1 > s)
        p1 <- which(ndata$time1 <= s & ndata$Stime > s)

        den <- KM(ndata$time1, ndata$event1, s)
        den2 <- sum(KMW(ndata$Stime, ndata$event)[p1])
        kmw1 <- KMW(ndata$Stime, ndata$event)

        for (k in 1: length(t)) {
          p2 <- which(ndata$time1 > s & ndata$time1 <= t[k] & ndata$Stime > t[k])
          p3 <- which(ndata$time1 <= s & ndata$Stime > t[k])
          n.p2 <- length(p2)
          n.p3 <- length(p3)
          res.ci[k, j, 1] <- KM(ndata$time1, ndata$event1, t[k]) / den
          ifelse(n.p2 == 0, res.ci[k, j, 2] <- 0, res.ci[k, j, 2] <- sum(kmw1[p2]) / den)
          res.ci[k, j, 3] <- 1 - res.ci[k, j, 1] - res.ci[k, j, 2]
          if(res.ci[k, j, 2] >  1) res.ci[k, j, 2] <- 1
          if(res.ci[k, j, 3] <  0) res.ci[k, j, 3] <- 0
          if (state1 > 0) {
            ifelse(n.p3 == 0, res.ci[k, j, 4] <- 0, res.ci[k, j, 4] <- sum(kmw1[p3]) / den2)
            res.ci[k, j, 5] <- 1 - res.ci[k, j, 4]
          }
        }
      }

      for (k in 1: length(t)) {
        p00.ci[k,1] <- quantile(res.ci[k,,1], (1 - conf.level) / 2)
        p00.ci[k,2] <- quantile(res.ci[k,,1], 1 - (1 - conf.level) / 2)
        p01.ci[k,1] <- quantile(res.ci[k,,2], (1 - conf.level) / 2)
        p01.ci[k,2] <- quantile(res.ci[k,,2], 1 - (1 - conf.level) / 2)
        p02.ci[k,1] <- quantile(res.ci[k,,3], (1 - conf.level) / 2)
        p02.ci[k,2] <- quantile(res.ci[k,,3], 1 - (1 - conf.level) / 2)
        if (state1 > 0) {
          p11.ci[k,1] <- quantile(res.ci[k,,4], (1 - conf.level) / 2)
          p11.ci[k,2] <- quantile(res.ci[k,,4], 1 - (1 - conf.level) / 2)
          p12.ci[k,1] <- quantile(res.ci[k,,5], (1 - conf.level) / 2)
          p12.ci[k,2] <- quantile(res.ci[k,,5], 1 - (1 - conf.level) / 2)
        }
      }
    }


    if(conf == TRUE){
      ci <- cbind(p00.ci, p01.ci, p02.ci, p11.ci, p12.ci)
      ci <- data.frame(ci)
      names(ci) <- c("p00.li.ci", "p00.ls.ci", "p01.li.ci", "p01.ls.ci", "p02.li.ci", "p02.ls.ci", "p11.li.ci", "p11.ls.ci", "p12.li.ci", "p12.ls.ci")
    }

    if(conf==TRUE)  result <- list(est=resu, CI=ci, conf.level=conf.level, s=s, t=t, conf=conf)
    else  result <- list(est=resu,  s=s, t=t, conf=conf)

    class(result) <- c("tpLIDA")
    return(invisible(result))
  }

