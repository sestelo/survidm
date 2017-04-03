tpPLDM <-
  function(object, s, conf = FALSE, n.boot = 199, conf.level = 0.95,
           cluster = FALSE, ncores = NULL)
  {
    if (missing(object))
      stop("Argument 'object' is missing, with no default")
   # if (!inherits(object, "survIDM")) stop("'object' must be of class 'survIDM'")
    if (missing(s))
      s <- 0

    ptimes <- which(object[[1]]$event1 == 1 | object[[1]]$event == 1)

    t1 <- object[[1]]$time1[object[[1]]$event1 == 1]
    t2 <- object[[1]]$Stime[object[[1]]$event == 1]
    t <- c(s, t1, t2)
    t <- t[t >= s]
    t <- sort(unique(t))

    if (any(t < 0)) stop("The values of 't' must be positive")
    if (s > max(object[[1]]$time1)) stop("The value of 's' is too large")
    if (s < 0) stop("'s' must be nonnegative")
    if (length(s) > 1) stop("Length of 's' must be 1")

    n <- length(object[[1]]$Stime)
    if(length(t) == 0) stop("Invalid values for 't'.")

    p0 <- which(object[[1]]$time1 > s)
    p1 <- which(object[[1]]$time1 <= s & object[[1]]$Stime > s)
    p1new <- which(object[[1]]$time1 <= s & object[[1]]$Stime > s & object[[1]]$event == 1)
    state1 <- length(p1new)

    p00 <- rep(NA, length(t))
    p01 <- rep(NA, length(t))
    p02 <- rep(NA, length(t))
    p11 <- rep(NA, length(t))
    p12 <- rep(NA, length(t))

    suppressWarnings(kmw0 <- PKMW(object[[1]]$time1[p0], object[[1]]$event1[p0]))
    suppressWarnings(kmw1 <- PKMW(object[[1]]$Stime[p0], object[[1]]$event[p0]))
    if(state1 > 0) suppressWarnings(kmw2 <- PKMW(object[[1]]$Stime[p1], object[[1]]$event[p1]))

    for (k in 1: length(t)) {
      q0 <- which(object[[1]]$time1[p0] <= t[k])
      q1 <- which(object[[1]]$Stime[p0] <= t[k])
      q2 <- which(object[[1]]$Stime[p1] <= t[k])
      nq0 <- length(q0)
      nq1 <- length(q1)
      nq2 <- length(q2)
      ifelse(nq0 == 0, p00[k] <- 1, p00[k] <- 1 - sum(kmw0[q0]))
      ifelse(nq1 == 0, p02[k] <- 0, p02[k] <- sum(kmw1[q1]))
      p01[k] <- 1 - p00[k] - p02[k]
      if(p01[k] <  0) p01[k] <- 0
      if (state1 > 0) {
        ifelse(nq2 == 0, p12[k] <- 0, p12[k] <- sum(kmw2[q2]))
        p11[k] <- 1 - p12[k]
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


    if (conf==TRUE) {
      res.ci <- array(NA, dim=c(length(t), n.boot, 5))
      n <- dim(object[[1]])[1]

      for (j in 1.:n.boot){
        xx <- sample.int(n, size = n, replace = TRUE)
        ndata <- object[[1]][xx,]
        p0 <- which(ndata$time1 > s)
        p1 <- which(ndata$time1 <= s & ndata$Stime > s)

        suppressWarnings(kmw0 <- PKMW(ndata$time1[p0], ndata$event1[p0]))
        suppressWarnings(kmw1 <- PKMW(ndata$Stime[p0], ndata$event[p0]))
        if(state1 > 0) suppressWarnings(kmw2 <- PKMW(ndata$Stime[p1], ndata$event[p1]))

        for (k in 1: length(t)) {
          q0 <- which(ndata$time1[p0] <= t[k])
          q1 <- which(ndata$Stime[p0] <= t[k])
          q2 <- which(ndata$Stime[p1] <= t[k])
          nq0 <- length(q0)
          nq1 <- length(q1)
          nq2 <- length(q2)

          ifelse(nq0 == 0, res.ci[k, j, 1] <- 1, res.ci[k, j, 1] <- 1 - sum(kmw0[q0]))
          ifelse(nq1 == 0, res.ci[k, j, 3] <- 0, res.ci[k, j, 3] <- sum(kmw1[q1]))
          res.ci[k, j, 2] <- 1 - res.ci[k, j, 1] - res.ci[k, j, 3]
          if(res.ci[k, j, 2] <  0) res.ci[k, j, 2] <- 0
          if (state1 > 0) {
            ifelse(nq2 == 0, res.ci[k, j, 5] <- 0, res.ci[k, j, 5] <- sum(kmw2[q2]))
            res.ci[k, j, 4] <- 1 - res.ci[k, j, 5]
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

    class(result) <- c("PLDM")
    return(invisible(result))
  }



