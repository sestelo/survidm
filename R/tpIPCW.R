tpIPCW <-
  function (object, s, t, z.name, z.value, bw = "dpik", window = "gaussian",
            method.weights = "NW", conf = FALSE, n.boot = 199, conf.level = 0.95,
            cluster = FALSE, ncores = NULL)
  {
    if (missing(object))
      stop("Argument 'object' is missing, with no default")
    # if (!inherits(object, "survIDM")) stop("'object' must be of class 'survIDM'")
    if (missing(s)) s <- 0
    obj <- object[[1]]

    ptimes <- which(obj$event1 == 1 | obj$event == 1)

    if (missing(t)) t <- obj$Stime[ptimes]
    if (any(t <= 0)) stop("The values of 't' must be positive")
    if (s > max(object[[1]]$time1)) stop("The value of 's' is too large")
    if (s < 0) stop("'s' must be nonnegative")
    if (length(s) > 1) stop("Length of 's' must be 1")

    t <- t[t >= s]
    t <- sort(unique(t))
    if(length(t) == 0) stop("Invalid values for 't'.")

    covar <- which(names(object[[1]]) == z.name)
    if (missing(z.value))
      z.value <- mean(object[[1]][, covar])
    n <- dim(obj)[1]
    delta4 <- rep(1, n)
    lenc <- dim(obj)[2]
    ifelse(lenc > 4, cov <- 1, cov <- 0)

    if(length(z.value) > 1){ stop("Argumment 'z' cannot have length greater than 1") }

    p0 <- which(obj$time1 > s)
    p1 <- which(obj$time1 <= s & obj$Stime > s)
    p1new <- which(obj$time1 <= s & obj$Stime > s & obj$event == 1)
    state1 <- length(p1new)
    n0 <- length(p0)
    n1 <- length(p1)

    if (bw == "dpik") {
      lbd2_0 <- dpik(x = obj[p0, covar])
      if(state1 > 0) lbd2_1 <- dpik(x = obj[p1, covar])
    }
    else if (bw == "np") {
      options(np.messages = FALSE)
      lbd2_0 <- npudensbw(dat = obj[p0, covar])$bw
      if(state1 > 0) lbd2_1 <- npudensbw(dat = obj[p1, covar])$bw
    }
    else {
      lbd2 <- bw
    }
    if (!is.numeric(bw) & !(bw %in% c("dpik", "np"))) {
      stop("Argument 'bw' have to be 'dpik', 'np' or a numeric.")
    }

    p00 <- rep(NA, length(t))
    p01 <- rep(NA, length(t))
    p02 <- rep(NA, length(t))
    p11 <- rep(NA, length(t))
    p12 <- rep(NA, length(t))

    be0 <- rep(0, n0)
    be1 <- rep(0, n0)
    be2 <- rep(0, n1)

    for (k in 1:n0) {
      be1[k] <- Beran(obj$Stime[p0], 1 - obj$event[p0], obj[p0, covar], delta4[p0], z.value, obj$Stime[p0][k],
                      kernel = window, bw = lbd2_0)
      be0[k] <- Beran(obj$time1[p0], 1 - obj$event1[p0], obj[p0, covar], delta4[p0], z.value, obj$time1[p0][k],
                      kernel = window, bw = lbd2_0)
    }
    if(state1 > 0) {
      for (k in 1:n1) {
        be2[k] <- Beran(obj$Stime[p1], 1 - obj$event[p1], obj[p1, covar], delta4[p1], z.value, obj$Stime[p1][k],
                        kernel = window, bw = lbd2_1)
      }
    }

    ifelse(method.weights == "NW", w0 <- NWW(obj[p0, covar], z.value, kernel = window, bw = lbd2_0),
           w0 <- LLW(obj[p0, covar], bw = lbd2_0, t1 = z.value))
    if(state1 > 0) {
      ifelse(method.weights == "NW", w1 <- NWW(obj[p1, covar], z.value, kernel = window, bw = lbd2_1),
             w1 <- LLW(obj[p1, covar], bw = lbd2_1, t1 = z.value))}

    for (k in 1: length(t)) {
      q0 <- which(obj$time1[p0] <= t[k] & obj$event1[p0] == 1)
      q1 <- which(obj$Stime[p0] <= t[k] & obj$event[p0] == 1)
      q2 <- which(obj$Stime[p1] <= t[k] & obj$event[p1] == 1)
      nq0 <- length(q0)
      nq1 <- length(q1)
      nq2 <- length(q2)
      ifelse(nq0 == 0, p00[k] <- 1, p00[k] <- 1 - sum(w0[q0]/be0[q0]))
      if (!is.na(p00[k])) { if(p00[k] <  0) p00[k] <- 0}
      ifelse(nq1 == 0, p02[k] <- 0, p02[k] <- sum(w0[q1]/be1[q1]))
      if (!is.na(p02[k])) { if(p02[k] > 1) p02[k] <- 1
      if(p02[k] < 0) p02[k] <- 0}
      p01[k] <- 1 - p00[k] - p02[k]
      if (!is.na(p01[k])) { if(p01[k] < 0) p01[k] <- 0
      if(p01[k] > 1) p01[k] <- 1}
      if (state1 > 0) {
        ifelse(nq2 == 0, p12[k] <- 0, p12[k] <- sum(w1[q2]/be2[q2]))
        if (!is.na(p12[k])) { if(p12[k] < 0) p12[k] <- 0
        if(p12[k] > 1) p12[k] <- 1}
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


    if (conf == TRUE) {
      res.ci <- array(NA, dim=c(length(t), n.boot, 5))
      n <- dim(obj)[1]

      for (j in 1:n.boot){
        xx <- sample.int(n, size = n, replace = TRUE)
        ndata <- object[[1]][xx,]
        p0 <- which(ndata$time1 > s)
        p1 <- which(ndata$time1 <= s & ndata$Stime > s)
        n0 <- length(p0)
        n1 <- length(p1)

        be0 <- rep(0, n0)
        be1 <- rep(0, n0)
        be2 <- rep(0, n1)
        for (k in 1:n0) {
          be1[k] <- Beran(ndata$Stime[p0], 1 - ndata$event[p0], ndata[p0, covar], delta4[p0], z.value, ndata$Stime[p0][k],
                          kernel = window, bw = lbd2_0)
          be0[k] <- Beran(ndata$time1[p0], 1 - ndata$event1[p0], ndata[p0, covar], delta4[p0], z.value, ndata$time1[p0][k],
                          kernel = window, bw = lbd2_0)
        }

        if(state1 > 0){
          for (k in 1:n1) {
            be2[k] <- Beran(ndata$Stime[p1], 1 - ndata$event[p1], ndata[p1, covar], delta4[p1], z.value, ndata$Stime[p1][k],
                            kernel = window, bw = lbd2_1)
          }
        }

        ifelse(method.weights == "NW", w0 <- NWW(ndata[p0, covar], z.value, kernel = window, bw = lbd2_0),
               w0 <- LLW(ndata[p0, covar], bw = lbd2_0, t1 = z.value))

        if(state1 > 0){
          ifelse(method.weights == "NW", w1 <- NWW(ndata[p1, covar], z.value, kernel = window, bw = lbd2_1),
                 w1 <- LLW(ndata[p1, covar], bw = lbd2_1, t1 = z.value))
        }

        for (k in 1: length(t)) {
          q0 <- which(ndata$time1[p0] <= t[k] & ndata$event1[p0] == 1)
          q1 <- which(ndata$Stime[p0] <= t[k] & ndata$event[p0] == 1)
          q2 <- which(ndata$Stime[p1] <= t[k] & ndata$event[p1] == 1)
          nq0 <- length(q0)
          nq1 <- length(q1)
          nq2 <- length(q2)

          ifelse(nq0 == 0, res.ci[k, j, 1] <- 1, res.ci[k, j, 1] <- 1 - sum(w0[q0]/be0[q0]))
          if (!is.na(res.ci[k, j, 1])) { if(res.ci[k, j, 1] <  0) res.ci[k, j, 1] <- 0}
          ifelse(nq1 == 0, res.ci[k, j, 3] <- 0, res.ci[k, j, 3] <- sum(w0[q1]/be1[q1]))
          if (!is.na(res.ci[k, j, 3])) { if(res.ci[k, j, 3] <  0) res.ci[k, j, 3] <- 0
          if(res.ci[k, j, 3] >  1) res.ci[k, j, 3] <- 1}
          res.ci[k, j, 2] <- 1 - res.ci[k, j, 1] - res.ci[k, j, 3]
          if (!is.na(res.ci[k, j, 2])) { if(res.ci[k, j, 2] <  0) res.ci[k, j, 2] <- 0
          if(res.ci[k, j, 2] >  1) res.ci[k, j, 2] <- 1}
          if (state1 > 0) {
            ifelse(nq2 == 0, res.ci[k, j, 5] <- 0, res.ci[k, j, 5] <- sum(w1[q2]/be2[q2]))
            if (!is.na(res.ci[k, j, 5])) { if(res.ci[k, j, 5] <  0) res.ci[k, j, 5] <- 0
            if(res.ci[k, j, 5] >  1) res.ci[k, j, 5] <- 1}
            res.ci[k, j, 4] <- 1 - res.ci[k, j, 5]
          }
        }
      }

      for (k in 1: length(t)) {
        p00.ci[k,1] <- quantile(res.ci[k,,1], (1 - conf.level) / 2, na.rm = TRUE)
        p00.ci[k,2] <- quantile(res.ci[k,,1], 1 - (1 - conf.level) / 2, na.rm = TRUE)
        p01.ci[k,1] <- quantile(res.ci[k,,2], (1 - conf.level) / 2, na.rm = TRUE)
        p01.ci[k,2] <- quantile(res.ci[k,,2], 1 - (1 - conf.level) / 2, na.rm = TRUE)
        p02.ci[k,1] <- quantile(res.ci[k,,3], (1 - conf.level) / 2, na.rm = TRUE)
        p02.ci[k,2] <- quantile(res.ci[k,,3], 1 - (1 - conf.level) / 2, na.rm = TRUE)
        if (state1 > 0) {
          p11.ci[k,1] <- quantile(res.ci[k,,4], (1 - conf.level) / 2, na.rm = TRUE)
          p11.ci[k,2] <- quantile(res.ci[k,,4], 1 - (1 - conf.level) / 2, na.rm = TRUE)
          p12.ci[k,1] <- quantile(res.ci[k,,5], (1 - conf.level) / 2, na.rm = TRUE)
          p12.ci[k,2] <- quantile(res.ci[k,,5], 1 - (1 - conf.level) / 2, na.rm = TRUE)
        }
      }
    }

    # callp <- paste("pij(s,t|", z.name, "=", z.value, ")", sep = "")


    if(conf == TRUE){
      ci <- cbind(p00.ci, p01.ci, p02.ci, p11.ci, p12.ci)
      ci <- data.frame(ci)
      names(ci) <- c("p00.li.ci", "p00.ls.ci", "p01.li.ci", "p01.ls.ci", "p02.li.ci", "p02.ls.ci", "p11.li.ci", "p11.ls.ci", "p12.li.ci", "p12.ls.ci")
    }

    if (state1 == 0) lbd2_1 <- NA
    if (conf == FALSE) {
      result <- list(est = resu, z.name = z.name, z.value = z.value, s = s,
                     t = t, bw = bw, window = window,
                     method.weights = method.weights, conf = conf, lbd = c(lbd2_0,lbd2_1))
    }

    if (conf == TRUE) {
      result <- list(est = resu, CI = ci, conf.level = conf.level, z.name = z.name,
                     z.value = z.value, bw = bw, window = window, method.weights = method.weights,
                     s = s, t = t, conf = conf, lbd = c(lbd2_0,lbd2_1))

    }

    class(result) <- c("IPCW", "tp")
    return(invisible(result))
  }

