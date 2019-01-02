tpLM <-
  function(object, s, conf = FALSE, n.boot = 199, conf.level = 0.95,
           conf.type = "log", cluster = FALSE, ncores = NULL)
  {
    if (missing(object))
      stop("Argument 'object' is missing, with no default")
    #  if (!inherits(object, "survIDM")) stop("'object' must be of class 'survIDM'")
    if (missing(s))
      s <- 0

    t1 <- object[[1]]$time1[object[[1]]$event1 == 1]
    t2 <- object[[1]]$Stime[object[[1]]$event == 1]
    t <- c(s, t1, t2)
    t <- t[t >= s]
    t <- sort(unique(t))

    if (s > max(object[[1]]$time1)) stop("The value of 's' is too large")
    if (s < 0) stop("'s' must be nonnegative")
    if (length(s) > 1) stop("Length of 's' must be 1")
    if (is.null(conf.type)) conf.type <- "bootstrap"
    if(conf.type == "linear") conf.type <- "plain"

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

    kmw0 <- KMW(object[[1]]$time1[p0], object[[1]]$event1[p0])
    kmw1 <- KMW(object[[1]]$Stime[p0], object[[1]]$event[p0])
    if(state1 > 0) kmw2 <- KMW(object[[1]]$Stime[p1], object[[1]]$event[p1])

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


    if (conf==TRUE & conf.type == "bootstrap") {
      res.ci <- array(NA, dim=c(length(t), n.boot, 5))
      n <- dim(object[[1]])[1]

      for (j in 1.:n.boot){
        xx <- sample.int(n, size = n, replace = TRUE)
        ndata <- object[[1]][xx,]
        p0 <- which(ndata$time1 > s)
        p1 <- which(ndata$time1 <= s & ndata$Stime > s)

        kmw0 <- KMW(ndata$time1[p0], ndata$event1[p0])
        kmw1 <- KMW(ndata$Stime[p0], ndata$event[p0])
        if(state1 > 0) kmw2 <- KMW(ndata$Stime[p1], ndata$event[p1])

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

    if (conf == TRUE & (conf.type == "plain" | conf.type == "log" | conf.type == "log-log")) {
      p01.res <- matrix(NA, nrow=length(t), ncol=n.boot)
      se.p01 <- rep(NA, length(t))
      n <- dim(object[[1]])[1]

      for (j in 1:n.boot){
        xx <- sample.int(n, size = n, replace = TRUE)
        ndata <- object[[1]][xx,]
        p0 <- which(ndata$time1 > s)
        kmw0 <- KMW(ndata$time1[p0], ndata$event1[p0])
        kmw1 <- KMW(ndata$Stime[p0], ndata$event[p0])

        for (k in 1: length(t)) {
          q0 <- which(ndata$time1[p0] <= t[k])
          q1 <- which(ndata$Stime[p0] <= t[k])
          nq0 <- length(q0)
          nq1 <- length(q1)

          ifelse(nq0 == 0, a <- 0, a <- sum(kmw0[q0]))
          ifelse(nq1 == 0, b <- 0, b <- sum(kmw1[q1]))
          p01.res[k, j] <- a - b
        }
      }
      for (k in 1: length(t)) { se.p01[k] <- sd(p01.res[k,])  }
      #print(se.p01)

      p0 <- which(object[[1]]$time1 > s)
      p1 <- which(object[[1]]$time1 <= s & object[[1]]$Stime > s)
      fit00 <- survfit(Surv(object[[1]]$time1[p0], object[[1]]$event1[p0]) ~ 1, conf.int = conf.level, conf.type=conf.type)
      fit0 <- survfit(Surv(object[[1]]$Stime[p0], object[[1]]$event[p0]) ~ 1, conf.int = conf.level, conf.type=conf.type)
      se.surv0 <- summary(fit0, times=t, extend=T)$std.err
      S0 <- summary(fit0, times=t, extend=T)$surv
      if (state1 > 0){
        fit1 <- survfit(Surv(object[[1]]$Stime[p1], object[[1]]$event[p1]) ~ 1, conf.int = conf.level, conf.type=conf.type)
        se.surv1 <- summary(fit1, times=t, extend=T)$std.err
        S1 <- summary(fit1, times=t, extend=T)$surv
      }

      p00.ci[,1] <- summary(fit00, times=t, extend=T)$lower
      p00.ci[,2] <- summary(fit00, times=t, extend=T)$upper
      if(state1 > 0){
        p11.ci[,1] <- summary(fit1, times=t, extend=T)$lower
        p11.ci[,2] <- summary(fit1, times=t, extend=T)$upper
      }

      if (conf.type == "plain" | conf.type == "linear") {
        A0 <- qnorm(0.5+0.5*conf.level)*se.surv0
        if (state1 > 0) A1 <- qnorm(0.5+0.5*conf.level)*se.surv1
        B <- qnorm(0.5+0.5*conf.level)*se.p01
      }

      if (conf.type == "log") {
        A0 <- qnorm(0.5+0.5*conf.level)*se.surv0/(1-S0)
        for (i in 1:length(S0)) { if (S0[i] == 1) A0[i] <- 0 }
        if (state1 > 0){ A1 <- qnorm(0.5+0.5*conf.level)*se.surv1/(1-S1)
        for (i in 1:length(S1)) { if (S1[i] == 1) A1[i] <- 0 }
        }
        B <- qnorm(0.5+0.5*conf.level)*se.p01/p01
        for (i in 1:length(p01)) { if (p01[i] == 0) B[i] <- 0 }
      }

      if (conf.type == "log-log") {
        A0 <- qnorm(0.5+0.5*conf.level)*se.surv0/((1-S0)*log(1-S0))
        for (i in 1:length(S0)) { if (S0[i] == 1) A0[i] <- 0 }
        if (state1 > 0) {A1 <- qnorm(0.5+0.5*conf.level)*se.surv1/((1-S1)*log(1-S1))
        for (i in 1:length(S1)) { if (S1[i] == 1) A1[i] <- 0 }
        }
        B <- qnorm(0.5+0.5*conf.level)*se.p01/(p01*log(p01))
        for (i in 1:length(p01)) { if (p01[i] == 0) B[i] <- 0 }
      }

      if (conf.type == "plain" | conf.type == "linear") {
        p02.ci[,1] <- 1 - summary(fit0, times=t, extend=T)$surv - A0
        p02.ci[,2] <- 1 - summary(fit0, times=t, extend=T)$surv + A0
        if(state1 > 0){
          p12.ci[,1] <- 1 - summary(fit1, times=t, extend=T)$surv - A1
          p12.ci[,2] <- 1 - summary(fit1, times=t, extend=T)$surv + A1
        }
        p01.ci[,1] <- p01 - B
        p01.ci[,2] <- p01 + B
      }

      if (conf.type == "log") {
        p02.ci[,1] <- (1 - summary(fit0, times=t, extend=T)$surv)*exp(-A0)
        p02.ci[,2] <- (1 - summary(fit0, times=t, extend=T)$surv)*exp(A0)
        if(state1 > 0){
          p12.ci[,1] <- (1 - summary(fit1, times=t, extend=T)$surv)*exp(-A1)
          p12.ci[,2] <- (1 - summary(fit1, times=t, extend=T)$surv)*exp(A1)
        }
        p01.ci[,1] <- p01*exp(-B)
        p01.ci[,2] <- p01*exp(B)
      }

      if (conf.type == "log-log") {
        p02.ci[,1] <- (1 - summary(fit0, times=t, extend=T)$surv)^exp(-A0)
        p02.ci[,2] <- (1 - summary(fit0, times=t, extend=T)$surv)^exp(A0)
        if(state1 > 0){
          p12.ci[,1] <- (1 - summary(fit1, times=t, extend=T)$surv)^exp(-A1)
          p12.ci[,2] <- (1 - summary(fit1, times=t, extend=T)$surv)^exp(A1)
        }
        p01.ci[,1] <- p01^exp(-B)
        p01.ci[,2] <- p01^exp(B)
      }
    }

    if(conf == TRUE){
      ci <- cbind(p00.ci, p01.ci, p02.ci, p11.ci, p12.ci)
      ci <- data.frame(ci)
      names(ci) <- c("p00.li.ci", "p00.ls.ci", "p01.li.ci", "p01.ls.ci", "p02.li.ci", "p02.ls.ci", "p11.li.ci", "p11.ls.ci", "p12.li.ci", "p12.ls.ci")
    }

    if(conf==TRUE)  result <- list(est=resu, CI=ci, conf.level=conf.level, s=s, t=t, conf=conf)
    else  result <- list(est=resu,  s=s, t=t, conf=conf)

    class(result) <- c("LM")
    return(invisible(result))
  }

