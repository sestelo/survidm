
markov.test <- function (formula, s, nm.method = "LM", data)
{

  if (missing(formula))
    stop("Argument 'formula' is missing with no default")
  if (class(formula) !="formula") stop("The argument 'formula' must be of class 'formula'")
  if (missing(data))
    stop("Argument 'data' is missing with no default")
  if (missing(s)) {
    s <- quantile(data$time1, 0.25)
    cat("Argument 's' is missing. First quartile of the sojourn time in the initial state has been considered for the graphical test","\n")
  }
  if(s==0) stop("Markov assumption is not relevant for the estimation of occupation probabilities (s==0)")

  if(nm.method == "AJ" | nm.method == "PAJ") stop("The chosen method is nor non-Markov")
  # if(nm.method != "LM" & nm.method != "LIDA") stop("'LM' and 'LIDA' are the only non-Markov methods available")

  term1 <- formula[[2]]
  var1 <- all.vars(formula[[2]])
  term2 <- formula[[3]]
  var2 <- all.vars(formula[[3]])
  if(length(var1)!=4) stop("Incorrect variables in 'formula'")

  pos1 <- match(all.vars(formula)[1:4], names(data))
  pos2 <- match(all.vars(formula)[-c(1:4)], names(data))

  mydata <- data[,c(pos1,pos2)]
  names(mydata)[1:4] <- c("time1", "event1", "Stime", "event")

  p01 <- which(mydata$Stime > mydata$time1)
  mydata12 <- mydata[p01, ]

  if(length(var2) == 0) {
    fmla2 <- as.formula(paste("Surv(time1, Stime, event) ~ time1"))
    fit12 <- survival::coxph(fmla2, data = mydata12)
  }

  ncov <- length(var2) #number of covariates
  covars <- var2 #vector with names of covariates

  if(length(var2) > 0){
    if(is.na(match("time1", var2))){
      fmla2 <- as.formula(paste("Surv(time1, Stime, event)~", paste(formula[3], collapse = "+")))
      fmla2 <- as.formula(paste("Surv(time1, Stime, event)~",
                                paste(c(var2,"time1"), collapse = "+")))
    }
    else   fmla2 <- as.formula(paste("Surv(time1, Stime, event)~", paste(formula[3], collapse = "+")))
    fit12 <- survival::coxph(fmla2, data = mydata12)
  }

  #results from AJ and the markov free estimator of the transition probability
  #without covariates
  resAJ <- tprob(survIDM(time1, event1, Stime, event) ~ 1, s = s,
                 method = "AJ", conf = FALSE, data = mydata12)
  resNM <- tprob(survIDM(time1, event1, Stime, event) ~ 1, s = s,
                 method = nm.method, conf = TRUE, data = mydata12)

  times <- resNM$t

  #if(nm.method == "LM" | nm.method == "LIDA") {
  est.aj.01 <- resAJ$est$p01
  est.nm.01 <- resNM$est$p01

  est.aj.02 <- resAJ$est$p02
  est.nm.02 <- resNM$est$p02

  est.aj.12 <- resAJ$est$p12
  est.nm.12 <- resNM$est$p12

  est.nm.12LCI <- resNM$CI[,9]
  est.nm.12UCI <- resNM$CI[,10]
  #}

  if(nm.method == "PLM" | nm.method == "LMAJ" | nm.method == "PLMAJ") {

    time.aj <- resAJ$t
    pos <- match(time.aj, times)
    est.aj.01 <- rep(0, length(times))
    est.aj.01[pos] <- resAJ$est$p01
    est.aj.02 <- rep(0, length(times))
    est.aj.02[pos] <- resAJ$est$p02
    est.aj.12 <- rep(0, length(times))
    est.aj.12[pos] <- resAJ$est$p12
    for (k in 2:length(times)){
      if(est.aj.01[k] == 0) est.aj.01[k] <- est.aj.01[k-1]
      if(est.aj.02[k] == 0) est.aj.02[k] <- est.aj.02[k-1]
      if(est.aj.12[k] == 0) est.aj.12[k] <- est.aj.12[k-1]
    }
  }

  TPestimates <- data.frame(aj01=est.aj.01, aj02=est.aj.02, aj12=est.aj.12,
                            nm01=est.nm.01, nm02= est.nm.02, nm12=est.nm.12,
                            nm12LCI=est.nm.12LCI, nm12UCI=est.nm.12UCI, times=times)

  res <-
    list(
      cox.markov.test = fit12,
      TPestimates = TPestimates,
      nm.method = nm.method,
      s = s
    )

  class(res) <- c("markov", "survIDM")
  return(res)
}
