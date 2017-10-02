

tpLMAJ <- function(object, s, conf = FALSE, conf.level = 0.95, conf.type = "log")
{
  if (missing(object))
    stop("Argument 'object' is missing, with no default")

 # if (!inherits(object, "survIDM")) stop("'object' must be of class 'survIDM'")
  if (missing(s))
    s <- 0
  data <- object[[1]]

  St1 <- which(data$time1 > s)
  St2 <- which(data$time1 <= s & data$Stime > s)
  data1 <- data[St1, ]
  data2 <- data[St2, ]
  dataS1 <- with(data1, survIDM(time1, event1, Stime, event))
  dataS2 <- with(data2, survIDM(time1, event1, Stime, event))

  resS1 <- tpAJ.aux(object=dataS1, s = 0, conf = conf, conf.level = conf.level,
                    conf.type = conf.type)
  resS2 <- tpAJ.aux(object=dataS2, s = s+0.000001, conf = conf,
                    conf.level = conf.level, conf.type = conf.type)



  if(conf == FALSE) {
    resS.probs <- resS2$probs
    resS.probs[1,] <- resS1$probs
    resS.all.probs <- joindata(resS1, resS2)
    }

  if(conf == TRUE) {
    resS.probs <- resS2$probs
    resS.probs[1:3,] <- resS1$probs
    resS.all.probs <- joindata(resS1, resS2, conf=TRUE)
    }

  newtimes <- sort(unique(c(resS1$times, resS2$times)))

  if(conf == TRUE) {

    # results:
    res <- list(
      # states information:
      s = s, t = resS1$t, states = resS1$states, ns = resS1$ns, tr.states = resS1$tr.states,
      conf.type = conf.type,
      # event times:
      times = newtimes,
      probs = resS.probs, all.probs = resS.all.probs,
      # posible transitions:
      p.trans = resS1$p.trans, conf = conf)

    #est
    suppressWarnings(aux <- matrix(res$all.probs[,1,], ncol = 5, nrow = length(res$times)))
    aux <- data.frame(t = res$times, aux)
    names(aux) <- c("t", "p00", "p01", "p02", "p11", "p12")

    #ci
    auxci <- res$all.probs[,2:3,]
    suppressWarnings(auxci <- data.frame(matrix(auxci, ncol = 10, nrow = length(res$times))))
    names(auxci) <- c("p00.li.ci", "p00.ls.ci", "p01.li.ci", "p01.ls.ci",
                      "p02.li.ci", "p02.ls.ci", "p11.li.ci", "p11.ls.ci",
                      "p12.li.ci", "p12.ls.ci")




    res <- list(est = aux,  CI = auxci, conf.level = conf.level,
                s = res$s, t = res$times, conf = conf, conf.type = res$conf.type)




  } #end if
  else{
    # results:
    res <- list(
      # states information:
      s = s, t = resS1$t, states = resS1$states, ns = resS1$ns, tr.states = resS1$tr.states,
      conf.type = conf.type,
      # event times:
      times = newtimes,
      # occupation or transition probabilities:
      probs = resS.probs, all.probs = resS.all.probs,
      #posible transitions:
      p.trans = resS1$p.trans, conf = conf)

    suppressWarnings(aux <- matrix(res$all.probs, ncol = 5,
                                   nrow = length(res$times)))
    aux <- data.frame(t = res$times, aux)
    names(aux) <- c("t", "p00", "p01", "p02", "p11", "p12")


    res <- list(est = aux,  s = res$s, t = res$times, conf = conf,
                conf.type = res$conf.type)

  }

 # res$call <- match.call()   #Change this?
  class(res) = "tpLMAJ"
  res
}


















# auxiliar function (the of of the package does not work)
# ---------------------------------------------------------

tpAJ.aux <- function(object, s, conf = FALSE, conf.level = 0.95, conf.type = "log")
{
  if (missing(object))
    stop("Argument 'object' is missing, with no default")

  if (!inherits(object, "survIDM")) stop("'object' must be of class 'survIDM'")
  if (missing(s))
    s <- 0
  #	data <- object[[1]]
  data <- object
  t <- max(data$Stime)
  n <- length(data[, 1])

  mint <- min(data$time1[data$event1 == 1])

  if (s < 0) stop("'s' must be nonnegative")
  if (s > 0 & s <= mint) {stop("Argument 's' must be 0 or greater than min(data$time1)=", mint)}

  #Choose the method for the confidence interval
  conf.type2 <- c("linear", "log", "log-log")
  q <- charmatch(conf.type, conf.type2, nomatch = 0) 	# return value 0 if there is no match
  if (q == 0) stop("conf.type should be 'linear', 'log' or 'log-log'")

  states <- c("1", "2", "3")			# states
  ns <- length(states)				# number of states
  tr.states <- states[!states == "3"]  	# transient states
  data$start.time <- 0				# start time is set to 0

  #initial probabilities for each initial state
  i.state <- integer(length(data[, 1]))
  for(i in 1:length(data[, 1])){
    if(data$start.time[i] == 0) i.state[i] <- 1
    if(data$start.time[i] == data$time1[i]) i.state[i] <- 2
    if(data$start.time[i] == data$Stime[i]) i.state[i] <- 3
  }

  i.state <- factor(i.state, levels = states, labels = states)
  initial.probs <- prop.table(table(i.state))

  # prepare data set to compute AJ method:
  ds.prep.AJ <- prepare.aj.data(data, states, tr.states)

  # reduces to event times:
  ds.event.AJ <- prepare.aj.event(ds.prep.AJ$dNs, ds.prep.AJ$Ys, ds.prep.AJ$sum_dNs, states, tr.states)
  event.times <- as.numeric(as.character(rownames(ds.event.AJ$dNs)))

  # Estimates for the AJ estimator:
  AJ.est <- fun.AJ(ns,states, ds.event.AJ$dNs, ds.event.AJ$Ys, ds.event.AJ$sum_dNs, s, t,  event.times, initial.probs)

  if(conf == TRUE) {
    # Variance of the AJ estimator:
    variances <- var.AJ(ns, states, AJ.est$dNs.id_tr, AJ.est$Ys.id_tr, AJ.est$sum_dNs.id_tr, AJ.est$TP.AJs, AJ.est$all.I.dA, tr.states)

    # Confidence Interval for the AJ estimator:
    ci <- ci.AJ(s,t, conf.level, conf.type, AJ.est$dNs.id_tr, AJ.est$TP.AJs, variances$cov.AJs, AJ.est$e.times.id_tr)

    # results:
    res <- list(
      # states information:
      s = s, t = AJ.est$t, states = states, ns = ns, tr.states = tr.states,
      conf.type = conf.type,
      # event times:
      times = AJ.est$e.times.id_tr,
      # occupation or transition probabilities:
      #est=AJ.est$probs, #all.est=AJ.est$TP.AJs,
      # confidence intervals:
      probs = ci$CI, all.probs = ci$all.CI,
      # posible transitions:
      p.trans = AJ.est$p.trans, conf = conf, AJ.est=AJ.est)
  } #end if
  else{
    # results:
    res <- list(
      # states information:
      s = s, t = AJ.est$t, states = states, ns = ns, tr.states = tr.states,
      conf.type = conf.type,
      # event times:
      times = AJ.est$e.times.id_tr,
      # occupation or transition probabilities:
      probs = round(AJ.est$probs, 7), all.probs = round(AJ.est$all.est, 7),
      #posible transitions:
      p.trans = AJ.est$p.trans, conf = conf)
  }

  res$call <- match.call()
  class(res) = "tpAJ"
  res
}

######





# function to joint the results of {p11, p12, p13} and {p22, p23}
# -----------------------------------------------------------------
joindata <- function(x, y, conf = FALSE)
{

  if(conf == FALSE){
    n1 <- dim(x$all.probs)[1]
    n2 <- dim(y$all.probs)[1]
    x1 <- x$all.probs[,1,]
    x2 <- y$all.probs[,1,]
    q <- match(rownames(x1), rownames(x2))
    nmatch <- sum(!is.na(q))
    m <- n1+n2-nmatch

    rn <- unique(c(rownames(y$all.probs),rownames(x$all.probs)))
    rn <- as.character(sort(as.numeric(rn)))
    cn <- colnames(y$all.probs[,1,])
    resdata <- array(NA, dim=c(m,1,5), dimnames=list(rn,"probs",cn))

    q1 <- match(rownames(x1), rownames(resdata))
    q2 <- match(rownames(x2), rownames(resdata))
    resdata[q1,1,1:3] <- x$all.probs[,1,1:3]
    resdata[q2,1,4:5] <- y$all.probs[,1,4:5]

    for (k in 2:m){
      for(j in 1:5){
        if (is.na(resdata[k,1,j])) resdata[k,1,j] <- resdata[k-1,1,j]
      }
    }
  }

  if(conf == TRUE){
    n1 <- dim(x$all.probs)[1]
    n2 <- dim(y$all.probs)[1]
    x1 <- x$all.probs[,1,]
    x2 <- y$all.probs[,1,]
    q <- match(rownames(x1), rownames(x2))
    nmatch <- sum(!is.na(q))
    m <- n1+n2-nmatch

    rn <- unique(c(rownames(y$all.probs),rownames(x$all.probs)))
    rn <- as.character(sort(as.numeric(rn)))
    cn <- colnames(y$all.probs[,1,])
    c3 <- c("probs","LCI","UCI","Var")
    resdata <- array(NA, dim=c(m,4,5), dimnames=list(rn,c3,cn))

    q1 <- match(rownames(x1), rownames(resdata))
    q2 <- match(rownames(x2), rownames(resdata))
    resdata[q1,1,1:3] <- x$all.probs[,1,1:3]
    resdata[q2,1,4:5] <- y$all.probs[,1,4:5]
    resdata[q1,2,1:3] <- x$all.probs[,2,1:3]
    resdata[q2,2,4:5] <- y$all.probs[,2,4:5]
    resdata[q1,3,1:3] <- x$all.probs[,3,1:3]
    resdata[q2,3,4:5] <- y$all.probs[,3,4:5]
    resdata[q1,4,1:3] <- x$all.probs[,4,1:3]
    resdata[q2,4,4:5] <- y$all.probs[,4,4:5]

    for (k in 2:m){
      for(i in 1:4){
        for(j in 1:5){
          if (is.na(resdata[k,i,j])) resdata[k,i,j] <- resdata[k-1,i,j]
        }
      }
    }
  }
  resdata
}

####
