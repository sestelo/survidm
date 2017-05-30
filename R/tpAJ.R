tpAJ <- function(object, s, conf = FALSE, conf.level = 0.95, conf.type = "log")
{
  if (missing(object))
    stop("Argument 'object' is missing, with no default")
 # if (!inherits(object, "survIDM")) stop("'object' must be of class 'survIDM'")
  if (missing(s))
    s <- 0
  data <- object[[1]]
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
      p.trans = AJ.est$p.trans, conf = conf)


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
      s = s, t = AJ.est$t, states = states, ns = ns, tr.states = tr.states,
      conf.type = conf.type,
      # event times:
      times = AJ.est$e.times.id_tr,
      # occupation or transition probabilities:
      probs = round(AJ.est$probs, 7), all.probs = round(AJ.est$all.est, 7),
      #posible transitions:
      p.trans = AJ.est$p.trans, conf = conf)


    suppressWarnings(aux <- matrix(res$all.probs, ncol = 5,
                                   nrow = length(res$times)))
    aux <- data.frame(t = res$times, aux)
    names(aux) <- c("t", "p00", "p01", "p02", "p11", "p12")


    res <- list(est = aux,  s = res$s, t = res$times, conf = conf,
                   conf.type = res$conf.type)
  }


  #res$call <- match.call()
  class(res) = "tpAJ"
  res
}

#############################################################################################
fun.AJ <- function(ns, states, dNs, Ys, sum_dNs, s, t, event.times, initial.probs)
{
  id_tr <- which(s < event.times & event.times <= t) ## location of those (s, t]
  l.id_tr <- length(id_tr)

  # to use in variance function:
  e.times.id_tr <- event.times[id_tr]

  ## Occupation Probability matrix
  OP <- matrix(NA, nrow = l.id_tr, ncol = length(states))
  rownames(OP) <- rownames(dNs)[id_tr]; colnames(OP) <- states
  all.dA <- all.I.dA <- array(dim = c(length(states), length(states), l.id_tr),
                              dimnames=list(rows=states,cols=states, time=rownames(dNs)[id_tr]))

  ## Transition Probability matrix
  cum.prod <- diag(ns)
  rownames(cum.prod) <- states

  TP.AJs <- array(dim=c(ns, ns, l.id_tr), dimnames = list(rows = states, cols = states,
                                                          dim = rownames(dNs)[id_tr]))

  dNs.id_tr <- dNs[id_tr,]
  Ys.id_tr <- Ys[id_tr,]
  sum_dNs.id_tr <- sum_dNs[id_tr,]

  for(i in 1:l.id_tr){
    I.dA <- diag(ns) ## creates trans matrix for current time
    dA <- matrix(0, nrow = ns, ncol = ns)
    colnames(I.dA) <- rownames(I.dA) <- colnames(dA) <- rownames(dA) <- states
    i_tr <- which(dNs.id_tr[i, , drop = FALSE] > 0)  	## indicator for kind of transition at time i
    dNs.event.names <- colnames(dNs.id_tr)[i_tr] 	## gets names of transitions (ie:  dN##)
    split_dNs.event <- strsplit(dNs.event.names, " ")    	## splits title of dN##
    st.start <- sapply(split_dNs.event, function(x) x[2])
    st.end <- sapply(split_dNs.event, function(x) x[3])  ## start & stop states as character strings
    i_tr.s <- matrix(as.character(c(st.start, st.end)), ncol = 2)
    i_tr.s2 <- matrix(as.character(c(st.start, st.start)), ncol = 2)

    dA[i_tr.s] <- dNs.id_tr[i, i_tr] / Ys.id_tr[i, paste("Y", st.start)]
    if (length(i_tr) == 1) { dA[st.start, st.start] <- -dNs.id_tr[i, i_tr]/Ys.id_tr[i, paste("Y", st.start)]}
    else {      dA[i_tr.s2] <- -rowSums(dA[st.start, ])     }

    I.dA <- I.dA + dA 				## I+dA (transition) matrix
    all.dA[, , i] <- dA     			## stores all dA matrices
    all.I.dA[, , i] <- I.dA 			## array for storing all tran matrices
    cum.prod <- cum.prod %*% I.dA
    TP.AJs[,,i] <- cum.prod
  } ## end of loop i

  if (s == 0){
    OP[i, ] <- initial.probs%*%TP.AJs[, , i] 	## state occupation probabilities
    op <- OP[i, ]
    p.transitions <- c("1 1", "1 2", "1 3")

    all.est <- array(0, dim = c(l.id_tr, 1, length(p.transitions)),
                     dimnames = list(rows = rownames(dNs)[id_tr], cols = "probs", trans = p.transitions))

    for(j in 1:length(p.transitions))   {
      idx <- unlist(strsplit(p.transitions[j], " "))
      all.est[ , 1, j] <- TP.AJs[idx[1], idx[2], ]
    }

    res <- list(probs = op, all.est = all.est, TP.AJs = TP.AJs, all.I.dA = all.I.dA, dNs.id_tr = dNs.id_tr,
                Ys.id_tr = Ys.id_tr, sum_dNs.id_tr = sum_dNs.id_tr, e.times.id_tr = e.times.id_tr,
                p.trans = p.transitions, t = t)
    return(res)
  } #end if
  else{
    p.transitions <- c("1 1", "1 2", "1 3", "2 2", "2 3")
    all.est <- array(0, dim = c(l.id_tr, 1, length(p.transitions)),
                     dimnames = list(rows = rownames(dNs)[id_tr], cols = "probs", trans = p.transitions))

    for(j in 1:length(p.transitions))	{
      idx <- unlist(strsplit(p.transitions[j]," "))
      all.est[, 1, j] <- TP.AJs[idx[1], idx[2], ]
    }

    res <- list(probs = cum.prod, all.est = all.est, TP.AJs = TP.AJs, all.I.dA = all.I.dA, dNs.id_tr = dNs.id_tr,
                Ys.id_tr = Ys.id_tr, sum_dNs.id_tr = sum_dNs.id_tr, e.times.id_tr = e.times.id_tr,
                p.trans = p.transitions, t = t)
    return(res)
  }
}

#############################################################################################
prepare.aj.data <- function(data, states, tr.states){
  ttime <- c(data$time1, data$Stime)
  times <- sort(unique(ttime)) # T_ik* indicates the instants of the k transitions

  # matrix with posible transitions:
  mat_w <- matrix(FALSE, nrow = 3, ncol = 3)
  mat_w[1, 2:3] <- TRUE
  mat_w[2, 3] <- TRUE

  # matrix with all transitions (including censoring transitions):
  mat_c <- matrix(TRUE, nrow = 3, ncol = 3)
  mat_c[2, 1] <- FALSE
  mat_c[3, 1:3] <- FALSE

  colnames(mat_w) <- rownames(mat_w) <- states
  colnames(mat_c) <- rownames(mat_c) <- states

  # output states contempling censoring (tr 11, tr 22)
  output.states.c <- lapply(1:dim(mat_c)[2], function(i){rownames(mat_c)[mat_c[i, ] == TRUE]})

  # into states whitout censoring (tr 11, tr 22)
  into.states <- lapply(1:dim(mat_w)[2], function(i){colnames(mat_w)[mat_w[, i] == TRUE]})

  # possible transitions (including censoring)
  to <- c(output.states.c[[1]], output.states.c[[2]])
  from <- c(rep(tr.states[[1]], length(output.states.c[[1]])), rep(tr.states[[2]], length(output.states.c[[2]])))
  transitions <- paste("tr", from, to)

  # risk sets for each non-absorbing state (names)
  ys <- paste("Y", tr.states)

  n = length(data$time1)
  m = length(times)

  # number of patientes with time1 = time k-transition (if event1=0 & event=0 --> 1 cens)
  g11 <- integer(m)
  for(i in 1:m){
    ss11 <- subset(data,data$event1==0 & data$event==0)
    g11[i] <- sum(ss11$time1==times[i])
  }

  # number of patients with time1<time k-transition (for transition 1->3)
  g13 <- integer(m)
  for (i in 1:m){
    ss13 <- subset(data, data$event1 == 1 & data$event == 1 & data$time1 == data$Stime)
    g13[i] <- sum(ss13$time1 == times[i])
  }

  # number of patients with time1<time k-transition (for transition 1->2)
  g12 <- integer(m)
  for (i in 1:m){
    ss12 <- subset(data, data$time1 < data$Stime)
    g12[i] <- sum(ss12$time1 == times[i])
  }

  # number of patients with time1<time k-transition (for transition 2->2)
  g22 <- integer(m)
  for (i in 1:m){
    ss22 <- subset(data, data$event1 == 1 & data$event == 0 & data$time1 < data$Stime)
    g22[i] <- sum(ss22$Stime == times[i])
  }

  # number of patients with time1<time k-transition (for transition 2->3)
  g23 <- integer(m)
  for (i in 1:m){
    ss23 <- subset(data, data$event1 == 1 & data$event == 1 & data$time1 < data$Stime)
    g23[i] <- sum(ss23$Stime == times[i])
  }

  ##  matrix of number of transitions:
  dNs <- cbind(g11, g12, g13, g22, g23)
  colnames(dNs) <- transitions
  rownames(dNs) <- times

  ## matrix of total number of transitions from each non-absorbing state:
  sum_n1 <- rowSums(dNs[, 2:3])
  sum_n2 <- dNs[, 5]
  sum_dNs <- cbind(sum_n1, sum_n2)

  data$start.time <- 0

  initial.state <- integer(length(data[, 1]))
  for(i in 1:length(data[, 1])){
    if(data$start.time[i] == 0) initial.state[i] <- 1
    else initial.state[i] <- 2
  }
  initial.state <- factor(initial.state, levels = tr.states, labels = tr.states)

  initial_risk_set <- table(initial.state)

  i <- 1
  n <- initial_risk_set[i]
  name <- paste("y", i)

  if(length(into.states[[i]]) > 0) in.st <- paste("tr", into.states[[i]], i)
  else in.st <- NULL

  if(length(output.states.c[[i]]) > 0) out.st<-paste("tr", i, output.states.c[[i]])
  else out.st <- NULL

  Y1 <- c(n, n + cumsum(rowSums(dNs[,in.st,drop=FALSE]))-cumsum(rowSums(dNs[,out.st,drop=FALSE])))

  i <- 2
  n <- initial_risk_set[i]
  name <- paste("y", i)

  if(length(into.states[[i]]) > 0) in.st <- paste("tr", into.states[[i]], i)
  else in.st <- NULL

  if(length(output.states.c[[i]]) > 0) out.st <- paste("tr", i, output.states.c[[i]])
  else out.st <- NULL

  Y2 <- c(n, n + cumsum(rowSums(dNs[,in.st,drop=FALSE]))-cumsum(rowSums(dNs[,out.st,drop=FALSE])))

  Ys <- cbind(Y1[1:m], Y2[1:m])

  rownames(dNs) <- rownames(sum_dNs) <- rownames(Ys) <- times
  colnames(dNs) <- transitions
  colnames(Ys) <- ys
  colnames(sum_dNs) <- paste("from", tr.states)

  split_dNs <- strsplit(colnames(dNs), "  ") ## string splits names
  split_Ys <- strsplit(colnames(Ys), " ")  ## string split names of Ys
  split_sum_dNs <- strsplit(colnames(sum_dNs), " ") ## string splits names of dNs

  ## censored columns in dNs, Ys and sum_dNs counting matrix, needed for D-S est
  cens.dNs.id <- which(sapply(split_dNs, function(x) x[1]=="tr 1 1"|x[1]=="tr 2 2"))
  cens.Ys.id <- which(sapply(split_Ys, function(x) x[2]%in%tr.states))
  cens.sum_dNs.id <- which(sapply(split_sum_dNs, function(x) x[2]%in%tr.states))

  K <- vector(length=nrow(dNs))
  dN.cens <- rowSums(dNs[, cens.dNs.id, drop=FALSE])
  Y.cens <- rowSums(Ys[, cens.Ys.id, drop=FALSE]) ## those at risk of being censored
  N.Y.cens <- ifelse(dN.cens/Y.cens=="NaN", 0, dN.cens/Y.cens)
  colnames(N.Y.cens) <- NULL
  H.t <- cumsum(N.Y.cens) ## calculating the hazard
  k <- exp(-H.t)
  K <- c(1, k[-length(k)])

  # weighted for right censoring
  dNs.w <- dNs/K  ## D-S dNs
  Ys.w <- Ys/K  ## D-S Ys
  sum_dNs.w <- sum_dNs/K

  res <- list(dNs = dNs.w, Ys = Ys.w, sum_dNs = sum_dNs.w, transitions = transitions)
  return(invisible(res))
}

#############################################################################################
prepare.aj.event <- function(dNs, Ys, sum_dNs, states, tr.states) {

  split_dNs <- strsplit(colnames(dNs), "  ") ## string splits names
  split_Ys <- strsplit(colnames(Ys), " ")  ## string split names of Ys
  split_sum_dNs <- strsplit(colnames(sum_dNs), " ") ## string splits names of dNs

  ## reducing dNs & Ys to just event times & noncens states
  ##  looks at noncensored columns
  event.dNs.id <- which(sapply(split_dNs, function(x) x[1]=="tr 1 2" | x[1]=="tr 1 3" | x[1]=="tr 2 3"))
  ## identifies times where transitions occur
  event.row.id <- which(apply(dNs[, event.dNs.id, drop=FALSE], 1, function(x) any(x>0)))
  dNs.event <- dNs[event.row.id, event.dNs.id, drop=FALSE] ## reduces dNs

  tr.states.event <- names(which(sapply(tr.states, function(x) length(x) > 0)))
  event.Ys.id <- which(sapply(split_Ys, function(x) x[2] %in% tr.states.event))
  Ys.event <- Ys[event.row.id, event.Ys.id, drop = FALSE] ## reduces Ys

  event.sum_dNs.id <- which(sapply(split_sum_dNs, function(x) x[2]%in%states))
  sum_dNs.event <- sum_dNs[event.row.id, event.sum_dNs.id, drop=FALSE]

  ans <- list(dNs=dNs.event, Ys=Ys.event, sum_dNs=sum_dNs.event)
  return(ans)
}

#############################################################################################
var.AJ <- function(ns, states, dNs.id_tr, Ys.id_tr, sum_dNs.id_tr, TP.AJs, all.I.dA, tr.states){

  transition.names <- paste(rep(states, ns), rep(states, each = ns))
  cov.dA <- cov.AJs <- array(0, dim = c(ns^2, ns^2, nrow(dNs.id_tr)),
                             dimnames = list(rows = transition.names, cols = transition.names,
                                             time = rownames(dNs.id_tr)))
  results.matrix <- array(0, dim = c(ns^2, ns^2)) 	## matrix to keep the results to build a cov.AJs matrix
  colnames(results.matrix) <- rownames(results.matrix) <- transition.names

  # identity matrices for Kronecker products
  bl.Id <- diag(ns^2)
  Id <- diag(ns)

  vp <- matrix(0, ns, ns)

  for (i in 1:nrow(dNs.id_tr)) { 	## loop through times

    ## VARIANCE OF A-J ( TRANS PROB MATRIX P(s, t) )
    ## loop on the blocks (g) - only needed for transient states
    for (g in tr.states) {
      ## find positioning of g in definition of states
      id_g <- which(states == g)
      ## covs: written componentwise, the recursive formula to calculate the variance of transitions in each block (g)
      covs <- matrix(0, nrow = ns, ncol = ns)
      colnames(covs) <- rownames(covs) <- states

      ## loop in the blocks
      for (j in 1:ns) {

        ## This just fills in upper diagonal matrix - use symmetry to fill in the rest
        for (r in j:ns) {
          state.j <- states[j]
          state.r <- states[r]
          g.Ys <- paste("Y", g)
          g.sum_dNs <- paste("from", g)

          if (Ys.id_tr[i, g.Ys] == 0) {  ## if Y_g = 0 then covariance = 0
            covs[j, r] <- 0
            next
          }

          if (state.j == g & state.r == g) {  ## g==j==r
            covs[j, r] <- (Ys.id_tr[i, g.Ys] - sum_dNs.id_tr[i, g.sum_dNs]) * sum_dNs.id_tr[i, g.sum_dNs] / Ys.id_tr[i, g.Ys]^3  }
          else if (state.j == g & state.r != g) {  ## g!=r
            name <- paste("tr", g, state.r)
            if (!name%in%colnames(dNs.id_tr)) next
            covs[j, r] <- -(Ys.id_tr[i, g.Ys] - sum_dNs.id_tr[i, g.sum_dNs])*dNs.id_tr[i, name] / Ys.id_tr[i, g.Ys]^3 }
          else if (state.j != g & state.r == g) {  ## g!=j
            name <- paste("tr", g, state.j)
            if (!name%in%colnames(dNs.id_tr)) next
            covs[j, r] <- -(Ys.id_tr[i, g.Ys] - sum_dNs.id_tr[i, g.sum_dNs])*dNs.id_tr[i, name]/Ys.id_tr[i, g.Ys]^3 }
          else { ## g!=j and g!=r
            name.r <- paste("tr", g, state.r)
            name.j <- paste("tr", g, state.j)
            if (!(name.j%in%colnames(dNs.id_tr) & name.r%in%colnames(dNs.id_tr))) next
            covs[j, r] <- (ifelse(j==r, 1, 0)*Ys.id_tr[i, g.Ys] - dNs.id_tr[i, name.j])*dNs.id_tr[i, name.r]/Ys.id_tr[i, g.Ys]^3
          } ## end of if/else statements

        } ## end of r loop
      } ## end of j loop

      covs[lower.tri(covs)] <- t(covs)[lower.tri(covs)]

      results.matrix[(seq(1, ns*(ns-1)+1, by=ns) + id_g - 1), (seq(1, ns*(ns-1)+1, by=ns) + id_g - 1)] <- covs

    }## end of g loop

    ## array holding var-cov matrix for I+dA matrix at each time (differential of NA estim)
    cov.dA[, , i] <- results.matrix

    if (i==1) { cov.AJs[, , i] <- bl.Id%*% cov.dA[, , i] %*% bl.Id   }
    else {
      cov.AJs[, , i] <- (t(all.I.dA[, , i]) %x% Id) %*% cov.AJs[, , i-1] %*%((all.I.dA[, , i]) %x% Id) +
        (Id %x% TP.AJs[, , i-1]) %*% cov.dA[, , i]  %*% (Id%x% t(TP.AJs[, , i-1]))
    }
  }
  return(list(TP.AJs = TP.AJs, cov.AJs = cov.AJs))
}

#############################################################################################
ci.AJ <- function(s, t, conf.level, conf.type = "linear", dNs.id_tr, TP.AJs, cov.AJs, e.times.id_tr)
{
  if(s == 0){
    p.transitions <- c("1 1", "1 2", "1 3")
    zsig <- qnorm(conf.level + (1 - conf.level) / 2)

    CI.tp <- array(0, dim = c(nrow(dNs.id_tr), 4, length(p.transitions)),
                   dimnames = list(rows = e.times.id_tr, cols = c("probs", "lower", "upper", "variance"),
                                   trans = p.transitions))

    # different transformations to built CI
    conf.type <- match.arg(conf.type, c("linear", "log", "log-log"))

    for (j in 1:length(p.transitions)) { 		## loop through possible transitions

      idx <- unlist(strsplit(p.transitions[j], " "))
      CI.tp[ , 1, j] <- P <- TP.AJs[idx[1], idx[2] , ]
      CI.tp[ , 4, j] <- var <- cov.AJs[p.transitions[j], p.transitions[j], ]


      switch(conf.type[1],
             "linear" = {	CI.tp[ , 2, j] <- P - zsig * sqrt(var)
             CI.tp[ , 3, j] <- P + zsig * sqrt(var)},
             "log" = {	CI.tp[ , 2, j] <- exp(log(P) - zsig * sqrt(var) / P)
             CI.tp[ , 3, j] <- exp(log(P) + zsig * sqrt(var) / P)},
             "log-log" = {CI.tp[ , 2, j] <- P^(exp(-zsig * (sqrt(var) / (P * log(P)))))
             CI.tp[ , 3, j] <- P^(exp(zsig * (sqrt(var) / (P * log(P)))))})

      CI.tp[ , 2, j] <- pmax(CI.tp[ , 2, j], 0)
      CI.tp[ , 3, j] <- pmin(CI.tp[ , 3, j], 1)
    } 	## end j loop

    CI.t <- matrix(0, nrow = length(p.transitions), ncol = 4)
    colnames(CI.t) <- c("probs", "lower", "upper", "variance")
    rownames(CI.t) <- p.transitions

    for(j in 1:length(p.transitions)){      CI.t[j, ] <- CI.tp[nrow(CI.tp[, , j]), , j]  }

    CI.tp <- round(CI.tp, 7)
    CI.t <-round(CI.t, 7)
  } #end if

  else{
    p.transitions <- c("1 1", "1 2", "1 3", "2 2", "2 3")
    zsig <- qnorm(conf.level + (1 - conf.level) / 2)

    CI.tp <- array(0, dim = c(nrow(dNs.id_tr), 4, length(p.transitions)),							#results for output
                   dimnames = list(rows = e.times.id_tr, cols = c("probs", "lower", "upper", "variance"),
                                   trans = p.transitions))

    # different transformations to built the Confidence Intervals
    conf.type <- match.arg(conf.type, c("linear", "log", "log-log"))

    for (j in 1:length(p.transitions)) { ## loop through possible transitions

      idx <- unlist(strsplit(p.transitions[j], " "))
      CI.tp[ , 1, j] <- P <- TP.AJs[idx[1], idx[2] , ]
      CI.tp[ , 4, j] <- var <- cov.AJs[p.transitions[j], p.transitions[j], ]

      switch(conf.type[1],
             "linear" = {	CI.tp[ , 2, j] <- P - zsig * sqrt(var)
             CI.tp[ , 3, j] <- P + zsig * sqrt(var)},
             "log" = {	CI.tp[ , 2, j] <- exp(log(P) - zsig * sqrt(var) / P)
             CI.tp[ , 3, j] <- exp(log(P) + zsig * sqrt(var) / P)},
             "log-log" = {CI.tp[ , 2, j] <- P^(exp(-zsig * (sqrt(var) / (P * log(P)))))
             CI.tp[ , 3, j] <- P^(exp(zsig * (sqrt(var) / (P * log(P)))))})

      CI.tp[ , 2, j] <- pmax(CI.tp[ , 2, j], 0)
      CI.tp[ , 3, j] <- pmin(CI.tp[ , 3, j], 1)

    } # end of j loop

    CI.t <- matrix(0, nrow = length(p.transitions), ncol = 4)
    colnames(CI.t) <- c("probs", "lower", "upper", "variance")
    rownames(CI.t) <- p.transitions

    for(j in 1:length(p.transitions)){ CI.t[j, ] <- CI.tp[nrow(CI.tp[, , j]), , j]    }

    CI.tp <- round(CI.tp, 7)
    CI.t <- round(CI.t, 7)
  }

  return(list(CI = CI.t, all.CI = CI.tp))
}
