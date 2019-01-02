tpPLMAJ <- function(object, s, conf = FALSE, conf.level = 0.95)
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
  dataS1 <- with(data1, TPmsm::survTP(time1, event1, Stime, event))
  dataS2 <- with(data2, TPmsm::survTP(time1, event1, Stime, event))

  resS1 <- TPmsm::transPAJ(object=dataS1, s = 0, t=Inf,
                           conf = conf, conf.level = conf.level)
  resS2 <- TPmsm::transPAJ(object=dataS2, s = s, t=Inf,
                           conf = conf, conf.level = conf.level)

  probs <- matrix(0, 3, 3)
  colnames(probs) <- c("1","2","3")
  rownames(probs) <- c("1","2","3")
  probs[1,] <- resS1$est[length(resS1$time),1:3]
  probs[2,] <- c(0,resS2$est[length(resS2$time),4:5])
  probs[3,] <- c(0,0,1)

  newtimes <- sort(unique(c(resS1$time,resS2$time)))
  newtimes <- newtimes[newtimes >= s]  #to include s

  resu <- joindataTP(resS1, resS2, s=s, conf = conf)

  #print(dim(resu))
  #print(resu)


  if(conf == TRUE) {
    # results:
    res <- list(
      # states information:
      s = s, t = max(newtimes),
      states = c("1", "2", "3"),
      ns = 3,
      tr.states = c("1", "2"),
      conf.type = "bootstrap",
      # event times:
      times = newtimes,
      probs = probs,
      all.probs = resu,
      # posible transitions:
      p.trans = c("1 1", "1 2", "1 3", "2 2", "2 3"),
      conf = conf)


    #est
    suppressWarnings(aux <- matrix(res$all.probs[,1,], ncol = 5,
                                   nrow = length(res$times)))
    aux <- data.frame(t = res$times, aux)
    names(aux) <- c("t", "p00", "p01", "p02", "p11", "p12")

    #ci
    auxci <- res$all.probs[,2:3,]
    suppressWarnings(auxci <- data.frame(matrix(auxci, ncol = 10,
                                                nrow = length(res$times))))
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
      s = s, t = max(newtimes),
      states = c("1", "2", "3"),
      ns = 3,
      tr.states = c("1", "2"),
      conf.type = "bootstrap",
      # event times:
      times = newtimes,
      # occupation or transition probabilities:
      probs = probs, all.probs = resu,
      #posible transitions:
      p.trans = c("1 1", "1 2", "1 3", "2 2", "2 3"), conf = conf)

    suppressWarnings(aux <- matrix(res$all.probs, ncol = 5,
                                   nrow = length(res$times)))
    aux <- data.frame(t = res$times, aux)
    names(aux) <- c("t", "p00", "p01", "p02", "p11", "p12")


    res <- list(est = aux,  s = res$s, t = res$times, conf = conf,
                conf.type = res$conf.type)



  }

  res$call <- match.call()
  class(res) = "tpPLMAJ"
  res
}





joindataTP <- function(x, y, conf = FALSE, s = NULL)
{
  newtimes <- sort(unique(c(x$time,y$time,s)))
  newtimes <- newtimes[newtimes >= s]
  nt <- length(newtimes)
  n1 <- dim(x$est)[1]
  n2 <- dim(y$est)[1]
  x1 <- x$time
  x1 <- x1[x1 >= s]
  x2 <- y$time
  x2 <- x2[x2 >= s]
  m <- nt

  rn <- newtimes
  rn <- as.character(rn)
  cn <- c("1 1","1 2","1 3","2 2","2 3")

  if(conf == FALSE){
    resdata <- array(NA, dim=c(m,1,5), dimnames=list(rn,"probs",cn))

    q1 <- match(x1, newtimes)
    q2 <- match(x2, newtimes)
    resdata[q1,1,1:3] <- x$est[x$time > s,1:3]
    resdata[q2,1,4:5] <- y$est[y$time >= s,4:5]

    resdata[1,1,] <- c(1,0,0,1,0)

    for (k in 2:m){
      for(j in 1:5){
        if (is.na(resdata[k,1,j])) resdata[k,1,j] <- resdata[k-1,1,j]
      }
    }
  }

  if(conf == TRUE){
    c3 <- c("probs","LCI","UCI","Var")
    resdata <- array(NA, dim=c(m,4,5), dimnames=list(rn,c3,cn))

    q1 <- match(x1, newtimes)
    q2 <- match(x2, newtimes)
    resdata[q1,1,1:3] <- x$est[x$time >= s,1:3]
    resdata[q2,1,4:5] <- y$est[y$time >= s,4:5]
    resdata[q1,2,1:3] <- x$inf[x$time >= s,1:3]
    resdata[q2,2,4:5] <- y$inf[y$time >= s,4:5]
    resdata[q1,3,1:3] <- x$sup[x$time >= s,1:3]
    resdata[q2,3,4:5] <- y$sup[y$time >= s,4:5]
    resdata[1,1,] <- c(1,0,0,1,0)
    resdata[1,2,] <- c(1,0,0,1,0)
    resdata[1,3,] <- c(1,0,0,1,0)

    for (k in 2:m){
      for(i in 1:4){
        for(j in 1:5){
          if (is.na(resdata[k,i,j])) resdata[k,i,j] <- resdata[k-1,i,j]
        }
      }
    }
  }

  #print(rn)
  p <- which(newtimes <=s)
  #print(resdata)
  #print(p)
  #resdata <- rbind(c(365, 1.0000000, 0.0000000, 0.0000000, 1.0000000, 0.0000000),resdata[-p,,])

  resdata
}


