

nevents <- function (dataidm, state.names=NULL)
{

  tmat <- matrix(NA, 3, 3)
  tmat[1, 2:3] <- 1:2
  tmat[2, 3] <- 3
  if (missing(state.names))
    state.names <- c("healthy", "illness", "death")
  else {
    if (length(state.names) != 3)
      stop("incorrect length of 'state.names' argument")
  }
  dimnames(tmat) <- list(from = state.names, to = state.names)

  tmat2 <- matrix(NA, 3, 3)
  n00 <- length(which(dataidm$event1 == 0))
  n01 <- length(which(dataidm$event1 == 1 & dataidm$time1 < dataidm$Stime))
  n02 <- length(which(dataidm$event1 == 1 & dataidm$event == 1 & dataidm$time1 == dataidm$Stime))
  n12 <- length(which(dataidm$event1 == 1 & dataidm$time1 < dataidm$Stime & dataidm$event == 1))

  n0 <- length(which(dataidm$event1 == 1 & dataidm$event == 0 & dataidm$time1 == dataidm$Stime))

  if(n0>0) warning(n0, " observations with ", "Stime = time1 from state ", state.names[2], " to state ",state.names[3]  )

  colnames(tmat2) <- c(state.names)
  rownames(tmat2) <- state.names

  #tmat2[1, ] <- c(n00,     n01, n02)
  tmat2[1, ] <- c(n00,     n01-length(which(dataidm$event1 == 1 & dataidm$time1==0)), n02) #retirar casos em que time1=0 e event1=1
  tmat2[2, ] <- c(0,        n01-n12, n12)
  tmat2[3, ] <- c(0,        0,   n02+n12)

  counts <- tmat2
  class(counts) <- "table"
  return(counts)
}
