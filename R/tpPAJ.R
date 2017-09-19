tpPAJ <- function(object,
                  s,
                  conf = FALSE,
                  conf.level = 0.95)
{
  if (missing(object))
    stop("Argument 'object' is missing, with no default")

  #  if (!inherits(object, "survIDM")) stop("'object' must be of class 'survIDM'")
  if (missing(s))
    s <- 0
  data <- object[[1]]

  TPobj <- with(data, TPmsm::survTP(time1, event1, Stime, event))
  res <-
    TPmsm::transPAJ(
      object = TPobj,
      s = s,
      t = Inf,
      conf = conf,
      conf.level = conf.level
    )
  #res <- list()

  probs <- matrix(0, 3, 3)
  colnames(probs) <- c("1", "2", "3")
  rownames(probs) <- c("1", "2", "3")
  probs[1, ] <- res$est[length(res$time), 1:3]
  probs[2, ] <- c(0, res$est[length(res$time), 4:5])
  probs[3, ] <- c(0, 0, 1)

  cn <- colnames(res$est)
  c3 <- "probs"
  rn <- as.character(res$time)
  m <- length(rn)
  ptrans <- c("1 1", "1 2", "1 3", "2 2", "2 3")

  if (conf == TRUE) {
    c3 <- c("probs", "LCI", "UCI", "Var")
    all.probs <- array(NA, dim = c(m, 4, 5), dimnames = list(rn, c3, cn))
    all.probs[, 1, 1:5] <- res$est
    all.probs[, 2, 1:5] <- res$inf
    all.probs[, 3, 1:5] <- res$sup
    # results:
    res <- list(
      # states information:
      s = s,
      t = max(res$time),
      states = c("1", "2", "3"),
      ns = 3,
      tr.states = c("1", "2"),
      conf.type = "bootstrap",
      # event times:
      times = res$time,
      probs = probs,
      all.probs = all.probs,
      # posible transitions:
      p.trans = ptrans,
      conf = conf
    )

    #est
    suppressWarnings(aux <- matrix(
      res$all.probs[, 1, ],
      ncol = 5,
      nrow = length(res$times)
    ))
    aux <- data.frame(t = res$times, aux)
    names(aux) <- c("t", "p00", "p01", "p02", "p11", "p12")

    #ci
    auxci <- res$all.probs[, 2:3, ]
    suppressWarnings(auxci <- data.frame(matrix(
      auxci, ncol = 10,
      nrow = length(res$times)
    )))
    names(auxci) <-
      c(
        "p00.li.ci",
        "p00.ls.ci",
        "p01.li.ci",
        "p01.ls.ci",
        "p02.li.ci",
        "p02.ls.ci",
        "p11.li.ci",
        "p11.ls.ci",
        "p12.li.ci",
        "p12.ls.ci"
      )




    res <- list(
      est = aux,
      CI = auxci,
      conf.level = conf.level,
      s = res$s,
      t = res$times,
      conf = conf,
      conf.type = res$conf.type
    )





  } #end if
  else{
    all.probs <- array(NA, dim = c(m, 1, 5), dimnames = list(rn, c3, cn))
    all.probs[, 1, 1:5] <- res$est
    # results:
    res <- list(
      # states information:
      s = s,
      t = max(res$time),
      states = c("1", "2", "3"),
      ns = 3,
      tr.states = c("1", "2"),
      conf.type = "",
      # event times:
      times = res$time,
      # occupation or transition probabilities:
      probs = probs,
      all.probs = all.probs,
      #posible transitions:
      p.trans = ptrans,
      conf = conf
    )

    suppressWarnings(aux <- matrix(
      res$all.probs,
      ncol = 5,
      nrow = length(res$times)
    ))
    aux <- data.frame(t = res$times, aux)
    names(aux) <- c("t", "p00", "p01", "p02", "p11", "p12")


    res <- list(
      est = aux,
      s = res$s,
      t = res$times,
      conf = conf,
      conf.type = res$conf.type
    )



  }
  return(res)
}
