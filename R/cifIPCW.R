cifIPCW <-
  function (object, t, z.name, z.value, bw = "dpik", window = "gaussian",
            method.weights = "NW", conf = FALSE, n.boot = 199, conf.level = 0.95,
            cluster = FALSE, ncores = NULL)
  {

    if (missing(object))
      stop("Argument 'object' is missing, with no default")
   # if (!inherits(object, "survIDM")) stop("'object' must be of class 'survIDM'")
    obj <- object[[1]]
    if(missing(t)) t <- sort(object[[1]]$Stime)
    t <- sort(unique(t))

    t <- c(0,t) # esto es lo nuevo de Luís!!! para que llegue a 0


    covar <- which(names(object[[1]]) == z.name)
    if (missing(z.value))
      z.value <- mean(object[[1]][, covar])
    if (bw == "dpik") {
      lbd2 <- dpik(x = object[[1]][, covar])
    }
    else if (bw == "np") {
      options(np.messages = FALSE)
      lbd2 <- npudensbw(dat = object[[1]][, covar])$bw
    }
    else {
      lbd2 <- bw
    }
    if (!is.numeric(bw) & !(bw %in% c("dpik", "np"))) {
      stop("Argument 'bw' have to be 'dpik', 'np' or a numeric.")
    }
    n <- dim(object[[1]])[1]
    delta4 <- rep(1, n)
    lenc <- dim(object[[1]])[2]
    ifelse(lenc > 4, cov <- 1, cov <- 0)
    #ntimes <- lenc%/%2


    if(length(z.value) > 1){
      stop("Argumment 'z' cannot have length greater than 1")
    }

    n.auxi <- length(t)
    res <- rep(0, n.auxi)

        be <- rep(0, n)
        for (k in 1:n) {
		#be[k] <- Beran(object[[1]]$Stime, 1 - object[[1]]$event,
                #         object[[1]][, covar], delta4, z.value, object[[1]]$Stime[k],
                #         kernel = window, bw = lbd2)
		be[k] <- Beran(object[[1]]$time1, 1 - object[[1]]$event1,
                          object[[1]][, covar], delta4, z.value, object[[1]]$time1[k],
                          kernel = window, bw = lbd2)
        }
        ifelse(method.weights == "NW", w1 <- NWW(object[[1]][, covar], z.value, kernel = window, bw = lbd2),
               w1 <- LLW(object[[1]][, covar], bw = lbd2, t1 = z.value))

        for (j in 1:n.auxi) {
		p1 <- which(object[[1]]$time1 <= t[j] & object[[1]]$time1 < object[[1]]$Stime)
		ifelse(any(be[p1] == 0), res[j] <- NA, res[j] <- sum(w1[p1]/be[p1]))

		if (res[j] > 1)  res[j] <- 1
		if (res[j] < 0)  res[j] <- 0
          }
        resu <- data.frame(cbind(t, res))
        names(resu) <- c("t", "estimate")


    ii <- duplicated(resu$estimate)
    t <- t[!ii]

    n.auxi <- length(t)
    res.li <- rep(0, n.auxi)
    res.ls <- rep(0, n.auxi)


    if (conf == TRUE) {
      simplebootIPCW <- function(object, n, n.auxi) {
        k <- 1
        res.ci <- matrix(0, nrow = n.auxi, ncol = k)
        xx <- sample.int(n, size = n, replace = TRUE)
        ndata <- object[[1]][xx, ]
        obj <- ndata
        ifelse(is.numeric(bw), lbd2 <- bw, lbd2 <- dpik(x = ndata[, covar]))

	be <- rep(0, n)
            for (i in 1:n) {
		#be[i] <- Beran(ndata$Stime, 1 - ndata$event,
                #             ndata[, covar], delta4, z.value, ndata$Stime[i],
                #             kernel = window, bw = lbd2)
		be[i] <- Beran(ndata$time1, 1 - ndata$event1,
                              ndata[, covar], delta4, z.value, ndata$time1[i],
                              kernel = window, bw = lbd2)
            }
            ifelse(method.weights == "NW", w1 <- NWW(ndata[, covar], z.value, kernel = window, bw = lbd2),
                   w1 <- LLW(ndata[, covar], bw = lbd2, t1 = z.value))

            for (j in 1:n.auxi) {
                p1 <- which(ndata$time1 <= t[j] & ndata$time1 < ndata$Stime)
                ifelse(any(be[p1] == 0), res.ci[j, k] <- NA, res.ci[j, k] <- sum(w1[p1]/be[p1]))

                if (res.ci[j, k] > 1) res.ci[j, k] <- 1
                if (res.ci[j, k] < 0) res.ci[j, k] <- 0
				  }
        return(res.ci)
      }

      if (isTRUE(cluster)) {
        if (is.null(ncores)) {
          num_cores <- detectCores() - 1
        }
        else {
          num_cores <- ncores
        }
        registerDoParallel(cores = num_cores)
        on.exit(stopImplicitCluster())
        suppressMessages(res.ci <- foreach(i = 1:n.boot, .combine = cbind) %dorng% simplebootIPCW(object, n, n.auxi))
      }
      else {
        suppressMessages(res.ci <- foreach(i = 1:n.boot, .combine = cbind) %do% simplebootIPCW(object, n, n.auxi))
      }
      for (k in 1:n.auxi) {
        res.li[k] <- quantile(res.ci[k, ], (1 - conf.level)/2, na.rm = TRUE)
        res.ls[k] <- quantile(res.ci[k, ], 1 - (1 - conf.level)/2, na.rm = TRUE)
      }
        resu <- data.frame(cbind(t, res[!ii], res.li, res.ls))
        names(resu) <- c("t", "estimate", paste("lower ",conf.level*100,"% CI", sep=""), paste("upper ",conf.level*100,"% CI", sep=""))
    }

    callp <- paste("CIF(t|", z.name, "=", z.value, ")", sep = "")

    if (conf == FALSE) {
      result <- list(est = resu, estimate = res[!ii], z.name = z.name,
                     z.value = z.value, t = t, bw = bw, window = window,
                     method.weights = method.weights, conf = conf, lbd = lbd2,
                     callp = callp)
    }

    if (conf == TRUE) {
      result <- list(est = resu, estimate = res[!ii], LCI = res.li,
                     UCI = res.ls, conf.level = conf.level, z.name = z.name,
                     z.value = z.value, bw = bw, window = window, method.weights = method.weights,
                     t = t, conf = conf, lbd = lbd2, callp = callp)

    }

    class(result) <- c("IPCW", "cif")
    return(invisible(result))
  }

