#' Nonparametric estimation of transition probabilities in
#' the illness-death model
#'
#' @description This function is used to obtain nonparametric estimates of
#' the transition probabilities in the illness-death model.
#'
#' @param formula A \code{formula} object, which must have a \code{survIDM}
#' object as the response on the left of the \code{~} operator and, if desired,
#' a term on the right. The term may be a qualitative or quantitative variable.
#' Without covariates, the right hand side should be \code{~ 1}.
#' @param s The first time for obtaining estimates for the transition
#' probabilities. If missing, 0 will be used.
#' @param method The method used to compute the transition probabilities.
#' Possible options are \code{"AJ"}, \code{"LIDA"} \code{"LDM"}, \code{"PLDM"} and
#' \code{"IPCW"}. Defaults to \code{"AJ"}. The \code{"IPCW"} method
#' is recommended to obtain conditional transition probabilities (i.e., with a
#' quantitative term on the right hand side of formula).
#' @param conf Provides pointwise confidence bands. Defaults to \code{FALSE}.
#' @param conf.level Level of confidence. Defaults to 0.95 (corresponding to 95\%).
#' @param conf.type Method to compute the confidence intervals.
#' Transformation applied to compute confidence intervals. Possible choices
#' are \code{"linear"}, \code{"log"}, \code{"log-log"} and \code{"bootstrap"}.
#' Default is \code{"linear"}.
#' @param n.boot The number of bootstrap replicates to compute the variance
#' of the non-Markovian estimator. Default is 199.
#' @param data A data.frame including at least four columns named
#' \code{time1}, \code{event1}, \code{Stime} and \code{event}, which correspond
#' to disease free survival time, disease free survival indicator, time to death
#' or censoring, and death indicator, respectively.
#' @param z.value The value of the covariate on the right hand side of formula
#' at which the transition probabilities are computed. For quantitative
#' covariates, i.e. of class integer and numeric.
#' @param bw A single numeric value to compute a kernel density bandwidth.
#' Use \code{"dpik"} for the \pkg{KernSmooth} package based selector or \code{"np"}
#' for the \code{'npudensbw'} function of the \pkg{np} package.
#' @param window A character string specifying the desired kernel.
#' See details below for possible options. Defaults to \code{"gaussian"}
#' where the gaussian density kernel will be used.
#' @param method.weights A character string specifying the desired weights method.
#' Possible options are \code{"NW"} for the Nadaraya-Watson weights and \code{"LL"}
#' for local linear weights. Defaults to \code{"NW"}.
#' @param cluster A logical value. If \code{TRUE} (default), the bootstrap procedure
#' for the confidence intervals is parallelized. Note that there are
#' cases (e.g., a low number of bootstrap repetitions) that \R will gain in
#' performance through serial computation. \R takes time to distribute tasks
#' across the processors also it will need time for binding them all together
#' later on. Therefore, if the time for distributing and gathering pieces
#' together is greater than the time need for single-thread computing,
#' it does not worth parallelize.
#' @param ncores An integer value specifying the number of cores to be used in
#' the parallelized procedure. If \code{NULL} (default), the number of cores
#' to be used is equal to the number of cores of the machine - 1.
#' @param na.rm A logical value indicating whether NA values should
#' be stripped in the computation.
#
#'
#' @details Possible options for argument window are \code{"gaussian"},
#' \code{"epanechnikov"}, \code{"tricube"}, \code{"boxcar"},
#' \code{"triangular"}, \code{"quartic"} or \code{"cosine"}.
#'
#'
#'
#' @return An object of class \code{"survIDM"} and one of the following
#' five classes: \code{"AJ"}, \code{"LIDA"}, \code{"LMD"}, \code{"PLDM"} and
#' \code{"tpIPCW"}. Objects are implemented as a list with elements:
#'
#' \item{est}{data.frame with estimates of the transition probabilities.}
#' \item{CI}{data.frame with the confidence intervals of the transition probabilities.}
#' \item{conf.level}{Level of confidence.}
#' \item{s}{The first time for obtaining estimates for the transition probabilities.}
#' \item{t}{The time for obtaining the estimates of transition probabilities.}
#' \item{conf}{logical; if \code{FALSE} (default) the pointwise confidence
#' bands are not given.}
#' \item{conf.type}{Type of the confidence interval.}
#' \item{callp}{The expression of the estimated probability.}
#' \item{Nlevels}{The number of levels of the covariate. Provides important
#' information when the covariate at the right hand side of formula
#' is of class factor.}
#' \item{levels}{The levels of the qualitative covariate
#' (if it is of class factor) on the right hand side of formula.}
#' \item{formula}{A formula object.}
#' \item{call}{A call object.}
#'
#'
#' @author Luis Meira-Machado and Marta Sestelo.
#'
#' @references
#' Aalen O. O., Johansen S. (1978) An Empirical Transition Matrix for
#' Nonhomogeneous Markov Chains Based on Censored Observations. Scandinavian
#' Journal of Statistics 5(3), 141--150.
#'
#' Meira-Machado L. F., de Una-Alvarez J. and Cadarso-Suarez C. (2006).
#' Nonparametric estimation of transition probabilities in a non-Markov
#' illness-death model. Lifetime Data Anal 12(3), 325--344.
#'
#' de Una-Alvarez J. and Meira-Machado L. (2015). Nonparametric estimation
#' of transition probabilities in a non-Markov illness-death model:
#' a comparative study. Biometrics 71, 364--375.
#'
#'
#' @examples
#' # Aalen-Johansen
#' # Occupation Probabilities Pj(t)=Pij(0,t)
#'
#' res <- tprob(survIDM(time1, event1, Stime, event) ~ 1, s = 0,
#' method = "AJ", conf = FALSE, data = colonIDM)
#'
#' summary(res, time=365*1:6)
#' plot(res)
#'
#'
#' # Transition Probabilities Pij(t)=Pij(365,t)
#'
#' # LIDA
#' res1 <- tprob(survIDM(time1, event1, Stime, event) ~ 1, s = 365,
#' method = "LIDA", conf = FALSE, data = colonIDM)
#'
#' summary(res1, time=365*1:6)
#' plot(res1)
#' plot(res1, trans="01", ylim=c(0,0.15))
#'
#' # Landmark (LDM)
#' res2 <- tprob(survIDM(time1, event1, Stime, event) ~ 1, s = 365,
#' method = "LDM", conf = FALSE, data = colonIDM)
#'
#' summary(res2, time=365*1:6)
#' plot(res2)
#'
#' # Presmoothed LDM
#' res3 <- tprob(survIDM(time1, event1, Stime, event) ~ 1, s = 365,
#' method = "PLDM", conf = FALSE, data = colonIDM)
#'
#' summary(res3, time=365*1:6)
#' plot(res3)
#'
#'
#'
#' # Conditional transition probabilities
#'
#' #with factor
#' res4 <- tprob(survIDM(time1, event1, Stime, event) ~ factor(sex), s = 365,
#' method = "AJ", conf = FALSE, data = colonIDM)
#'
#' summary(res4, time=365*1:6)
#' plot(res4, trans="02", ylim=c(0,0.5))
#'
#' # with continuous covariate (IPCW)
#' res5 <- tprob(survIDM(time1, event1, Stime, event) ~ age, s = 365,
#' method = "IPCW", z.value = 48, conf = FALSE, data = colonIDM,
#' bw = "dpik", window = "gaussian", method.weights = "NW")
#'
#' summary(res5, time=365*1:6)
#' plot(res5)
#'
#'
#' # Confidence intervals
#' res6 <- tprob(survIDM(time1, event1, Stime, event) ~ 1, s = 365,
#' method = "AJ", conf = TRUE, conf.level = 0.95,
#' conf.type = "log", data = colonIDM)
#'
#' summary(res6, time=365*1:7)
#' plot(res6)





tprob <- function(formula, s, method = "AJ", conf = FALSE, conf.level = 0.95,
                  conf.type = "linear", n.boot = 199, data, z.value, bw = "dpik",
                  window = "gaussian", method.weights = "NW", cluster = FALSE,
                  ncores = NULL, na.rm = TRUE){


  if (missing(formula)) stop("A formula argument is required")
  if (missing(s)) stop("argument 's' is missing, with no default")

  if (!(method %in% c("AJ", "LIDA", "LDM", "PLDM", "IPCW"))){
    stop("Possible methods are 'AJ', 'LIDA', 'LDM', 'PLDM' and 'IPCW'." )
  }



  # formula
  fmla <- eval(formula, parent.frame())
  Terms <- terms(fmla)
  mf <- Terms[[2]]
  object <- with(data = data, eval(mf))
  if (!inherits(object, "survIDM")) stop("Response must be a survIDM object")
  object <- list(data = object) # new since survCS doesn't return a list
  obj_data <- object[[1]]

  X <- Terms[[3]] #covariate
  if(length(attr(terms(formula),"term.labels")) > 1)
    stop("only one covariate is supported")
  Class <- class(with(data = data, eval(Terms[[3]])))
  if (Class != "numeric" & Class != "integer" & Class != "factor")
    stop("covariate must be one of the following classes 'numeric', 'integer' or 'factor'")
  xval <- with(data = data, eval(X))
  lencov <- length(xval)
  lencov2 <- length(attr(terms(formula),"term.labels"))
  if(lencov2 != 0) Xval <- with(data = data, eval(Terms[[3]]))
  if(lencov != dim(obj_data)[1] & lencov2 != 0) stop("length of the covariate does not match")

  #lenc <- dim(object[[1]])[2]
  #ntimes <- lenc%/%2
  #if (length(x) != ntimes-1) {
  #  cat("The number of consecutive event times in 'survCS' is", ntimes, ". The length of 'x' should be", ntimes-1,"\n")
  #    stop("The length of 'x' is not supported for the selected 'object'")
  # }




  # without covariates
  if (length(attr(terms(formula), "term.labels")) == 0) {  #AJ, LIDA, LDM, PLDM without covariate

    if (!(method %in% c("AJ", "LIDA", "LDM", "PLDM"))){
      stop("The model does not include covariates. Possible methods are 'AJ', 'LIDA', 'LDM' and 'PLDM'." )
    }



    # AJ method  (no bootstrap)
    if (method == "AJ"){
      res <- tpAJ(object = object, s = s, conf = conf,
                  conf.level = conf.level, conf.type = conf.type)
      class(res) <- c("AJ", "survIDM")

      add2t <- s
      add2est <- c(s, 1, 0, 0, 1, 0)
      add2CI <- c(1, 1, 0, 0, 0, 0, 1, 1, 0, 0)
      res$t <- c(add2t, res$t)
      res$est <- rbind(add2est, res$est)
      res$CI <- rbind(add2CI, res$CI)

       if (s == 0){
       res$est[, 5:6] <- NA
       res$CI[, 7:10] <- NA
       }

    }


    # LIDA method
    if (method == "LIDA"){

      if(conf == TRUE & conf.type != "bootstrap") {
        warning("This method only allows bootstrap confidence intervals.")
      }

      res <- tpLIDA(object = object, s = s, conf = conf,
                    conf.level = conf.level, n.boot = n.boot,
                    cluster = cluster, ncores = ncores)
      class(res) <- c("LIDA", "survIDM")
    }


    # LDM method
    if (method == "LDM"){
      res <- tpLDM(object = object, s = s, conf = conf,
                   conf.level = conf.level, conf.type = conf.type,
                   n.boot = n.boot, cluster = cluster, ncores = ncores)
      class(res) <- c("LDM", "survIDM")
    }


    # PLDM method
    if (method == "PLDM"){

      if(conf == TRUE & conf.type != "bootstrap") {
        warning("This method only allows bootstrap confidence intervals.")
      }

      res <- tpPLDM(object = object, s = s, conf = conf,
                    conf.level = conf.level, n.boot = n.boot,
                    cluster = cluster, ncores = ncores)
      class(res) <- c("PLDM", "survIDM")
    }

    # in order to have the same output
    x.nlevels <- 1
    levels <- NULL

  } # end methods without covariate







  # numeric or integer covariate
  if(length(attr(terms(formula),"term.labels")) != 0 & (Class == "numeric" | Class == "integer")) {#IPCW

   if(method != "IPCW") {
     warning("With continuous covariates, the used method is 'IPCW'.")
   }

    if(conf == TRUE & conf.type != "bootstrap") {
      warning("This method only allows bootstrap confidence intervals.")
    }


    obj1 <- object
    obj1[[1]] <- cbind(obj1[[1]], xval)
    obj1[[1]] <- na.omit(obj1[[1]])
    colnames(obj1[[1]]) <- c(colnames(object[[1]]), attr(terms(formula),"term.labels"))

    res <- tpIPCW(object = obj1, s = s,
                   z.name = attr(terms(formula),"term.labels"),
                   z.value = z.value, bw = bw, window = window,
      method.weights = method.weights, conf = conf, n.boot = n.boot,
      conf.level = conf.level, cluster = cluster, ncores = ncores)


    class(res) <- c("tpIPCW", "survIDM")
    callp <- paste("pij(s=",s,",t|", attr(terms(formula),"term.labels"), "=", z.value, ")", sep = "")

    # in order to have the same output
    x.nlevels <- 1
    levels <- NULL


  } # end method with numeric or integer covariate







  # factor covariate
  if (length(attr(terms(formula),"term.labels")) > 0 & Class == "factor") {  #LDM/PLMD/KMW by levels of the covariate

    if (!(method %in% c("AJ", "LIDA", "LDM", "PLDM"))){
      stop("A factor is included in the model. Possible methods are 'AJ', 'LIDA', 'LDM' and 'PLDM'." )
    }


    x.nlevels <- nlevels(with(data=data, eval(formula[[3]])))
    levels <- levels(with(data=data, eval(formula[[3]])))

    estim <- list()
    ci <- list()


    for (k in 1:x.nlevels) {
      v.level <- levels(with(data=data, eval(formula[[3]])))[k]
      p<- which(Xval == v.level)
      obj<- object
      obj$data <- object$data[p,]


      #------------------------------

      # AJ method  (no bootstrap)
      if (method == "AJ"){
        res <- tpAJ(object = obj, s = s, conf = conf,
                    conf.level = conf.level, conf.type = conf.type)
        class(res) <- c("AJ", "survIDM")

        add2t <- s
        add2est <- c(s, 1, 0, 0, 1, 0)
        add2CI <- c(1, 1, 0, 0, 0, 0, 1, 1, 0, 0)
        res$t <- c(add2t, res$t)
        res$est <- rbind(add2est, res$est)
        res$CI <- rbind(add2CI, res$CI)

        if (s == 0){
          res$est[, 5:6] <- NA
          res$CI[, 7:10] <- NA
        }

      }


      # LIDA method
      if (method == "LIDA"){

        if(conf == TRUE & conf.type != "bootstrap") {
          warning("This method only allows bootstrap confidence intervals.")
        }

        res <- tpLIDA(object = obj, s = s, conf = conf,
                      conf.level = conf.level, n.boot = n.boot,
                      cluster = cluster, ncores = ncores)
        class(res) <- c("LIDA", "survIDM")
      }


      # LDM method
      if (method == "LDM"){
        res <- tpLDM(object = obj, s = s, conf = conf,
                     conf.level = conf.level, conf.type = conf.type,
                     n.boot = n.boot, cluster = cluster, ncores = ncores)
        class(res) <- c("LDM", "survIDM")
      }


      # PLDM method
      if (method == "PLDM"){

        if(conf == TRUE & conf.type != "bootstrap") {
          warning("This method only allows bootstrap confidence intervals.")
        }


        res <- tpPLDM(object = obj, s = s, conf = conf,
                      conf.level = conf.level, n.boot = n.boot,
                      cluster = cluster, ncores = ncores)
        class(res) <- c("PLDM", "survIDM")
      }

      #------------------------------

      # if (conf == TRUE) {
      #   resu <- data.frame(cbind(res$y, res$estimate, res$LCI, res$UCI))
      #   names(resu) <- c("y", "estimate", paste("lower ",conf.level*100,"% CI", sep=""), paste("upper ",conf.level*100,"% CI", sep=""))
      # }
      #
      # if (conf == FALSE) {
      #   resu <- data.frame(cbind(res$y, res$estimate))
      #   names(resu) <- c("y","estimate")
      # }

      estim[[paste(levels[k])]] <- res$est
      ci[[paste(levels[k])]] <- res$CI

    } # ends the level's loop
    res$est <- estim
    res$CI <- ci

  } # end method with factor covariate






  #callp <- paste("P(T>y|", text3, ")", sep = "")
  callp <- paste("pij(s=",s,",t)", sep = "")
  res$callp <- callp
  res$Nlevels <- x.nlevels
  res$levels <- levels
  res$formula <- formula
  res$call <- match.call()

  return(invisible(res))


}
