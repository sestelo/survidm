#' @title Nonparametric estimation of the Cumulative Incident Functions
#' in the illness-death model
#'
#' @description This function is used to obtain nonparametric estimates of
#' the cumulative incidence probabilities in the illness-death model. They
#' represent the probability of one individual's being or having been in
#' state j at time t.
#'
#' @param formula A \code{formula} object, which must have a \code{survIDM}
#' object as the response on the left of the \code{~} operator and, if desired,
#' a term on the right. The term may be a qualitative or quantitative variable.
#' Without covariates, the right hand side should be \code{~ 1}.
#' @param s The first time for obtaining estimates for the cumulative
#' incidence functions. If missing, 0 will be used.
#' @param data A data.frame including at least four columns named
#' \code{time1}, \code{event1}, \code{Stime} and \code{event}, which correspond
#' to disease free survival time, disease free survival indicator, time to death
#' or censoring, and death indicator, respectively.
#' @param conf Provides pointwise confidence bands. Defaults to \code{FALSE}.
#' @param n.boot The number of bootstrap replicates to compute the variance
#' of the estimator. Default is 199.
#' @param conf.level Level of confidence. Defaults to 0.95 (corresponding to 95\%).
#' @param z.value The value of the covariate on the right hand side of formula
#' at which the cumulative incidence probabilities are computed. For quantitative
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
#'
#'
#'
#' @details Possible options for argument window are \code{"gaussian"},
#' \code{"epanechnikov"}, \code{"tricube"}, \code{"boxcar"},
#' \code{"triangular"}, \code{"quartic"} or \code{"cosine"}.
#'
#' @return An object of class \code{"survIDM"} and one of the following
#' two classes: \code{"CIF"}, and
#' \code{"cifIPCW"}. Objects are implemented as a list with elements:
#'
#' \item{est}{data.frame with estimates of the cumulative incidence probabilities.}
#' \item{CI}{data.frame with the confidence intervals of the cumulative
#' incidence probabilities.}
#' \item{conf.level}{Level of confidence.}
#' \item{s}{The first time for obtaining estimates for the cumulative
#' incidence probabilities.}
#' \item{t}{The time for obtaining the estimates of cumulative incidence probabilities.}
#' \item{conf}{logical; if \code{FALSE} (default) the pointwise confidence
#' bands are not given.}
#' \item{callp}{The expression of the estimated probability.}
#' \item{Nlevels}{The number of levels of the covariate. Provides important
#' information when the covariate at the right hand side of formula
#' is of class factor.}
#' \item{levels}{The levels of the qualitative covariate
#' (if it is of class factor) on the right hand side of formula.}
#' \item{formula}{A formula object.}
#' \item{call}{A call object.}
#'
#' @author Luis Meira-Machado and Marta Sestelo.
#'
#' @references
#' Geskus, R.B. (2011). Cause-specific cumulative incidence estimation
#' and the fine and gray model under both left truncation and right censoring.
#' Biometrics, 67, 39--49.
#'
#' Kalbeisch, J. D. and Prentice R. L. (1980) The statistical analysis of
#' failure time data. John Wiley & Sons, New York.
#'
#' @examples
#'
#' # Cumulative Incidence Function (CIF)
#' res <- CIF(survIDM(time1, event1, Stime, event) ~ 1, data = colonIDM,
#' conf = FALSE)
#' res
#' summary(res, time=365*1:7)
#' plot(res, ylim=c(0, 0.6))
#'
#'
#' # CIF for those in State 1 at time s=365, Y(s)=0
#' res1 <- CIF(survIDM(time1, event1, Stime, event) ~ 1, data = colonIDM,
#' s = 365, conf = FALSE)
#' summary(res1, time=365*1:7)
#' plot(res1, ylim=c(0, 0.6))
#'
#'
#' # Conditional CIF (with a factor)
#' res2 <- CIF(survIDM(time1, event1, Stime, event) ~ factor(sex),
#' data = colonIDM, s = 365, conf = FALSE)
#' summary(res2, time=365*1:5)
#' plot(res2)
#'
#'
#' # Conditional CIF (with continuous covariate)
#' res3 <- CIF(survIDM(time1, event1, Stime, event) ~ age, data = colonIDM,
#' z.value = 56, conf = FALSE)
#' summary(res3, time=365*1:6)
#' plot(res3)
#'




CIF <- function(formula, s, data, conf = FALSE, n.boot = 199,
                conf.level = 0.95, z.value, bw = "dpik", window = "gaussian",
                method.weights = "NW", cluster = FALSE, ncores = NULL){


  if (missing(formula)) stop("A formula argument is required")
  if (missing(s))  s <- 0  #stop("argument 's' is missing, with no default")


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



  # without covariates
  if (length(attr(terms(formula), "term.labels")) == 0) {  #cif without covariate

      res <- CIF_ini(object = object, s = s, conf = conf,
                    conf.level = conf.level, n.boot = n.boot,
                    cluster = cluster, ncores = ncores)
      class(res) <- c("CIF", "survIDM")

  } # end methods without covariate



  # numeric or integer covariate (no argumento s)
  if(length(attr(terms(formula),"term.labels")) != 0 & (Class == "numeric" | Class == "integer")) {#IPCW

    if(s != 0) {warning(paste("If the formula includes a continuous covariate,
                             's' argument will not be used "))}

    obj1 <- object
    obj1[[1]] <- cbind(obj1[[1]], xval)
    obj1[[1]] <- na.omit(obj1[[1]])
    colnames(obj1[[1]]) <- c(colnames(object[[1]]), attr(terms(formula),"term.labels"))

    res <- cifIPCW(object = obj1,
                  z.name = attr(terms(formula),"term.labels"),
                  z.value = z.value, bw = bw, window = window,
                  method.weights = method.weights, conf = conf, n.boot = n.boot,
                  conf.level = conf.level, cluster = cluster, ncores = ncores)


    class(res) <- c("cifIPCW", "survIDM")
    #callp <- paste("CIF(s=",s,",t|", attr(terms(formula),"term.labels"), "=", z.value, ")", sep = "")

  } # end method with numeric or integer covariate


  # factor covariate
  if (length(attr(terms(formula),"term.labels")) > 0 & Class == "factor") {  #CIF by levels of the covariate

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



      res <- CIF_ini(object = obj, s = s, conf = conf,
                 conf.level = conf.level, n.boot = n.boot,
                 cluster = cluster, ncores = ncores)
      class(res) <- c("CIF", "survIDM")



      estim[[paste(levels[k])]] <- res$est
      ci[[paste(levels[k])]] <- res$CI

    } # ends the level's loop
    res$est <- estim
    res$CI <- ci
  }else{
    x.nlevels = 1
    levels = NULL
  } # end method with factor covariate


  #callp <- paste("P(T>y|", text3, ")", sep = "")
  callp <- paste("CIF(s=",s,",t)", sep = "")
  res$callp <- callp
  res$Nlevels <- x.nlevels
  res$levels <- levels
  res$formula <- formula
  res$call <- match.call()

  return(invisible(res))


}
