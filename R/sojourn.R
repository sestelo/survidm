#' @title Nonparametric estimation of the Sojourn time distributions in the
#' recurrence state in the illness-death model
#'
#' @description This function is used to obtain nonparametric estimates of
#' of the sojourn probabilities in the recurrence state in the illness-death
#' model.
#'
#' @param formula A \code{formula} object, which must have a \code{survIDM}
#' object as the response on the left of the \code{~} operator and, if desired,
#' a term on the right. The term may be a qualitative or quantitative variable.
#' Without covariates, the right hand side should be \code{~ 1}.
#' @param data A data.frame including at least four columns named
#' \code{time1}, \code{event1}, \code{Stime} and \code{event}, which correspond
#' to disease free survival time, disease free survival indicator, time to death
#' or censoring, and death indicator, respectively.
#' @param conf Provides pointwise confidence bands. Defaults to \code{FALSE}.
#' @param n.boot The number of bootstrap replicates to compute the variance
#' of the non-Markovian estimator. Default is 199.
#' @param conf.level Level of confidence. Defaults to 0.95 (corresponding to 95\%).
#' @param z.value The value of the covariate on the right hand side of formula
#' at which the sojourn probabilities are computed. For quantitative
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
#' @param method The method used to compute the sojourn estimates.
#' Possible options are \code{"LM"} and \code{"Satten-Datta"}.
#' Defaults to \code{"LM"}.
#' @param presmooth - A logical value. If \code{TRUE}, the presmoothed landmark
#' estimator of the sojourn function is computed. Only valid for \code{method = "LM"}.
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
#'#'
#' @details Possible options for argument window are \code{"gaussian"},
#' \code{"epanechnikov"}, \code{"tricube"}, \code{"boxcar"},
#' \code{"triangular"}, \code{"quartic"} or \code{"cosine"}.
#'
#' @return An object of class \code{"survIDM"} and one of the following
#' two classes: \code{"soj"} (Sojourn Time Distribution), and
#' \code{"sojIPCW"} (Inverse Probability of Censoring Weighting for the Sojourn Time Distribution). Objects are implemented as a list with elements:
#'
#' \item{est}{data.frame with estimates of the sojourn probabilities.}
#' \item{CI}{data.frame with the confidence intervals of the sojourn
#' probabilities.}
#' \item{conf.level}{Level of confidence.}
#' \item{t}{The time for obtaining the estimates of sojourn probabilities.}
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
#' @author Luis Meira-Machado, Marta Sestelo and Gustavo Soutinho.
#'
#' @references
#' Satten, G.A. and Datta, S. (2002) Marginal estimation for
#' multi-stage models: waiting time distributions and competing risks
#' analyses. Statistics in Medicine, 21, 3--19.
#'
#'
#' @examples
#'
#' res <- sojourn(survIDM(time1, event1, Stime, event) ~ 1,
#' data = colonIDM, conf = FALSE, conf.level = 0.95)
#' res
#' summary(res, time=365*1:6)
#' plot(res)
#'
#' res1 <- sojourn(survIDM(time1, event1, Stime, event) ~ 1,
#' data = colonIDM, conf = FALSE, conf.level = 0.95, method = "LM",
#' presmooth = TRUE)
#' res1
#' summary(res1, time=365*1:6)
#' plot(res1)
#'
#'
#' # not run:
#' #res2 <- sojourn(survIDM(time1, event1, Stime, event) ~ 1,
#' #data = colonIDM, conf = FALSE, conf.level = 0.95, method = "Satten-Datta")
#' #res2
#'
#'
#' # with a factor
#' res3 <- sojourn(survIDM(time1, event1, Stime, event) ~ factor(sex),
#' data = colonIDM, conf = FALSE, conf.level = 0.95)
#' res3
#' summary(res3, time=365*1:6)
#' plot(res3)
#'
#' # with a qualitative covariate
#' res4 <- sojourn(survIDM(time1, event1, Stime, event) ~ age, data = colonIDM,
#' z.value = 56, conf = FALSE, conf.level = 0.95)
#' res4
#' summary(res4, time=365*1:6)
#' plot(res4)
#'

sojourn <- function(formula, data, conf = FALSE, n.boot = 199,
                    conf.level = 0.95, z.value, bw = "dpik", window = "gaussian",
                    method.weights = "NW", method = "LM", presmooth = FALSE,
                    cluster = FALSE, ncores = NULL){


  if (missing(formula)) stop("A formula argument is required")
  # if (missing(s)) stop("argument 's' is missing, with no default")


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
  if (length(attr(terms(formula), "term.labels")) == 0) {  #sojourn without covariate

    res <- sojourn_ini(object = object, conf = conf,
                       conf.level = conf.level, n.boot = n.boot,
                       cluster = cluster, ncores = ncores, method = method,
                       presmooth = presmooth)
    res$levels <- NULL
    class(res) <- c("soj", "survIDM")

  } # end methods without covariate


  # numeric or integer covariate
  if(length(attr(terms(formula),"term.labels")) != 0 & (Class == "numeric" | Class == "integer")) {#IPCW

    obj1 <- object
    obj1[[1]] <- cbind(obj1[[1]], xval)
    obj1[[1]] <- na.omit(obj1[[1]])
    colnames(obj1[[1]]) <- c(colnames(object[[1]]), attr(terms(formula),"term.labels"))

    res <- sojIPCW(object = obj1,
                   z.name = attr(terms(formula),"term.labels"),
                   z.value = z.value, bw = bw, window = window,
                   method.weights = method.weights, conf = conf, n.boot = n.boot,
                   conf.level = conf.level, cluster = cluster, ncores = ncores)


    class(res) <- c("sojIPCW", "survIDM")
    #callp <- paste("CIF(s=",s,",t|", attr(terms(formula),"term.labels"), "=", z.value, ")", sep = "")

  } # end method with numeric or integer covariate


  # factor covariate
  if (length(attr(terms(formula),"term.labels")) > 0 & Class == "factor") {  #sojourn by levels of the covariate

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



      res <- sojourn_ini(object = obj, conf = conf,
                         conf.level = conf.level, n.boot = n.boot,
                         cluster = cluster, ncores = ncores, method = method,
                         presmooth = presmooth)
      class(res) <- c("soj", "survIDM")



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
  callp <- paste("sojourn(t)", sep = "")
  res$callp <- callp
  res$Nlevels <- x.nlevels
  res$levels <- levels
  res$formula <- formula
  res$call <- match.call()

  return(invisible(res))

}
