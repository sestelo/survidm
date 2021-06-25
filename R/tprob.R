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
#' Possible options are \code{"AJ"}, \code{"LIDA"} \code{"LM"}, \code{"PLM"},
#' \code{"LMAJ"}, \code{"PLMAJ"}, \code{"PAJ"}, \code{"IPCW"} and \code{"breslow"}.
#' Defaults to \code{"AJ"}. The \code{"IPCW"} method
#' is recommended to obtain conditional transition probabilities (i.e., with a
#' quantitative term on the right hand side of formula). The \code{"breslow"} method
#' is based on a Cox's regression model (Cox, 1972) fitted marginally to each allowed
#' transition, with the corresponding baseline hazard function estimated by the
#' Breslow's method (Breslow, 1972).
#' @param conf Provides pointwise confidence bands. Defaults to \code{FALSE}.
#' @param conf.level Level of confidence. Defaults to 0.95 (corresponding to 95\%).
#' @param conf.type Method to compute the confidence intervals. Depends on the
#' choice of the estimation method of the transition probabilities. For
#' Aalen-Johansen type estimators (\code{"AJ"}, \code{"LMAJ"}, \code{"PAJ"} and
#' \code{"PLMAJ"}) possible choices are \code{"linear"}, \code{"log"} and
#' \code{"log-log"}. Default method is \code{"log"}. The \code{"linear"} option provides
#' the standard intervals curve +-k *se(curve), where k is determined from \code{"conf.int"}.
#' The \code{"log"} option calculates the intervals based on the cumulative hazard or
#' -log(survival). The \code{"log-log"} option uses the log hazard function or
#' log(-log(survival)). For the remaining estimation methods (\code{"LIDA"},
#' \code{"LM"}, \code{"PLM"}, \code{"IPCW"} and \code{"breslow"}) the percentile
#' bootstrap which resamples each datum with probability 1/n is used.
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
#' \code{"triangular"}, \code{"quartic"} or \code{"cosine"}. The LIDA estimator
#' was labelled according to the acronym of the Lifetime Data Analysis journal
#' in which the estimator was described for the first time (Meira-Machado,
#' U?a-?lvarez and Cadarso-Su?rez, 2006).
#'

#' Possible methods are:
#' \itemize{
#' \item{\code{AJ} }{Aalen-Johansen estimator}
#' \item{\code{PAJ} }{Presmoothed Aalen-Johansen estimator}
#' \item{\code{LIDA} }{LIDA estimator}
#' \item{\code{LM} }{Landmark approach estimator}
#' \item{\code{PLM} }{Presmoothed Landmark approach estimator}
#' \item{\code{LMAJ} }{Landmark approach Aalen-Johansen estimator}
#' \item{\code{PLDAJ} }{Presmoothed Landmark approach Aalen-Johansen estimator}
#' \item{\code{tpIPCW} }{Inverse Probability of Censoring Weighting for Transition Probabilities}
#' \item{\code{tpBreslow} }{Breslow method}
#' }
#'
#' @return An object of class \code{"survIDM"} and one of the following
#' five classes: \code{"AJ"}, \code{"LIDA"}, \code{"LM"}, \code{"PLM"},
#' \code{"LMAJ"}, \code{"PLMAJ"}, \code{"PAJ"},
#' \code{"tpIPCW"} and \code{"tpBreslow"}. Objects are implemented as a list with elements:
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
#' @author Luis Meira-Machado, Marta Sestelo and Gustavo Soutinho.
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
#' Cox, DR (1972). Regression models and life tables (with discussion). Journal
#' of the Royal Statistical Society, Series B 34, 187-200.
#'
#' Breslow, N. (1972). Discussion of paper by dr cox. Journal of Royal Statistical
#' Society, Series B 34, 216-217.
#'
#'
#' @examples
#' \dontrun{
#' # Aalen-Johansen
#  Occupation Probabilities Pj(t)=Pij(0,t)

#' res <- tprob(survIDM(time1, event1, Stime, event) ~ 1, s = 0,
#' method = "AJ", conf = FALSE, data = colonIDM)

#' summary(res, time=365*1:6)
#' plot(res)

#' # Transition Probabilities Pij(t)=Pij(365,t)

#' # LIDA
#' res1 <- tprob(survIDM(time1, event1, Stime, event) ~ 1, s = 365,
#'              method = "LIDA", conf = FALSE, data = colonIDM)
#'
#' summary(res1, time=365*1:6)
#' plot(res1)
#' plot(res1, trans="01", ylim=c(0,0.15))

#' # Landmark (LM)
#' res2 <- tprob(survIDM(time1, event1, Stime, event) ~ 1, s = 365,
#'               method = "LM", conf = FALSE, data = colonIDM)

#' summary(res2, time=365*1:6)
#' plot(res2)

#' # Presmoothed LM
#' res3 <- tprob(survIDM(time1, event1, Stime, event) ~ 1, s = 365,
#'               method = "PLM", conf = TRUE, data = colonIDM)

#' summary(res3, time=365*1:6)
#' autoplot(res3, interactive = TRUE)

#' # Conditional transition probabilities

#' # With factor
#' res4 <- tprob(survIDM(time1, event1, Stime, event) ~ factor(sex), s = 365,
#'               method = "AJ", conf = TRUE, data = colonIDM)
#' summary(res4, time=365*1:6)
#' plot(res4, trans="02", ylim=c(0,0.5))

#' res5 <- tprob(survIDM(time1, event1, Stime, event) ~ rx, s =365,
#'               method = "breslow", z.value='Lev', conf = TRUE, data =colonIDM)

#' summary(res5, time=365*1:6)
#' plot(res5,trans="02", ylim=c(0,0.5))


#' # with continuous covariate (IPCW and Breslow Method)
#' res6 <- tprob(survIDM(time1, event1, Stime, event) ~ age, s = 365,
#'               method = "IPCW", z.value = 48, conf = FALSE, data = colonIDM,
#'               bw = "dpik", window = "gaussian", method.weights = "NW")

#' summary(res6, time=365*1:6)
#' plot(res6)

#' res7 <- tprob(survIDM(time1, event1, Stime, event) ~ age, s =365,
#'               method = "breslow", z.value=60, conf = FALSE, data =colonIDM)

#' summary(res7, time=365*1:6)
#' autoplot(res7, interactive=TRUE)

#' res8 <- tprob(survIDM(time1, event1, Stime, event) ~ age, s =365,
#'               method = "breslow", conf.type='bootstrap', z.value=60, conf = TRUE, data =colonIDM)

#' summary(res8, time=365*1:6)
#' plot(res8)

#' res9 <- tprob(survIDM(time1, event1, Stime, event) ~ rx, s =365,
#'               method = "breslow", conf.type='bootstrap',  conf = TRUE, data =colonIDM)

#' summary(res9, time=365*1:6)
#' plot(res9, trans="02", ylim=c(0,0.5))

#' # more than a covariate (Breslow Method)
#' res10<- tprob(survIDM(time1, event1, Stime, event) ~ nodes + factor(rx), s =365,
#'               method = "breslow", conf = TRUE, data =colonIDM)

#' summary(res10,t=365*1:5)
#' autoplot(res10)

#' res11<- tprob(survIDM(time1, event1, Stime, event) ~ nodes + factor(rx), s =365,
#'               method = "breslow", z.value=c(10,'Obs'), conf = TRUE, data =colonIDM)
#' summary(res11,t=365*1:5)
#' autoplot(res11)

#' # more than a covariate for Non Linear Models (Breslow Method)
#' res12<- tprob(survIDM(time1, event1, Stime, event) ~ pspline(age)+ nodes + factor(rx), s =365,
#'               method = "breslow", conf = TRUE, data =colonIDM)

#' summary(res12,t=365*1:5)
#' autoplot(res12)

#' # Confidence intervals
#' res13 <- tprob(survIDM(time1, event1, Stime, event) ~ 1, s = 365,
#'                method = "AJ", conf = TRUE, n.boot = 5, conf.level = 0.95,
#'                conf.type = "log", data = colonIDM)

#' summary(res13, time=365*1:7)
#' autoplot(res13)
#' }


tprob<-function(formula,
                s,
                method = "AJ",
                conf = FALSE,
                conf.level = 0.95,
                conf.type = "log",
                n.boot = 199,
                data,
                z.value,
                bw = "dpik",
                window = "gaussian",
                method.weights = "NW",
                cluster = FALSE,
                ncores = NULL,
                na.rm = TRUE) {


  if (missing(formula))
    stop("A formula argument is required")
  if (missing(s))
    #stop("argument 's' is missing, with no default")
    s<-0

  if (!(method %in% c("AJ", "LIDA", "LM", "PLM", "IPCW", "LMAJ", "PLMAJ", "PAJ", "breslow"))) {
    stop(
      "Possible methods are 'AJ', 'LIDA', 'LM', 'PLM', 'LMAJ', 'PAJ', 'PLMAJ', 'breslow' and 'IPCW'."
    )
  }

  fmla <- eval(formula, parent.frame())
  Terms <- terms(fmla)

  mf <- Terms[[2]]
  object <- with(data = data, eval(mf))

  object <-
    list(data = object)
  obj_data <- object[[1]]


  if (length(attr(terms(formula), "term.labels")) <= 1 & !('pspline' %in% substr(attr(terms(formula), "term.labels"),
                                                                                 1,7))){

    X <- Terms[[3]]

    Class <- class(with(data = data, eval(Terms[[3]])))


    if (Class != "numeric" & Class != "integer" & Class != "factor")

      stop("covariate must be one of the following classes 'numeric', 'integer' or 'factor'")

    xval <- with(data = data, eval(X))

    lencov <- length(xval)

    lencov2 <- length(attr(terms(formula), "term.labels"))

    if (lencov2 != 0)
      Xval <- with(data = data, eval(Terms[[3]]))

    if (lencov != dim(obj_data)[1] &
        lencov2 != 0)
      stop("length of the covariate does not match")

    #--------------------------------------------------------------------------------------------------------
    # without covariates
    if (length(attr(terms(formula), "term.labels")) == 0) {
      #AJ, LIDA, LM, PLM without covariate

      if (!(method %in% c("AJ", "LIDA", "LM", "PLM", "LMAJ", "PAJ", "PLMAJ"))) {
        stop(
          "The model does not include covariates. Possible methods are 'AJ', 'LIDA', 'LM', 'PLM',
            'LMAJ', 'PAJ' and 'PLMAJ'."
        )
      }

      # AJ method  (no bootstrap)
      if (method == "AJ") {
        res <- tpAJ(
          object = object,
          s = s,
          conf = conf,
          conf.level = conf.level,
          conf.type = conf.type
        )
        class(res) <- c("AJ", "survIDM")

        add2t <- s
        add2est <- c(s, 1, 0, 0, 1, 0)
        add2CI <- c(1, 1, 0, 0, 0, 0, 1, 1, 0, 0)
        res$t <- c(add2t, res$t)
        res$est <- rbind(add2est, res$est)
        res$CI <- rbind(add2CI, res$CI)

        if (s == 0) {
          res$est[, 5:6] <- NA
          res$CI[, 7:10] <- NA
        }

      }


      # LIDA method
      if (method == "LIDA") {
        if (conf == TRUE & conf.type != "bootstrap") {
          warning("This method only allows bootstrap confidence intervals.")
        }

        res <- tpLIDA(
          object = object,
          s = s,
          conf = conf,
          conf.level = conf.level,
          n.boot = n.boot,
          cluster = cluster,
          ncores = ncores
        )
        class(res) <- c("LIDA", "survIDM")
      }


      # LM method
      if (method == "LM") {
        res <- tpLM(
          object = object,
          s = s,
          conf = conf,
          conf.level = conf.level,
          conf.type = conf.type,
          n.boot = n.boot,
          cluster = cluster,
          ncores = ncores
        )
        class(res) <- c("LM", "survIDM")
      }


      # PLM method
      if (method == "PLM") {
        if (conf == TRUE & conf.type != "bootstrap") {
          warning("This method only allows bootstrap confidence intervals.")
        }

        res <- tpPLM(
          object = object,
          s = s,
          conf = conf,
          conf.level = conf.level,
          n.boot = n.boot,
          cluster = cluster,
          ncores = ncores
        )
        class(res) <- c("PLM", "survIDM")
      }


      # LMAJ method
      if (method == "LMAJ") {
        res <- tpLMAJ(
          object = object,
          s = s,
          conf = conf,
          conf.level = conf.level
        )
        class(res) <- c("LMAJ", "survIDM")
      }



      # PAJ method
      if (method == "PAJ") {
        res <- tpPAJ(
          object = object,
          s = s,
          conf = conf,
          conf.level = conf.level
        )
        class(res) <- c("PAJ", "survIDM")
      }


      # PLMAJ method
      if (method == "PLMAJ") {
        res <- tpPLMAJ(
          object = object,
          s = s,
          conf = conf,
          conf.level = conf.level
        )
        class(res) <- c("PLMAJ", "survIDM")
      }


      # in order to have the same output
      x.nlevels <- 1
      levels <- NULL

    } # end methods without covariate

    #--------------------------------------------------------------------------------------------------------

    # numeric or integer covariate:

    if (length(attr(terms(formula), "term.labels")) != 0 &
        (Class == "numeric" | Class == "integer")) {


      if (method != c("IPCW") & method != c("breslow")) {
        warning("With continuous covariates, the used method is 'IPCW' or 'breslow.")
      }


      obj1 <- object

      obj1[[1]] <- cbind(obj1[[1]], xval)

      obj1[[1]] <- na.omit(obj1[[1]])

      colnames(obj1[[1]]) <-

        c(colnames(object[[1]]), attr(terms(formula), "term.labels"))


      if (method == c("IPCW")) {


        if (conf == TRUE & conf.type != "bootstrap") {
          warning("This method only allows bootstrap confidence intervals.")
        }


        res <- tpIPCW(
          object = obj1,
          s = s,
          z.name = attr(terms(formula), "term.labels"),
          z.value = z.value,
          bw = bw,
          window = window,
          method.weights = method.weights,
          conf = conf,
          n.boot = n.boot,
          conf.level = conf.level,
          cluster = cluster,
          ncores = ncores
        )

        class(res)

        class(res) <- c("tpIPCW", "survIDM")

        callp <-
          paste("pij(s=",
                s,
                ",t|",
                attr(terms(formula), "term.labels"),
                "=",
                z.value,
                ")",
                sep = "")

        #> callp
        #[1] "pij(s=365,t|age=48)"

        x.nlevels <- 1

        levels <- NULL

      } #fim method IPCW


      if (method == c("breslow")) {


        #if (!exists("z.value"))

        if (missing(z.value))

          z.value<-NULL

        if (is.null(z.value))

          z.value<-mean(xval)


        if (missing(conf.type)){
          conf.type<-'log'

        }

        if(conf.type == "linear"){

          conf.type<-'plain'
        }


        res <- tpBreslow(

          object = obj1,
          s = s,
          z.name = attr(terms(formula), "term.labels"),
          z.value = z.value,
          conf = conf,
          conf.type=conf.type,
          conf.level=conf.level,
          parte2=''

        )


        class(res) <- c("tpBreslow", "survIDM")

        callp <-
          paste("pij(s=",
                s,
                ",t|",
                attr(terms(formula), "term.labels"),
                "=",
                z.value,
                ")",
                sep = "")

        # in order to have the same output
        x.nlevels <- 1
        levels <- NULL

      }


    }


    #-----------------------------------------------------------------------------------------

    # factor covariate
    if (length(attr(terms(formula), "term.labels")) > 0 & Class == "factor") {
      #LM/PLMD/KMW by levels of the covariate


      if (!(method %in% c("AJ", "LIDA", "LM", "PLM", "LMAJ", "PAJ", "PLMAJ","breslow"))) {

        stop(
          "A factor is included in the model. Possible methods are 'AJ', 'LIDA', 'LM', 'PLM',
            'LMAJ', 'PAJ', 'PLMAJ' and 'breslow'."
        )
      }


      if (method %in% c("AJ", "LIDA", "LM", "PLM", "LMAJ", "PAJ", "PLMAJ")) {

        x.nlevels <- nlevels(with(data = data, eval(formula[[3]])))
        levels <- levels(with(data = data, eval(formula[[3]])))

        estim <- list()
        ci <- list()

        for (k in 1:x.nlevels) {

          v.level <- levels(with(data = data, eval(formula[[3]])))[k]
          p <- which(Xval == v.level)
          obj <- object
          obj$data <- object$data[p, ]

          #------------------------------

          # AJ method  (no bootstrap)

          if (method == "AJ") {
            res <- tpAJ(
              object = obj,
              s = s,
              conf = conf,
              conf.level = conf.level,
              conf.type = conf.type
            )
            class(res) <- c("AJ", "survIDM")

            add2t <- s
            add2est <- c(s, 1, 0, 0, 1, 0)
            add2CI <- c(1, 1, 0, 0, 0, 0, 1, 1, 0, 0)
            res$t <- c(add2t, res$t)
            res$est <- rbind(add2est, res$est)
            res$CI <- rbind(add2CI, res$CI)

            if (s == 0) {
              res$est[, 5:6] <- NA
              res$CI[, 7:10] <- NA
            }

          }


          # LIDA method
          if (method == "LIDA") {
            if (conf == TRUE & conf.type != "bootstrap") {
              warning("This method only allows bootstrap confidence intervals.")
            }

            res <- tpLIDA(
              object = obj,
              s = s,
              conf = conf,
              conf.level = conf.level,
              n.boot = n.boot,
              cluster = cluster,
              ncores = ncores
            )
            class(res) <- c("LIDA", "survIDM")
          }


          # LM method
          if (method == "LM") {
            res <- tpLM(
              object = obj,
              s = s,
              conf = conf,
              conf.level = conf.level,
              conf.type = conf.type,
              n.boot = n.boot,
              cluster = cluster,
              ncores = ncores
            )
            class(res) <- c("LM", "survIDM")
          }


          # PLM method
          if (method == "PLM") {
            if (conf == TRUE & conf.type != "bootstrap") {
              warning("This method only allows bootstrap confidence intervals.")
            }


            res <- tpPLM(
              object = obj,
              s = s,
              conf = conf,
              conf.level = conf.level,
              n.boot = n.boot,
              cluster = cluster,
              ncores = ncores
            )
            class(res) <- c("PLM", "survIDM")
          }

          # LMAJ method

          if (method == "LMAJ") {
            res <- tpLMAJ(
              object = obj,
              s = s,
              conf = conf,
              conf.level = conf.level
            )
            class(res) <- c("LMAJ", "survIDM")

          }


          # PAJ method
          if (method == "PAJ") {
            res <- tpPAJ(
              object = obj,
              s = s,
              conf = conf,
              conf.level = conf.level
            )
            class(res) <- c("PAJ", "survIDM")
          }


          # PLMAJ method
          if (method == "PLMAJ") {
            res <- tpPLMAJ(
              object = obj,
              s = s,
              conf = conf,
              conf.level = conf.level
            )
            class(res) <- c("PLMAJ", "survIDM")
          }


          estim[[paste(levels[k])]] <- res$est
          ci[[paste(levels[k])]] <- res$CI

        } # ends the level's loop

        res$est <- estim

        res$CI<-ci


      }

      if (method == c("breslow")) {

        obj1 <- object

        obj1[[1]] <- cbind(obj1[[1]], xval)

        obj1[[1]] <- na.omit(obj1[[1]])

        colnames(obj1[[1]]) <-

          c(colnames(object[[1]]), attr(terms(formula), "term.labels"))

        estim <- list()

        ci <- list()

        if (missing(conf.type)){
          conf.type<-'log'
        }



        #if (!exists("z.value"))

        if (missing(z.value))

          z.value<-NULL

        if (is.null(z.value))

          z.value<-levels(with(data = data, eval(formula[[3]])))

        if(conf.type == "linear"){

          conf.type<-'plain'
        }

        for (k in 1:length(z.value)) {

          #k<-1

          res <- tpBreslow(

            object = obj1,
            s = s,
            z.name = attr(terms(formula), "term.labels"),
            z.value = z.value[k],
            conf = conf,
            conf.type=conf.type,
            conf.level=conf.level,
            parte2=''

          )


          class(res) <- c("tpBreslow", "survIDM")

          estim[[paste(z.value[k])]] <- res$est

          ci[[paste(z.value[k])]] <- res$CI

        } # ends the level's loop

        res$est <- estim

        res$CI <- ci



        if(length(z.value)==1){
          res$est<-as.data.frame(res$est)
          res$CI<-as.data.frame(res$CI)

        }

        callp <-
          paste("pij(s=",
                s,
                ",t|",
                attr(terms(formula), "term.labels"),
                "=",
                z.value,
                ")",
                sep = "")

        # in order to have the same output
        #x.nlevels <- length(levels(xval))
        x.nlevels <- length(z.value)
        levels <- levels(xval)
        #x.nlevels <- 1
        #levels <- NULL
      }


    } # end method with factor covariate


  }else{

    cov<-attr(terms(formula), "term.labels")

    covF<-rep(NA, length(cov))

    obj1 <- object

    for(i in 1:length(cov)){

      if(substr(cov[i],1,6)=='factor'){
        covF[i]<-substr(cov[i],8,nchar(cov[i])-1)
      }else

        if(substr(cov[i],1,7)=='pspline'){
          covF[i]<-substr(cov[i],9,nchar(cov[i])-1)
        }else{

          covF[i]<-cov[i]
        }

      xval <- with(data = data,data[,which(colnames(data)==covF[i])])

      obj1[[1]] <- cbind(obj1[[1]], xval)

    }

    obj1[[1]] <- na.omit(obj1[[1]])

    colnames(obj1[[1]])[5:(5+length(covF)-1)] <-covF

    #if (!exists("z.value"))

    if (missing(z.value))

      z.value<-NULL


    if (method == c("breslow")) {

      if (missing(conf.type)){

        conf.type<-'log'

      }

      if(conf.type == "linear"){

        conf.type<-'plain'
      }

      res <- tpBreslow(

        object = obj1,
        s = s,
        z.name =colnames(obj1[[1]])[5:(5+length(covF)-1)],
        z.value = z.value,
        conf = conf,
        conf.type=conf.type,
        conf.level=conf.level,
        parte2=attr(terms(formula), "term.labels")

      )

      #dim(res$est)
      #table(res$est$t %in% times2)
      #times1 %in%  res$est$t

      #res$est[1:20,]
      #res$CI[1:20,1:4]
      #res$est[(length(times2)-10):length(times2),]
      #res$CI[(length(times2)-10):length(times2),]

      class(res) <- c("tpBreslow", "survIDM")

      #summary(res,times=c(365*2:6,2826))

      callp <-
        paste("pij(s=",
              s,
              ",t|",
              attr(terms(formula), "term.labels"),
              "=",
              z.value,
              ")",
              sep = "")

      # in order to have the same output
      x.nlevels <- 1
      levels <- NULL

    }

  }#end more than a covariate

  #callp <- paste("P(T>y|", text3, ")", sep = "")
  callp <- paste("pij(s=", s, ",t)", sep = "")
  res$callp <- callp
  res$Nlevels <- x.nlevels
  res$levels <- levels
  res$formula <- formula
  res$call <- match.call()
  options(warn=-1)
  suppressWarnings(res)
  return(invisible(res))

}


