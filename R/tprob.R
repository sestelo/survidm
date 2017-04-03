tprob <- function(formula, s, t, method = "AJ", conf = TRUE, conf.level = 0.95,
                  conf.type = "log", n.boot = 200, data, z.value, bw = "dpik",
                  window = "gaussian", method.weights = "NW", cluster = FALSE,
                  ncores = NULL, na.rm = TRUE){


  if (missing(formula)) stop("A formula argument is required")
  if (missing(s)) stop("argument 's' is missing, with no default")

  if (!(method %in% c("AJ", "LIDA", "LDM", "PLDM"))){
    stop("Possible methods are 'AJ', 'LIDA', 'LDM' and 'PLDM'." )
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
  if(length(attr(terms(formula),"term.labels")) > 1) stop("only one covariate is supported")
  Class <- class(with(data = data, eval(Terms[[3]])))
  if (Class != "numeric" & Class != "integer" & Class != "factor") stop("covariate must be one of the following classes 'numeric', 'integer' or 'factor'")
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

    # AJ method  (no bootstrap)
    if (method == "AJ"){
      res <- tpAJ(object = object, s = s, t = t, conf = conf,
                  conf.level = conf.level, conf.type = conf.type)
      class(res) <- c("AJ", "survIDM")
    }


    # LIDA method
    if (method == "LIDA"){
      res <- tpLIDA(object = object, s = s, t = t, conf = conf,
                    conf.level = conf.level, n.boot = n.boot,
                    cluster = cluster, ncores = ncores)
      class(res) <- c("LIDA", "survIDM")
    }


    # LDM method
    if (method == "LDM"){
      res <- tpLDM(object = object, s = s, t = t, conf = conf,
                   conf.level = conf.level, conf.type = conf.type,
                   n.boot = n.boot, cluster = cluster, ncores = ncores)
      class(res) <- c("LDM", "survIDM")
    }


    # PLDM method
    if (method == "PLDM"){
      res <- tpPLDM(object = object, s = s, t = t, conf = conf,
                    conf.level = conf.level, n.boot = n.boot,
                    cluster = cluster, ncores = ncores)
      class(res) <- c("PLDM", "survIDM")
    }
  } # end methods without covariate








  # numeric or integer covariate
  if(length(attr(terms(formula),"term.labels")) != 0 & (Class == "numeric" | Class == "integer")) {#IPCW



    # falta esto

  } # end method with numeric or integer covariate




  #
  #
  #
  #
  # # factor covariate
  # if (length(attr(terms(formula),"term.labels")) > 0 & Class == "factor") {  #LDM/PLMD/KMW by levels of the covariate
  #
  #   x.nlevels <- nlevels(with(data=data, eval(formula[[3]])))
  #   levels <- levels(with(data=data, eval(formula[[3]])))
  #
  #   estim <- c()
  #   LCI <- c()
  #   UCI <- c()
  #
  #   for (k in 1:x.nlevels) {
  #     v.level <- levels(with(data=data, eval(formula[[3]])))[k]
  #     p<- which(Xval == v.level)
  #     obj<- object
  #     obj$data <- object$data[p,]
  #
  #
  #     # AJ method  (no bootstrap)
  #     if (method == "AJ"){
  #       res <- tpAJ(object = obj, s = s, t = t, conf = conf,
  #                   conf.level = conf.level, conf.type = conf.type)
  #       class(res) <- c("AJ", "survIDM")
  #     }
  #
  #
  #     # LIDA method
  #     if (method == "LIDA"){
  #       res <- tpLIDA(object = obj, s = s, t = t, conf = conf,
  #                     conf.level = conf.level, conf.type = conf.type,
  #                     n.boot = n.boot, cluster = cluster, ncores = ncores)
  #       class(res) <- c("LIDA", "survIDM")
  #     }
  #
  #
  #     # LDM method
  #     if (method == "LDM"){
  #       res <- tpLDM(object = obj, s = s, t = t, conf = conf,
  #                    conf.level = conf.level, conf.type = conf.type,
  #                    n.boot = n.boot, cluster = cluster, ncores = ncores)
  #       class(res) <- c("LDM", "survIDM")
  #     }
  #
  #
  #     # PLDM method
  #     if (method == "PLDM"){
  #       res <- tpPLDM(object = obj, s = s, t = t, conf = conf,
  #                     conf.level = conf.level, conf.type = conf.type,
  #                     n.boot = n.boot, cluster = cluster, ncores = ncores)
  #       class(res) <- c("PLDM", "survIDM")
  #     }
  #
  #     if (conf == TRUE) {
  #       resu <- data.frame(cbind(res$y, res$estimate, res$LCI, res$UCI))
  #       names(resu) <- c("y", "estimate", paste("lower ",conf.level*100,"% CI", sep=""), paste("upper ",conf.level*100,"% CI", sep=""))
  #     }
  #
  #     if (conf == FALSE) {
  #       resu <- data.frame(cbind(res$y, res$estimate))
  #       names(resu) <- c("y","estimate")
  #     }
  #
  #     estim[[paste(levels[k])]] <- resu
  #
  #
  #
  #
  #
  #
  #
  # } # end method with factor covariate
  #
  #
  #
  return(res)

}
