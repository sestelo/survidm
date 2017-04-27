cif <- function(formula, t, s, data, conf = FALSE, n.boot = 199,
                conf.level = 0.95, z.value, bw = "dpik", window = "gaussian",
                method.weights = "NW", cluster = FALSE, ncores = NULL){


  if (missing(formula)) stop("A formula argument is required")
  if (missing(s)) stop("argument 's' is missing, with no default")


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

      res <- CIF_ini(object = object, t = t, s = s, conf = conf,
                    conf.level = conf.level, n.boot = n.boot,
                    cluster = cluster, ncores = ncores)
      class(res) <- c("CIF", "survIDM")

  } # end methods without covariate



  # numeric or integer covariate
  if(length(attr(terms(formula),"term.labels")) != 0 & (Class == "numeric" | Class == "integer")) {#IPCW

    obj1 <- object
    obj1[[1]] <- cbind(obj1[[1]], xval)
    obj1[[1]] <- na.omit(obj1[[1]])
    colnames(obj1[[1]]) <- c(colnames(object[[1]]), attr(terms(formula),"term.labels"))

    res <- cifIPCW(object = obj1, t = t,
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



      res <- CIF_ini(object = obj, t = t, s = s, conf = conf,
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
