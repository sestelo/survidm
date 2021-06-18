coxidm <- function (formula, data, semiMarkov = FALSE) {

  object <- data

  if (missing(formula))
    stop("Argument 'formula' is missing with no default")
  if (missing(data))
    stop("Argument 'data' is missing with no default")
  if (class(formula) != "formula")
    stop("Argument 'formula' must be of class 'formula'")
  # if (all(class(colonIDM) != "survIDM"))
  #  stop("Argument 'object' must be of class 'survIDM'")

  ncov <- length(all.vars(formula)) - 4 #number of covariates
  covars <- all.vars(formula)[-c(1:4)] #vector with names of covariates

  term2 <- formula[[3]] #right hand side of formula
  mydata <- object
  #print(names(mydata))

  pos1 <- match(all.vars(formula)[1:4], names(data))
  pos2 <- match(all.vars(formula)[-c(1:4)], names(data))

  newdata <- data[,c(pos1,pos2)]
  names(newdata)[1:4] <- c("time1","event1","Stime","event")

  mydata <- newdata

  #0->2
  p02 <- which(mydata$event == 1 & mydata$Stime == mydata$time1)
  s02 <- rep(0, length(mydata$time1))
  s02[p02] <- 1
  fmla0 <- as.formula(paste("Surv(mydata[,1],s02)~",
                            paste(formula[3], collapse = "+")))

  fit02 <- survival::coxph(fmla0, data = mydata)


  #fit02.anova <- survival::anova(fit02)
  fit02.anova <- anova(fit02)


  #fit02.zph <- survival::cox.zph(fit02)
  fit02.zph <- cox.zph(fit02)


  #0->1
  p01 <- which(mydata$Stime > mydata$time1)
  s01 <- rep(0, length(mydata$time1))
  s01[p01] <- 1
  fmla1 <-
    as.formula(paste("Surv(mydata[,1],s01)~", paste(formula[3], collapse =
                                                      "+")))

  fit01 <- survival::coxph(fmla1, data = mydata)

  #fit01.anova <- survival::anova(fit01)
  fit01.anova <- anova(fit01)


  #fit01.zph <- survival::cox.zph(fit01)
  fit01.zph <- cox.zph(fit01)


  npar <- length(fit01$coef)

  #1->2	CMM vs CsMM
  mydata12 <- mydata[p01, ]
  if (semiMarkov == FALSE) {
    fmla2 <-
      as.formula(paste(
        "Surv(time1, Stime,event)~",
        paste(formula[3], collapse = "+")
      ))
  }
  else {
    fmla2 <-
      as.formula(paste(
        "Surv(Stime - time1,event)~",
        paste(formula[3], collapse = "+")
      ))
  }

  fit12 <- survival::coxph(fmla2, data = mydata12)


  #fit12.anova <- survival::anova(fit12)
  fit12.anova <- anova(fit12)

  #fit12.zph <- survival::cox.zph(fit12)
  fit12.zph <- cox.zph(fit12)


  term01<-termplot(fit01,  se = T, col.term = 1, col.se = 2, plot=FALSE)
  term02<-termplot(fit02,  se = T, col.term = 1, col.se = 2, plot=FALSE)
  term12<-termplot(fit12,  se = T, col.term = 1, col.se = 2, plot=FALSE)


  res <-
    list(
      coxmm01 = fit01,
      coxmm02 = fit02,
      coxmm12 = fit12,

      coxmm01.anova = fit01.anova,
      coxmm02.anova = fit02.anova,
      coxmm12.anova = fit12.anova,

      coxmm01.zph = fit01.zph,
      coxmm02.zph = fit02.zph,
      coxmm12.zph = fit12.zph,

      term01=term01,
      term02=term02,
      term12=term12,

      ncov = ncov,
      fmla = formula,
      npar = npar,
      semiMarkov = semiMarkov
    )
  #class(res) <- c("coxph", "survIDM")
  class(res) <- "cmm"
  return(res)
  #return(invisible(res))
}
