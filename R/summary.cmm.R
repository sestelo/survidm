summary.cmm <- function(object, type=NULL, conf.level = 0.95, ...){

  if(is.null(type)){

    #model <- object$coxmm01
    tab01 <- summary(object$coxmm01, conf.int = conf.level)$conf.int
    tab01 <- cbind(tab01, summary(object$coxmm01)$coef[, 5])
    tab01 <- cbind(summary(object$coxmm01)$coef[, 1], tab01)
    aa <- rownames(tab01)
    tab01 <- data.frame(tab01)
    tab01 <- tab01[, -3]
    tab01 <- as.matrix(tab01)
    rownames(tab01) <- aa
    colnames(tab01) <-
      c(
        "coef",
        "exp(coef)",
        paste("lower", conf.level),
        paste("upper", conf.level),
        "Pr(>|z|)"
      )

    #model <- object$coxmm02
    tab02 <- summary(object$coxmm02, conf.int = conf.level)$conf.int
    tab02 <- cbind(tab02, summary(object$coxmm02)$coef[, 5])
    tab02 <- cbind(summary(object$coxmm02)$coef[, 1], tab02)
    aa <- rownames(tab02)
    tab02 <- data.frame(tab02)
    tab02 <- tab02[, -3]
    tab02 <- as.matrix(tab02)
    rownames(tab02) <- aa
    colnames(tab02) <-
      c(
        "coef",
        "exp(coef)",
        paste("lower", conf.level),
        paste("upper", conf.level),
        "Pr(>|z|)"
      )

    tab12 <- summary(object$coxmm12, conf.int = conf.level)$conf.int
    tab12 <- cbind(tab12, summary(object$coxmm12)$coef[, 5])
    tab12 <- cbind(summary(object$coxmm12)$coef[, 1], tab12)
    aa <- rownames(tab12)
    tab12 <- data.frame(tab12)
    tab12 <- tab12[, -3]
    tab12 <- as.matrix(tab12)
    rownames(tab12) <- aa
    colnames(tab12) <-
      c(
        "coef",
        "exp(coef)",
        paste("lower", conf.level),
        paste("upper", conf.level),
        "Pr(>|z|)"
      )

    cat("Cox Markov Model: transition 0 -> 1", "\n")
    cat("\n")
    print(tab01)

    cat("\n")
    cat("\n")
    cat("Cox Markov Model: transition 0 -> 2", "\n")
    cat("\n")
    print(tab02)

    cat("\n")
    cat("\n")
    if (object$semiMarkov == FALSE)
      cat("Cox Markov Model: transition 1 -> 2", "\n")
    else
      cat("Cox semi-Markov Model: transition 1 -> 2", "\n")
    cat("\n")
    print(tab12)

    res <-
      list(
        cmm.idm.01 = object$coxmm01,
        cmm.idm.02 = object$coxmm02,
        cmm.idm.12 = object$coxmm12
      )
    return(invisible(res))

  }else{

    if(type=='anova'){

      cat("Cox Markov Model: transition 0 -> 1", "\n")
      cat("\n")
      print(object$coxmm01.anova)

      cat("\n")
      cat("\n")
      cat("Cox Markov Model: transition 0 -> 2", "\n")
      cat("\n")
      print(object$coxmm02.anova)

      cat("\n")
      cat("\n")
      if (object$semiMarkov == FALSE)
        cat("Cox Markov Model: transition 1 -> 2", "\n")
      else
        cat("Cox semi-Markov Model: transition 1 -> 2", "\n")
      cat("\n")
      print(object$coxmm12.anova)

    }

    if(type=='ph'){


      cat("Cox Markov Model: transition 0 -> 1", "\n")
      cat("Test the Proportional Hazards Assumption", "\n")
      cat("\n")
      print(object$coxmm01.zph)

      cat("\n")
      cat("\n")
      cat("Cox Markov Model: transition 0 -> 2", "\n")
      cat("Test the Proportional Hazards Assumption", "\n")
      cat("\n")
      print(object$coxmm02.zph)

      cat("\n")
      cat("\n")
      if (object$semiMarkov == FALSE){
        cat("Cox Markov Model: transition 1 -> 2", "\n")
        cat("Test the Proportional Hazards Assumption", "\n")
      }else{

        cat("Cox semi-Markov Model: transition 1 -> 2", "\n")
        cat("Test the Proportional Hazards Assumption", "\n")
      }

      cat("\n")
      print(object$coxmm12.zph)
    }

    if(type!='ph' & type!='anova'){

      cat("Possible option for type are 'anova' or 'ph'")

    }

  }

}
