summary.survIDM <- function(object, times = NULL, ...){


  if (inherits(object, "survIDM")) {

    if (class(object)[1] == "data.frame"){ # for the output of survIDM
      summary(object)
    }else{

      if (is.null(times)) { #the whole of the times

        cat("\n")
        cat("Estimation of", object$callp, "\n")
        cat("\n")

        if (object$Nlevels > 1) { # for levels
          res <- list()
          for (i in 1:object$Nlevels) {
            v.level <- object$levels[i]
            cat("   ", attr(terms(object$formula),"term.labels"), "=", v.level, "\n")
            aux <- object$est[[v.level]]
            names(aux) <- c("t", "00", "01", "02", "11", "12")
            print(aux, row.names = FALSE)
            res$est[[v.level]] <- aux
            if(object$conf == TRUE){
              lci <- data.frame(t = object$est[[v.level]][,1], object$CI[[v.level]][,c(1,3,5,7,9)])
              names(lci) <- c("t", "00", "01", "02", "11", "12")
              uci <- data.frame(t = object$est[[v.level]][,1], object$CI[[v.level]][,c(2,4,6,8,10)])
              names(uci) <- c("t", "00", "01", "02", "11", "12")
              cat("\n")
              cat((1-object$conf.level)/2*100,"%", "\n", sep="")
              cat("\n")
              print(lci, row.names = FALSE)
              cat("\n")
              cat((1-(1-object$conf.level)/2)*100,"%", "\n", sep="")
              cat("\n")
              print(uci, row.names = FALSE)
              res$LCI[[v.level]] <- lci
              res$UCI[[v.level]] <- uci
            }
          } #end  levels

        }else{

          if(object$conf == FALSE) {
            res <- list(est = object$est)
            names(res$est) <- c("t", "00", "01", "02", "11", "12")
            print(res$est, row.names = FALSE)
          }else{
            res <- list(est = object$est,
                        LCI = data.frame(t = object$t, object$CI[,c(1,3,5,7,9)]),
                        UCI = data.frame(t = object$t, object$CI[,c(2,4,6,8,10)]))
            names(res$est) <- c("t", "00", "01", "02", "11", "12")
            names(res$LCI) <- c("t", "00", "01", "02", "11", "12")
            names(res$UCI) <- c("t", "00", "01", "02", "11", "12")
            print(res$est, row.names = FALSE)
            cat("\n")
            cat((1-object$conf.level)/2*100,"%", "\n", sep="")
            cat("\n")
            print(res$LCI, row.names = FALSE)
            cat("\n")
            cat((1-(1-object$conf.level)/2)*100,"%", "\n", sep="")
            cat("\n")
            print(res$UCI, row.names = FALSE)

          }

        }
      }else{ #times no null
        if (object$Nlevels > 1) {

          # to control the times argument
          # -----------------------------
          pp <- list()
          for (i in 1:object$Nlevels) {
            pp[[i]] <- sapply(times, function(x)ifelse(x >= min(object$est[[i]][,1]) & x <= max(object$est[[i]][,1]), 1, NA))
          }

          if (all(is.na(unlist(pp)))) {
            stop(paste("At least one element of the 'times' vector has to be in the range of 't' "))
          }

          if (any(is.na(unlist(pp)))) {
            warning(paste("'times' must be in the range of 't' (for each level)" ))
          }
          #  -----
          cat("\n")
          cat("Estimation of", object$callp, "\n")
          cat("\n")


          res <- list()
          for (i in 1:object$Nlevels) {
            v.level <- object$levels[i]
            cat("   ", attr(terms(object$formula),"term.labels"), "=", v.level, "\n")
            ii <- sapply(times, function(x)ifelse( x >= min(object$est[[v.level]][,1]) & x <= max(object$est[[v.level]][,1]),which.max(object$est[[v.level]][,1][object$est[[i]][,1] <= x]), NA))

            if (all(is.na(ii))) {
              aux <- data.frame(times, matrix(NA, nrow = length(times), ncol = dim(object$est[[v.level]])[2] - 1))
            }else{
              aux <- data.frame(times, object$est[[v.level]][ii, -1])
            }
            names(aux) <- c("t", "00", "01", "02", "11", "12")
            print(aux, row.names = FALSE)
            res$est[[v.level]] <- aux
            if(object$conf == TRUE){
              if (all(is.na(ii))) {
                lci <- data.frame(times, matrix(NA, nrow = length(times), ncol = dim(object$est[[v.level]])[2] - 1))
                uci <- data.frame(times, matrix(NA, nrow = length(times), ncol = dim(object$est[[v.level]])[2] - 1))
              }else{
                lci <- data.frame(t = times, object$CI[[v.level]][ii, c(1,3,5,7,9)])
                uci <- data.frame(t = times, object$CI[[v.level]][ii, c(2,4,6,8,10)])
              }
              names(lci) <- c("t", "00", "01", "02", "11", "12")
              names(uci) <- c("t", "00", "01", "02", "11", "12")
              cat("\n")
              cat((1-object$conf.level)/2*100,"%", "\n", sep="")
              cat("\n")
              print(lci, row.names = FALSE)
              cat("\n")
              cat((1-(1-object$conf.level)/2)*100,"%", "\n", sep="")
              cat("\n")
              print(uci, row.names = FALSE)
              res$LCI[[v.level]] <- lci
              res$UCI[[v.level]] <- uci
              cat("\n")
            }
            #--

          }
        }else{ # starts with no levels
          ii <- sapply(times, function(x)ifelse( x >= min(object$t) & x <= max(object$t),
                                                 which.max(object$t[object$t <= x]), NA))
          if (all(is.na(ii))) {
            stop(paste("At least one element of the 'times' vector has to be between",min(object$t), "and", max(object$t)))
          }

          if (any(is.na(ii))) {
            warning(paste("'times' must be between",min(object$t), "and", max(object$t)))
          }

          cat("\n")
          cat("Estimation of", object$callp, "\n")
          cat("\n")

          if(object$conf == FALSE){
            res <- list(est = data.frame(times, object$est[ii, -1]))
            names(res$est) <- c("t", "00", "01", "02", "11", "12")
            print(res$est, row.names = FALSE)
          }else{

            res <- list(est = data.frame(times, object$est[ii, -1]),
                        LCI = data.frame(t = times, object$CI[ii, c(1,3,5,7,9)]),
                        UCI = data.frame(t = times, object$CI[ii, c(2,4,6,8,10)]))
            names(res$est) <- c("t", "00", "01", "02", "11", "12")
            names(res$LCI) <- c("t", "00", "01", "02", "11", "12")
            names(res$UCI) <- c("t", "00", "01", "02", "11", "12")
            print(res$est, row.names = FALSE)
            cat("\n")
            cat((1-object$conf.level)/2*100,"%", "\n", sep="")
            cat("\n")
            print(res$LCI, row.names = FALSE)
            cat("\n")
            cat((1-(1-object$conf.level)/2)*100,"%", "\n", sep="")
            cat("\n")
            print(res$UCI, row.names = FALSE)
            cat("\n")
          }




          #res <- data.frame(times, object$est[ii, -1])
          #names(res) <- names(object$est)
          #print(res, row.names = FALSE)
        }

      }
      class(res) <- c("summary.surv")
      return(invisible(res))
    }
  }else{
    stop("Argument x must be either survIDM object.")
  }
}


