summary.survIDM <- function(object, times = NULL, ...){


  if (inherits(object, "survIDM")) {

    if (class(object)[1] == "data.frame"){
      summary(object)
    }else{

      if (is.null(times)) { #the whole of the times

        cat("\n")
        cat(object$callp, "\n")
        cat("\n")

        if (object$Nlevels > 1) { # for levels
          for (i in 1:object$Nlevels) {
            v.level <- object$levels[i]
            cat("   ", attr(terms(object$formula),"term.labels"), "=", v.level, "\n")
            print(object$est[[i]], row.names = FALSE)
            cat("\n")
          }
        }else{
          print(object$est, row.names = FALSE)
        }
        res <- ifelse(conf == FALSE,
                      list(est = object$est),
                      list(est = object$est,
                           LCI = data.frame(t = object$t, object$CI[,c(1,3,5,7,9)]),
                           UCI = data.frame(t = object$t, object$CI[,c(2,4,6,8,10)])))



      }else{ #times no null
        if (object$Nlevels > 1) {

          # to control the times argument
          # -----------------------------
          pp <- list()
          for (i in 1:object$Nlevels) {
            pp[[i]] <- sapply(times, function(x)ifelse(x >= min(object$est[[i]]$y) & x <= max(object$est[[i]]$y), 1, NA))
          }
          if (all(is.na(unlist(pp)))) {
            stop(paste("At least one element of the 'times' vector has to be in the range of 'y' "))
          }

          if (any(is.na(unlist(pp)))) {
            warning(paste("'times' must be in the range of 'y' (for each level)" ))
          }
          #--------------------

          res <- list()
          for (i in 1:object$Nlevels) {
            v.level <- object$levels[i]
            ii <- sapply(times, function(x)ifelse( x >= min(object$est[[i]]$y) & x <= max(object$est[[i]]$y),
                                                   which.max(object$est[[i]]$y[object$est[[i]]$y <= x]), NA))

            if (all(is.na(ii))) {
              aux <- data.frame(times, matrix(NA, nrow = length(times), ncol = dim(object$est[[i]])[2] - 1))
            }else{
              aux <- data.frame(times, object$est[[i]][ii,-1 ])
            }
            names(aux) <- names(object$est[[i]])
            res[[paste(object$levels[i])]] <- aux
            cat("   ", attr(terms(object$formula),"term.labels"), "=", v.level, "\n")
            print(res[[i]], row.names = FALSE)
            cat("\n")
          }
        }else{
          ii <- sapply(times, function(x)ifelse( x >= min(object$est$y) & x <= max(object$est$y),
                                                 which.max(object$est$y[object$est$y <= x]), NA))
          if (all(is.na(ii))) {
            stop(paste("At least one element of the 'times' vector has to be between",min(object$est$y), "and", max(object$est$y)))
          }

          if (any(is.na(ii))) {
            warning(paste("'times' must be between",min(object$est$y), "and", max(object$est$y)))
          }
          res <- data.frame(times, object$est[ii, -1])
          names(res) <- names(object$est)
          print(res, row.names = FALSE)
        }
      }

      class(res) <- "summary.surv"
      return(invisible(res))

    }
  }else{
    stop("Argument x must be either survCS object.")
  }
}


