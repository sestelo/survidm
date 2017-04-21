plot.survIDM <- function(x = object, y = NULL, trans = "all", conf = NULL, type = NULL,
                        conftype = NULL, col = 1:6, confcol = 1:6, lty = 1, conflty = 2,
                        xlab = "Time", ylab = "Survival", ylim = NULL, xlim = NULL, ...) {

  if (inherits(x, "survIDM")) {

    if (class(x)[1] == "data.frame") {
      plot(x)

    }else{


    #tp
    #----------------------
      trans2 = trans
      tp <- c("00", "01", "02", "11", "12")
      if(trans == "all") {trans2 = tp}
      ii <- trans2 == tp
      itp <- 2:6
      itpCI <- c(1, 3, 5, 7, 9)
    #--------------

      object <- x

      if (object$Nlevels != length(col))
        col <- rep(col, times = object$Nlevels)
      if (object$Nlevels != length(confcol))
        confcol <- rep(confcol, times = object$Nlevels)

      if (class(object)[1] != "AJ" & class(object)[1] != "LIDA" &
          class(object)[1] != "LDM" & class(object)[1] != "PLDM" &
          class(object)[1] != "IPCW") {
        stop("The argumment 'Object' must be of one of the following classes
     'AJ', 'LIDA', 'LDM', 'PLDM', 'IPCW'")
      }



      if (is.null(type))
        type <- "s"
      if (is.null(conftype))
        conftype <- "s"


      if (is.null(conf)) {
        ci <- object$conf
      } else {
        if (conf == TRUE & object$conf == FALSE) {
          stop("The surv object does not contain confidence intervals")
        }
        if (conf == TRUE & object$conf == TRUE)
          ci <- TRUE
        if (conf == FALSE)
          ci <- FALSE

      }



      if (is.null(ylim)) ylim <- c(0, 1)

      ob <- object$est
      obCI <- object$CI

      if (is.null(xlim) & object$Nlevels > 1) {
        xlim <- c(min(sapply(ob, function(x) min(x[, 1]), simplify = TRUE)),
                  max(sapply(ob, function(x) max(x[, 1]), simplify = TRUE)))
      }


      #if (is.null(ylim)) ylim <- c(0,1)



       if (object$Nlevels == 1) {

        matplot(ob[, 1], ob[, itp[ii]], type = type, col = col, xlab = xlab,
             ylab = ylab, lty = lty, ylim = ylim, xlim = xlim, ...)

        if (ci == TRUE) {
          matlines(x = ob[, 1], y = obCI[, itpCI[ii]], type = conftype,
                lty = conflty, col = confcol, ...)
          matlines(x = ob[, 1], y = obCI[, itpCI[ii] + 1], type = conftype,
                lty = conflty, col = confcol, ...)
        }

       if(trans == "all") legend("topright", c("00", "01", "02", "11", "12")[itp - 1], col = col, lty = lty)

      } else { # more than 1 level

        if (trans == "all") {
          stop(paste("The argumment 'trans' can't be 'all' if the factor", attr(terms(object$formula),"term.labels"),
                                        "is included in the formula, you must select one of the transition probabilities."))
        }


        plot(ob[[1]][, 1], ob[[1]][, 2], type = "n", xlab = xlab,
             ylab = ylab, ylim = ylim, xlim = xlim, ...)

        for (i in 1:object$Nlevels) {


          lines(ob[[i]][, 1], ob[[i]][, itp[ii]], type = type, col = col[i],
                lty = lty, ...)
          if (ci == TRUE) {
            lines(x = ob[[i]][, 1], y = obCI[[i]][, itpCI[ii]], type = conftype,
                  lty = conflty, col = confcol[i], ...)
            lines(x = ob[[i]][, 1], y = obCI[[i]][, itpCI[ii] + 1], type = conftype,
                  lty = conflty, col = confcol[i], ...)

          }
        }
        legend("topright", object$levels, col = col, lty = lty)
      }
    }
  }else{
    stop("Argument x must be either survCS object.")
  }
}






