plot.survIDM <- function(x = object, y = NULL, trans = "all", func = "distribution",
                         conf = NULL, type = NULL,conftype = NULL, col = 1:6,
                         confcol = 1:6, lty = 1, conflty = 2, xlab = "Time (years)",
                         ylab = NULL, ylim = NULL, xlim = NULL, ...) {

  object <- x

  if (inherits(x, "survIDM")) {


    if (class(x)[1] == "data.frame") {
      plot(x)
    }


    if (!class(object)[1] %in% c("AJ", "LIDA", "LM", "PLM", "tpIPCW", "CIF",
                                 "cifIPCW", "soj", "sojIPCW", "LMAJ",
                                 "PLMAJ", "PAJ")) {
      stop("The argumment 'Object' must be of one of the following classes
           'AJ', 'LIDA', 'LM', 'PLM', 'LMAJ', 'PLMAJ', 'PAJ', 'tpIPCW',
           'CIF', 'cifIPCW', 'soj' or 'sojIPCW'")
    }



    # for all
    #-----------------------------------------------

    object <- x

    if (object$Nlevels != length(col))
      col <- rep(col, times = object$Nlevels)
    if (object$Nlevels != length(confcol))
      confcol <- rep(confcol, times = object$Nlevels)


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


    #---------------------------






      #tp
      #----------------------

      if (class(object)[1] %in% c("AJ", "LIDA", "LM", "PLM", "tpIPCW",
                                  "LMAJ", "PLMAJ", "PAJ")) {

        if(is.null(ylab) & class(object)[1] != "tpIPCW")
          ylab <- bquote(paste(p[ij], "(", .(x$s), ",t)"))

        if(is.null(ylab) & class(object)[1] == "tpIPCW")
          ylab <- bquote(paste(p[ij], "(", .(x$s), ",t|", .(x$z.name),")"))

        #-------------
        trans2 = trans
        tp <- c("00", "01", "02", "11", "12")
        if(trans == "all") {trans2 = tp}
        ii <- trans2 == tp
        itp <- 2:6
        itpCI <- c(1, 3, 5, 7, 9)
        #--------------

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
      } #end plot for tp



    #cif
      #---------------------

      if (class(object)[1] %in%  c("CIF", "cifIPCW")) {
        if (is.null(ylab) &  class(object)[1] == "CIF") ylab <- "CIF(t)"

        if (is.null(ylab) &  class(object)[1] == "cifIPCW")
          ylab <- bquote(paste("CIF(t|", .(x$z.name), ")"))


        if(class(object)[1] == "cifIPCW") object$s <- 0


        if (object$Nlevels == 1) {

          if(object$s != 0){
            ob <- ob[, -2]
            #obCI <- obCI[, -c(1:2)]
          }


          if(class(object)[1] == "cifIPCW" & ci == TRUE) {
            obCI <- ob[, 3:4] # in order to corerct the out of cifIPCW
          }


          matplot(ob[, 1], ob[, 2], type = type, col = col, xlab = xlab,
                  ylab = ylab, lty = lty, ylim = ylim, xlim = xlim, ...)

          if (ci == TRUE) {

            if(object$s != 0){
              #ob <- ob[, -2]
              obCI <- obCI[, -c(1:2)]
            }

            matlines(x = ob[, 1], y = obCI[, 1], type = conftype,
                     lty = conflty, col = confcol, ...)
            matlines(x = ob[, 1], y = obCI[, 2], type = conftype,
                     lty = conflty, col = confcol, ...)
          }

        }else{ # more than 1 level


          plot(ob[[1]][, 1], ob[[1]][, 2], type = "n", xlab = xlab,
               ylab = ylab, ylim = ylim, xlim = xlim, ...)

          for (i in 1:object$Nlevels) {

            if(object$s != 0){
              ob [[i]]<- ob[[i]][, -2]
             # obCI[[i]] <- obCI[[i]][, -c(1:2)]
            }

            lines(ob[[i]][, 1], ob[[i]][, 2], type = type, col = col[i],
                  lty = lty, ...)
            if (ci == TRUE) {
              if(object$s != 0){
                #ob [[i]]<- ob[[i]][, -2]
                obCI[[i]] <- obCI[[i]][, -c(1:2)]
              }

              lines(x = ob[[i]][, 1], y = obCI[[i]][, 1], type = conftype,
                    lty = conflty, col = confcol[i], ...)
              lines(x = ob[[i]][, 1], y = obCI[[i]][, 2], type = conftype,
                    lty = conflty, col = confcol[i], ...)
              }
          }
          legend("bottomright", object$levels, col = col, lty = lty)



        }



      } #ends for CIF



    # soj
      #--------------------------------

      if (class(object)[1] %in%  c("soj", "sojIPCW")) {
        if (is.null(ylab) &  class(object)[1] == "soj") ylab <- "Sojourn(t)"

        if (is.null(ylab) &  class(object)[1] == "sojIPCW")
          ylab <- bquote(paste("Sojourn(t|", .(x$z.name), ")"))





        if (object$Nlevels == 1) {


       #   if(class(object)[1] == "cifIPCW") obCI <- ob[, 3:4] # in order to corerct the out of cifIPCW
          if(func == "survival"){ob[, 2] <- 1 - ob[, 2]}
          if(func == "survival" & ci == TRUE){obCI <- 1 - obCI}

          matplot(ob[, 1], ob[, 2], type = type, col = col, xlab = xlab,
                  ylab = ylab, lty = lty, ylim = ylim, xlim = xlim, ...)

          if (ci == TRUE) {
            matlines(x = ob[, 1], y = obCI[, 1], type = conftype,
                     lty = conflty, col = confcol, ...)
            matlines(x = ob[, 1], y = obCI[, 2], type = conftype,
                     lty = conflty, col = confcol, ...)
          }

        }else{ # more than 1 level


          plot(ob[[1]][, 1], ob[[1]][, 2], type = "n", xlab = xlab,
               ylab = ylab, ylim = ylim, xlim = xlim, ...)

          for (i in 1:object$Nlevels) {

            if(func == "survival"){ob[[i]][, 2] <- 1 - ob[[i]][, 2]}
            if(func == "survival" & ci == TRUE){obCI[[i]] <- 1 - obCI[[i]]}

            lines(ob[[i]][, 1], ob[[i]][, 2], type = type, col = col[i],
                  lty = lty, ...)
            if (ci == TRUE) {
              lines(x = ob[[i]][, 1], y = obCI[[i]][, 1], type = conftype,
                    lty = conflty, col = confcol[i], ...)
              lines(x = ob[[i]][, 1], y = obCI[[i]][, 2], type = conftype,
                    lty = conflty, col = confcol[i], ...)
            }
          }
          legend("bottomright", object$levels, col = col, lty = lty)



        }








      } #ends for sojourn




  }else{
    stop("Argument x must be either survIDM object.")
  }
}







