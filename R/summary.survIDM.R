
summary.survIDM <- function(object, times = NULL, ...){

  #object<-res #vem de tprob
  #times<-c(2*365,3*365,4*365)
  #object$conf<- FALSE

  if (inherits(object, "survIDM")) {

    if (class(object)[1] == "data.frame"){ # for the output of survIDM
      summary(object)
    }else{

      #---------------------------------------------------------
      #1      # para sojourn
      #------------------
      if (class(object)[1] %in%  c("soj", "sojIPCW")) {

        if (object$Nlevels > 1) {

          if(!is.null(times)){
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
          }


          cat("\n")
          cat("Estimation of", object$callp, "\n")
          cat("\n")


          res <- list()
          for (i in 1:object$Nlevels) {
            v.level <- object$levels[i]
            cat("   ", attr(terms(object$formula),"term.labels"), "=", v.level, "\n")

            if(is.null(times)) {times <- object$est[[v.level]][,1]}

            ii <- sapply(times, function(x)ifelse( x >= min(object$est[[v.level]][,1]) & x <= max(object$est[[v.level]][,1]),which.max(object$est[[v.level]][,1][object$est[[i]][,1] <= x]), NA))



            if (all(is.na(ii))) {
              aux <- data.frame(times, matrix(NA, nrow = length(times), ncol = dim(object$est[[v.level]])[2] - 1))
            }else{
              aux <- data.frame(times, object$est[[v.level]][ii, -1])
            }
            names(aux) <- c("t", "sojourn")
            print(aux, row.names = FALSE)
            res$est[[v.level]] <- aux
            if(object$conf == TRUE){
              if (all(is.na(ii))) {
                lci <- data.frame(times, matrix(NA, nrow = length(times), ncol = dim(object$est[[v.level]])[2] - 1))
                uci <- data.frame(times, matrix(NA, nrow = length(times), ncol = dim(object$est[[v.level]])[2] - 1))
              }else{
                lci <- data.frame(t = times, object$CI[[v.level]][ii, c(1)])
                uci <- data.frame(t = times, object$CI[[v.level]][ii, c(2)])
              }
              names(lci) <- c("t", "sojourn")
              names(uci) <- c("t", "sojourn")
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

          }
        }else{ # starts with no levels

          if(is.null(times)) {times <- object$est[,1]}


          ii <- sapply(times, function(x)ifelse( x >= min(object$est[,1]) & x <= max(object$est[,1]),
                                                 which.max(object$est[,1][object$est[,1] <= x]), NA))
          if (all(is.na(ii))) {
            stop(paste("At least one element of the 'times' vector has to be between",min(object$est[,1]), "and", max(object$est[,1])))
          }

          if (any(is.na(ii))) {
            warning(paste("'times' must be between",min(object$t), "and", max(object$t)))
          }


          cat("\n")
          cat("Estimation of", object$callp, "\n")
          cat("\n")

          if(object$conf == FALSE){

            if(class(object)[1] == "cifIPCW") {

              res <- list(est = data.frame(times, object$est[ii, 2]))
              names(res$est) <- c("t", "CIF")


            }else{
              res <- list(est = data.frame(times, object$est[ii, -1]))
              names(res$est) <- c("t", "sojourn")

            }



            print(res$est, row.names = FALSE)
          }else{



            if(class(object)[1] == "cifIPCW") {
              res <- list(est = data.frame(times, object$est[ii, 2]),
                          LCI = data.frame(t = times, object$LCI[ii]),
                          UCI = data.frame(t = times, object$UCI[ii]))

              names(res$est) <- c("t", "sojourn")
              names(res$LCI) <- c("t", "sojourn")
              names(res$UCI) <- c("t", "sojourn")

            }else{

              res <- list(est = data.frame(times, object$est[ii, -1]),
                          LCI = data.frame(t = times, object$CI[ii, c(1)]),
                          UCI = data.frame(t = times, object$CI[ii, c(2)]))

              names(res$est) <- c("t", "sojourn")
              names(res$LCI) <- c("t", "sojourn")
              names(res$UCI) <- c("t", "sojourn")
            }


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


        }



      } # ends sojourn


      #---------------------------------------------------------
      #2    # para CIF
      #------------------
      if (class(object)[1] %in%  c("CIF", "cifIPCW")) {


        if (class(object)[1] == "cifIPCW") object$s <- 0

        if (object$Nlevels > 1) {

          if(!is.null(times)){
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
          }


          cat("\n")
          #cat("Estimation of", object$callp, "\n")
          if (object$s == 0) {
            cat("Estimation of CIF(t)")
          }else{
            cat("Estimation of CIF(t|Y(s)=0)")
          }
          cat("\n")


          res <- list()
          for (i in 1:object$Nlevels) {
            v.level <- object$levels[i]
            cat("   ", attr(terms(object$formula),"term.labels"), "=", v.level, "\n")

            if(is.null(times)) {times <- object$est[[v.level]][,1]}

            ii <- sapply(times, function(x)ifelse( x >= min(object$est[[v.level]][,1]) & x <= max(object$est[[v.level]][,1]),which.max(object$est[[v.level]][,1][object$est[[i]][,1] <= x]), NA))



            if (all(is.na(ii))) {
              aux <- data.frame(times, matrix(NA, nrow = length(times), ncol = dim(object$est[[v.level]])[2] - 1))
            }else{
              aux <- data.frame(times, object$est[[v.level]][ii, -1])
            }
            names(aux) <- c("t", "CIF", "CIF(t|Y(s)=0)")


            if(object$s == 0) {
              aux <- aux[, -3]
            }else{
              aux <- aux[, -2]
            }

            print(aux, row.names = FALSE)
            res$est[[v.level]] <- aux
            if(object$conf == TRUE){
              if (all(is.na(ii))) {
                lci <- data.frame(times, matrix(NA, nrow = length(times), ncol = dim(object$est[[v.level]])[2] - 1))
                uci <- data.frame(times, matrix(NA, nrow = length(times), ncol = dim(object$est[[v.level]])[2] - 1))
              }else{
                lci <- data.frame(t = times, object$CI[[v.level]][ii, c(1,3)])
                uci <- data.frame(t = times, object$CI[[v.level]][ii, c(2,4)])
              }
              names(lci) <- c("t", "CIF", "CIF(t|Y(s)=0)")
              names(uci) <- c("t", "CIF", "CIF(t|Y(s)=0)")
              if(object$s == 0) {
                lci <- lci[, -3]
                uci <- uci[, -3]
              }else{
                lci <- lci[, -2]
                uci <- uci[, -2]
              }
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



          if(is.null(times)) {times <- object$est[,1]}


          ii <- sapply(times, function(x)ifelse( x >= min(object$est[,1]) & x <= max(object$est[,1]),
                                                 which.max(object$est[,1][object$est[,1] <= x]), NA))
          if (all(is.na(ii))) {
            stop(paste("At least one element of the 'times' vector has to be between",min(object$est[,1]), "and", max(object$est[,1])))
          }

          if (any(is.na(ii))) {
            warning(paste("'times' must be between",min(object$t), "and", max(object$t)))
          }


          cat("\n")

          if (object$s == 0) {
            cat("Estimation of CIF(t)")
          }else{
            cat("Estimation of CIF(t|Y(s)=0)")
          }

          cat("\n")

          if(object$conf == FALSE){

            if(class(object)[1] == "cifIPCW") {

              res <- list(est = data.frame(times, object$est[ii, 2]))
              names(res$est) <- c("t", "CIF")


            }else{
              res <- list(est = data.frame(times, object$est[ii, -1]))
              names(res$est) <- c("t", "CIF", "CIF(t|Y(s)=0)")

              if(object$s == 0) {
                res$est <- res$est[, -3]
              }else{
                res$est <- res$est[, -2]
              }

            }
            print(res$est, row.names = FALSE)
          }else{



            if(class(object)[1] == "cifIPCW") {
              res <- list(est = data.frame(times, object$est[ii, 2]),
                          LCI = data.frame(t = times, object$LCI[ii]),
                          UCI = data.frame(t = times, object$UCI[ii]))

              names(res$est) <- c("t", "CIF")
              names(res$LCI) <- c("t", "CIF")
              names(res$UCI) <- c("t", "CIF")

            }else{

              res <- list(est = data.frame(times, object$est[ii, -1]),
                          LCI = data.frame(t = times, object$CI[ii, c(1,3)]),
                          UCI = data.frame(t = times, object$CI[ii, c(2,4)]))

              names(res$est) <- c("t", "CIF", "CIF(t|Y(s)=0)")
              names(res$LCI) <- c("t", "CIF", "CIF(t|Y(s)=0)")
              names(res$UCI) <- c("t", "CIF", "CIF(t|Y(s)=0)")

              if(object$s == 0) {
                res$est <- res$est[, -3]
                res$LCI <- res$LCI[, -3]
                res$UCI <- res$UCI[, -3]
              }else{
                res$est <- res$est[, -2]
                res$LCI <- res$LCI[, -2]
                res$UCI <- res$UCI[, -2]
              }


            }


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


        }

      } # ends CIF

      #----------------------


      #---------------------------------------------------------
      #3

      # tprob
      #------------------
      if (class(object)[1] %in% c("AJ", "LIDA", "LM", "PLM", "tpIPCW",
                                  "LMAJ", "PLMAJ", "PAJ","tpBreslow")) {

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

                lci[,-1][lci[,-1]<0]<-0
                uci[,-1][uci[,-1]>1]<-1
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


              res$LCI[,-1][res$LCI[,-1]<0]<-0
              res$UCI[,-1][res$UCI[,-1]>1]<-1

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
          if (object$Nlevels > 1) {  #rx varia de 1 a 3

            # to control the times argument
            # -----------------------------
            pp <- list()

            #dim(object$est)

            for (i in 1:object$Nlevels) {
              #i<-1
              #length(object$est[[i]])
              #[,1][,1]

              object$est
              class(object$est[[i]])
              pp[[i]] <- sapply(times, function(x)
                ifelse(x >= min(object$est[[i]]) & x <= max(object$est[[i]]), 1, NA))
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

                  lci[,-1][lci[,-1]<0]<-0
                  uci[,-1][uci[,-1]>1]<-1

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

            #error due to $Obs

            ii <- sapply(times, function(x)ifelse( x >= min(object$est[,1]) & x <= max(object$est[,1]),
                                                   which.max(object$est[,1][object$est[,1] <= x]), NA))
            if (all(is.na(ii))) {
              stop(paste("At least one element of the 'times' vector has to be between",min(object$est[,1]), "and", max(object$est[,1])))
            }

            if (any(is.na(ii))) {
              warning(paste("'times' must be between",min(object$est[,1]), "and", max(object$est[,1])))
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

              res$LCI[,-1][res$LCI[,-1]<0]<-0
              res$UCI[,-1][res$UCI[,-1]>1]<-1

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
      }

      #----------------------

      class(res) <- c("summary.surv")
      return(invisible(res))
    }

  }else{
    stop("Argument x must be either survIDM object.")
  }
}


