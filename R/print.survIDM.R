print.survIDM <- function(x, ...){

  #x<-res

  if(inherits(x, "survIDM") & class(x)[1] =='markov'){

    #cat("Call:\n")

    #print("Markov Test")

    print(x$cox.markov.test)



  }else{


    if (inherits(x, "survIDM")) {

      if(class(x)[1] == "data.frame") {

        print(x)

      }else{


        cat("Call:\n")
        print(x$call)
        cat("\nMethod:\n")

        if(class(x)[1] == "AJ") method <- "Aalen-Johansen estimator"
        if(class(x)[1] == "LIDA") method <- "LIDA estimator"
        if(class(x)[1] == "LM") method <- "Landmark approach estimator"
        if(class(x)[1] == "PLM") method <- "Presmoothed Landmark approach estimator"
        if(class(x)[1] == "tpIPCW") method <- "Inverse Probability of Censoring Weighting for Transition Probabilities"
        if(class(x)[1] == "CIF") method <- "Cumulative Incidence Function"
        if(class(x)[1] == "cifIPCW") method <- "Inverse Probability of Censoring Weighting for the Cumulative Incidence Function"
        if(class(x)[1] == "soj") method <- "Sojourn Time Distribution"
        if(class(x)[1] == "sojIPCW") method <- "Inverse Probability of Censoring Weighting for the Sojourn Time Distribution"
        if(class(x)[1] == "LMAJ") method <- "Landmark approach Aalen-Johansen estimator"
        if(class(x)[1] == "PLMAJ") method <- "Presmoothed Landmark approach Aalen-Johansen estimator"
        if(class(x)[1] == "PAJ") method <- "Presmoothed Aalen-Johansen estimator"
        if(class(x)[1] == "PAJ") method <- "Presmoothed Aalen-Johansen estimator"
        if(class(x)[1] == "tpBreslow") method <- "Breslow Method"
        print(method)
      }

    }else{
      stop("Argument x must be either survIDM object.")
    }

  }
}
