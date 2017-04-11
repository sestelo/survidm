print.survIDM <- function(x, ...){

  if (inherits(x, "survIDM")) {

    if(class(x)[1] == "data.frame") {

      print(x)

    }else{


      cat("Call:\n")
      print(x$call)
      cat("\nMethod:\n")

      if(class(x)[1] == "AJ") method <- "Aalen-Johansen estimator"
      if(class(x)[1] == "LIDA") method <- "LIDA estimator"
      if(class(x)[1] == "LDM") method <- "Landmark approach estimator"
      if(class(x)[1] == "PLDM") method <- "Presmoothed Landmark approach estimator"
      if(class(x)[1] == "IPCW") method <- "Inverse Probability of Censoring Weighting"

      print(method)
    }

  }else{
    stop("Argument x must be either survIDM object.")
  }
}
