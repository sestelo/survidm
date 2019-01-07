#' @title Count number of observed transitions.
#'
#' @description Given a dataset of class "survIDM", this function counts the number of observed transitions
#' in the multi-state model.
#'
#' @param dataidm A dataframe including at least four columns named
#' \code{time1}, \code{event1}, \code{Stime} and \code{event}, which correspond
#' to disease free survival time, disease free survival indicator, time to death
#' or censoring, and death indicator, respectively.
#' @param state.names Names for the transition states. If \code{NULL} (default),
#' transition states are named by \code{"healthy"}, \code{"illness"} and \code{"death"}.
#'
#' @details The colums of the dataset needs to have the format of class "survIDM", which holds
#' the transition matrix of the multi-state model.
#'
#'
#' @examples
#'
#' nevents(colonIDM)
#'
#' nevents(colonIDM, c('State1','State2', 'State3'))
#'
#' @usage nevents(dataidm, state.names=NULL)
#'
#' @author Luis Meira-Machado, Marta Sestelo and Gustavo Soutinho.
#'
#' @references
#' L. Meira-Machado, J. de Una-Alvarez, C. Cadarso-Suarez, and P. Andersen. Multi-state models for the
#' analysis of time to event data. Statistical Methods in Medical Research, 18:195-222, 2009.
#'
#' J. de Una-Alvarez and L. Meira-Machado. Nonparametric estimation of transition probabilities in
#' the non-markov illness-death model: A comparative study. Biometrics, 71(2):364-375, 2015.
#'
#' L. Meira-Machado and M. Sestelo. Estimation in the progressive illness-death model: A nonexhaustive
#' review. Biometrical Journal, 2018.

nevents <- function (dataidm, state.names=NULL)
{

  tmat <- matrix(NA, 3, 3)
  tmat[1, 2:3] <- 1:2
  tmat[2, 3] <- 3
  if (missing(state.names))
    state.names <- c("healthy", "illness", "death")
  else {
    if (length(state.names) != 3)
      stop("incorrect length of 'state.names' argument")
  }
  dimnames(tmat) <- list(from = state.names, to = state.names)

  tmat2 <- matrix(NA, 3, 3)
  n00 <- length(which(dataidm$event1 == 0))
  n01 <- length(which(dataidm$event1 == 1 & dataidm$time1 < dataidm$Stime))
  n02 <- length(which(dataidm$event1 == 1 & dataidm$event == 1 & dataidm$time1 == dataidm$Stime))
  n12 <- length(which(dataidm$event1 == 1 & dataidm$time1 < dataidm$Stime & dataidm$event == 1))

  n0 <- length(which(dataidm$event1 == 1 & dataidm$event == 0 & dataidm$time1 == dataidm$Stime))

  if(n0>0) warning(n0, " observations with ", "Stime = time1 from state ", state.names[2], " to state ",state.names[3]  )

  colnames(tmat2) <- c(state.names)
  rownames(tmat2) <- state.names
  tmat2[1, ] <- c(n00,     n01, n02)
  tmat2[2, ] <- c(0,        n01-n12, n12)
  tmat2[3, ] <- c(0,        0,   n02+n12)

  counts <- tmat2
  class(counts) <- "table"
  return(counts)
}
