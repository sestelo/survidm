#' @title Create a survIDM object
#'
#' @description Creates a "survIDM" object, usually used as a response
#' variable in a model formula.
#'
#' @param time1 First time or censoring time.
#' @param event1 Indicator of the first time; 0 if the first time is censored
#' and 1 otherwise.
#' @param Stime The total time of the process.
#' @param event Censoring indicator of the survival time of the process; 0 if
#' the total time is censored and 1 otherwise.
#' @param ... Other options.
#'
#'
#' @details Arguments in this function must be introduced in the following
#' order: \code{time1}, \code{event1}, \code{Stime} and \code{event}, where
#' \code{time1} and \code{Stime} are ordered event times and
#' \code{event1} and \code{event} their corresponding indicator statuses.
#'
#' @return An object of class "survIDM". "survIDM" objects are implemented
#' as a single dataframe.
#'
#' @author Luis Meira-Machado and Marta Sestelo.
#'
#' @examples
#' with(colonIDM, survIDM(time1, event1, Stime, event))


survIDM <- function(time1, event1, Stime, event, ...)
{
  if (missing(time1))
    stop("Argument 'time1' is missing, with no default")
  if (missing(event1))
    stop("Argument 'event1' is missing, with no default")
  if (missing(Stime))
    stop("Argument 'Stime' is missing, with no default")
  if (missing(event))
    stop("Argument 'event' is missing, with no default")

  data <- list(time1 = as.double(time1), event1 = as.integer(event1),
               Stime = as.double(Stime), event = as.integer(event), ...)

  names(data)[1:4] <- c(substitute(time1), substitute(event1),
                        substitute(Stime), substitute(event))
  datalen <- length(data)

  attr(data, "row.names") <- as.integer(1:length(time1))
  attr(data, "class") <- "data.frame"
  #object <- list(data = na.omit(data))
  object <- na.omit(data)
  attr(object, "class") <- c("data.frame", "survIDM")
  return(object)
}

