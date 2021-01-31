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
#' \code{time1} and \code{Stime} are the sojourn time in the initial state and
#' the total time, respectively. \code{event1} and \code{event} denote their
#' corresponding indicator statuses. This function checks the following
#' conditions: (i) the arguments \code{time1} and \code{Stime} must be numeric
#' and nonnegative; \code{event1} and \code{event} must be 0 or 1 if numeric
#' and TRUE or FALSE if logical. \code{Stime} must be greater or equal to
#' argument arguments \code{time1}. \code{Stime} and \code{time1} must be
#' equal when argument \code{event1} equals 0 or FALSE. Argument \code{event}
#' must be equal to 0 or FALSE when argument \code{event1} equals 0 or FALSE.
#' When arguments \code{Stime} and \code{time1} are equal and argument
#' \code{event1} equals 1 or TRUE, argument \code{event} must be equal to
#' 1 or TRUE.
#'
#' @return An object of class "survIDM". "survIDM" objects are implemented
#' as a single dataframe.
#'
#' @author Luis Meira-Machado, Marta Sestelo and Gustavo Soutinho.
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
  if ( !is.numeric(time1) ) stop("Argument 'time1' is not numeric")
  if ( !( is.logical(event1) | is.numeric(event1) ) ) stop("Argument event1 must be logical or numeric")
  if ( !is.numeric(Stime) ) stop("Argument 'Stime' is not numeric")
  if ( !( is.logical(event) | is.numeric(event) ) ) stop("Argument event must be logical or numeric")
  len <- length(time1)
  if ( len != length(event1) | len != length(Stime) | len != length(event) ) stop("Arguments 'time1', 'event1', 'Stime' and 'event' must have the same length")

  if ( any( (event1 != 0 & event1 != 1) | (event1 != FALSE & event1 != TRUE) ) ) stop("Argument 'event1' must be 0 or 1 if numeric and TRUE or FALSE if logical")
  if ( any( (event != 0 & event != 1) | (event != FALSE & event != TRUE) ) ) stop("Argument 'event' must be 0 or 1 if numeric and TRUE or FALSE if logical")
  if ( any(time1 < 0 | Stime < 0) ) stop("Arguments 'time1' and 'Stime' must be greater or equal 0")
  if ( any(Stime < time1) ) stop("Argument 'Stime' must be greater or equal to argument 'time1'")
  if ( any(!event1 & Stime != time1) ) stop("Arguments 'Stime' and 'time1' must be equal when argument 'event1' equals 0 or FALSE")
  if ( any(!event1 & event) ) stop("Argument 'event' must be equal to 0 or FALSE when argument 'event1' equals 0 or FALSE")
  if ( any(time1 == Stime & event1 & !event) ) {
    pos<-which(time1==Stime & event1==1 & event==0)
    cat("Stime=time1 and event1=event=1 for the following observations","\n")
    print(pos)}
  if ( any(time1 == Stime & event1 & !event) ) {stop("When arguments 'Stime' and 'time1' are equal and argument 'event1' equals 1 or TRUE, argument 'event' must equal 1 or TRUE")}

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
