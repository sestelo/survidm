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

