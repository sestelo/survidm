PKM <-
function(time,status,t)
{
 if (missing(time)) stop("Argument 'time' is missing with no default")
 if (missing(status)) stop("Argument 'status' is missing with no default")
 if (missing(t)) stop("Argument 't' is missing with no default")
p <- which(time <= t)
res <- 1 - sum(PKMW(time,status)[p])
return(res)
}
