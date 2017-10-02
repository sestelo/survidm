sojourn_ini <- function(object, t, conf = FALSE, n.boot = 199, conf.level = 0.95,
                    cluster = FALSE, ncores = NULL, method = "LM",
                    presmooth = FALSE)
{

	if (missing(object))
		stop("Argument 'object' is missing, with no default")
#	if (!inherits(object, "survIDM")) stop("'object' must be of class 'survIDM'")

	if (missing(t)) {
		t <- object[[1]]$Stime
		ptimes <- which(object[[1]]$event == 1 & object[[1]]$time1 < object[[1]]$Stime)
		t <- t[ptimes]
		           }

	if(missing(t) & method == "Datta-Satten"){
		t <- object[[1]]$Stime - object[[1]]$time1
		ptimes <- which(object[[1]]$event == 1 & object[[1]]$time1 < object[[1]]$Stime)
		t <- t[ptimes]
		t <- t[t > 0]
	}



	if (any(t <= 0)) stop("The values of 't' must be positive")
	t <- sort(unique(t))
	n <- length(t)
	if(length(t) == 0) stop("Invalid values for 't'.")

	t <- c(0,t) # esto es lo nuevo de Luís!!! para que llegue a 0
	n <- length(t)


	soj <- rep(NA, length(t))

	if (method == "LM"){
	p0 <- which(object[[1]]$time1 <  object[[1]]$Stime)
	kmw <- KMW(object[[1]]$Stime[p0], object[[1]]$event[p0])
	if (presmooth == TRUE) kmw <- PKMW(object[[1]]$Stime[p0], object[[1]]$event[p0])

	for (k in 1: length(t)) {
		p <- which(object[[1]]$Stime[p0] - object[[1]]$time1[p0] <= t[k])
		soj[k] <- sum(kmw[p])
                                     }
				     }  #method LM

	if (method == "Satten-Datta"){
		for (k in 1: length(t)) { soj[k] <- cond(object[[1]]$time1, object[[1]]$Stime, object[[1]]$event1, object[[1]]$event, t[k])}
				                }

	resu <- data.frame(cbind(t, soj))
	names(resu) <- c("t", "sojourn")

	soj.ci <- matrix(NA, length(t), 2)

    if (conf == TRUE) {
        soj.ci0 <- matrix(NA, nrow= length(t), ncol = n.boot)
        n <- dim(object[[1]])[1]

	for (j in 1:n.boot){
	  cat("", "\r")

	  cat(" Bootstrap sample number: ", j, "\r")

	xx <- sample.int(n, size = n, replace = TRUE)
        ndata <- object[[1]][xx,]
        p0 <- which(ndata$time1 < ndata$Stime)
	kmw <- KMW(ndata$Stime[p0], ndata$event[p0])
	if (presmooth == TRUE) kmw <- PKMW(ndata$Stime[p0], ndata$event[p0])

	if (method == "LM"){
	for (k in 1: length(t)) {
		p <- which(ndata$Stime[p0] - ndata$time1[p0] <= t[k])
		soj.ci0[k,j] <- sum(kmw[p])
                                     }
				     } #method LM
	if (method == "Satten-Datta"){
		for (k in 1: length(t)) { soj.ci0[k,j] <- cond(ndata$time1, ndata$Stime, ndata$event1, ndata$event, t[k])}
				                }
        }

cat("", "\r")

      for (k in 1: length(t)) {
        soj.ci[k,1] <- quantile(soj.ci0[k,], (1 - conf.level) / 2, na.rm = TRUE)
        soj.ci[k,2] <- quantile(soj.ci0[k,], 1 - (1 - conf.level) / 2, na.rm = TRUE)
				   }
        }

	   if(conf == TRUE){
	ci <- cbind(soj.ci)
	ci <- data.frame(ci)
        names(ci) <- c("soj.li.ci", "soj.ls.ci")
	                         }

	if(conf==TRUE)  result <- list(est=resu, CI=ci, conf.level=conf.level, t=t, conf=conf)
	else  result <- list(est=resu, t=t, conf=conf)

    class(result) <- c("soj")
    return(invisible(result))

}

############################################################
N <- function(time1, time, event1, status, u){
  time2 <- time - time1
  v <- which(time2 <= u & status == 1 & time2 > 0)
  if(length(v) == 0) res <- 0
  if(length(v) != 0){
	v <- sort(unique(time2[v]))
	m <- length(v)
	p <- vector(mode = "list", length = m)
	res <- vector(length = m)
	n <- length(time)
	G1<-vector(length = n)
	for (j in 1:n){ G1[j] <- KM(time, 1 - status, t = time[j])}

	for (k in 1:m) {
		   p[[k]] <- which(time2 <= v[k] & time2 > 0  & status == 1)
		   if (length(p[[k]]) == 0) res[k] <- 0
		   else res[k] <- sum(1/G1[p[[k]]])
                 }
res <- c(res[1], diff(res))}
return(res)
}


Y <- function(time1, time, event1, status, u)
{
  time2 <- time - time1
  v <- which(time2 <= u & status == 1 & time2 > 0)
  if(length(v) == 0) res <- 0
  if(length(v) != 0){
	v <- sort(unique(time2[v]))
	m <- length(v)
	p <- vector(mode = "list", length = m)
	res <- vector(length = m)
	for (k in 1:m) {
                p[[k]] <- which(time2 >= v[k] & time2 > 0)
		m1 <- length(p[[k]])
                if (m1 == 0)  res[k] <- 1
		if (m1 > 0) {
                        G1 <- vector(length = m1)
		        for (j in 1:m1){ G1[j] <- KM(time, 1 - status, t = (time1[p[[k]][j]] + v[k]))}
			res[k] <- sum(1 / G1)
                               }
                           }
			}
return(res)
}

cond <- function(time1, time, event1, status, u)
{
    D <- N(time1, time, event1, status, u)/Y(time1, time, event1, status, u)
    dd <- which(D == "NaN")
    D[dd] <- 0
    dd <- which(D == "Inf")
    D[dd] <- 0
    ddd <- which(D > 1)
    D[ddd] <- 1
    res <- 1 - prod(1 - D)
return(res)
}
