CIF_ini <- function(object, t, s, conf = FALSE, n.boot = 199, conf.level = 0.95,
                cluster = FALSE, ncores = NULL)
{

	if (missing(object))
		stop("Argument 'object' is missing, with no default")
#	if (!inherits(object, "survIDM")) stop("'object' must be of class 'survIDM'")
  #object <- list(object)
	if (missing(t)) {
		t <- object[[1]]$time1
		ptimes <- which(object[[1]]$event1 == 1)
		t <- t[ptimes]
		           }
	if (missing(s)) s <- mean(object[[1]]$time1)

	if (any(t <= 0)) stop("The values of 't' must be positive")
	t <- sort(unique(t))
	n <- length(t)
	if(length(t) == 0) stop("Invalid values for 't'.")


	t <- c(0,t) # esto es lo nuevo de Luís!!! para que llegue a 0
	n <- length(t) # esto es lo nuevo de Luís!!! para que llegue a 0

	cif <- rep(NA, length(t))
	condcif <- rep(NA, length(t))

	kmw <- KMW(object[[1]]$time1, object[[1]]$event1)
	p0 <- which(object[[1]]$time1 > s)
	kmw1 <- KMW(object[[1]]$time1[p0], object[[1]]$event1[p0])

	for (k in 1: length(t)) {
		p <- which(object[[1]]$time1 < object[[1]]$Stime & object[[1]]$time1 < t[k])
		q <- which(object[[1]]$time1[p0] < object[[1]]$Stime[p0] & object[[1]]$time1[p0] < t[k])
		cif[k] <- sum(kmw[p])
		condcif[k] <- sum(kmw1[q])
                                     }

	resu <- data.frame(cbind(t, cif, condcif))
	names(resu) <- c("t", "CIF", "CIF(t|Y(s)=0)")

	cif.ci <- matrix(NA, length(t), 2)
	condcif.ci <- matrix(NA, length(t), 2)

    if (conf == TRUE) {
        cif.ci0 <- matrix(NA, nrow= length(t), ncol = n.boot)
	condcif.ci0 <- matrix(NA, nrow= length(t), ncol = n.boot)
        n <- dim(object[[1]])[1]

	for (j in 1:n.boot){
	xx <- sample.int(n, size = n, replace = TRUE)
        ndata <- object[[1]][xx,]
        p0 <- which(ndata$time1 > s)

	kmw <- KMW(ndata$time1, ndata$event1)
	kmw1 <- KMW(ndata$time1[p0], ndata$event1[p0])

	for (k in 1: length(t)) {
		p <- which(ndata$time1 < ndata$Stime & ndata$time1 < t[k])
		q <- which(ndata$time1[p0] < ndata$Stime[p0] & ndata$time1[p0] < t[k])
		cif.ci0[k,j] <- sum(kmw[p])
		condcif.ci0[k,j] <- sum(kmw1[q])
                                     }
        }

      for (k in 1: length(t)) {
        cif.ci[k,1] <- quantile(cif.ci0[k,], (1 - conf.level) / 2)
        cif.ci[k,2] <- quantile(cif.ci0[k,], 1 - (1 - conf.level) / 2)
        condcif.ci[k,1] <- quantile(condcif.ci0[k,], (1 - conf.level) / 2)
        condcif.ci[k,2] <- quantile(condcif.ci0[k,], 1 - (1 - conf.level) / 2)
				   }
        }

	   if(conf == TRUE){
	ci <- cbind(cif.ci, condcif.ci)
	ci <- data.frame(ci)
        names(ci) <- c("cif.li.ci", "cif.ls.ci", "condcif.li.ci", "condcif.ls.ci")

}

	if(conf==TRUE)  result <- list(est=resu, CI=ci, conf.level=conf.level, s=s, t=t, conf=conf)
	else  result <- list(est=resu,  s=s, t=t, conf=conf)

    class(result) <- c("CIF")
    return(invisible(result))


}
