sojIPCW <- function(object, t, z.name, z.value, bw = "dpik", window = "gaussian",
            method.weights = "NW", conf = FALSE, n.boot = 199, conf.level = 0.95,
	    cluster = FALSE, ncores = NULL)

{

	if (missing(object))
		stop("Argument 'object' is missing, with no default")
#	if (!inherits(object, "survIDM")) stop("'object' must be of class 'survIDM'")
	obj <- object[[1]]

	if (missing(t)) {
		t <- obj$Stime
		ptimes <- which(obj$event == 1 & obj$time1 < obj$Stime)
		t <- t[ptimes]
		           }

	covar <- which(names(obj) == z.name)
	if (missing(z.value))
	z.value <- mean(obj[, covar])

	n <- dim(obj)[1]
	delta4 <- rep(1, n)
	lenc <- dim(obj)[2]

	if (any(t <= 0)) stop("The values of 't' must be positive")
	t <- sort(unique(t))
	m <- length(t)
	if(length(t) == 0) stop("Invalid values for 't'.")


	t <- c(0,t) # esto es lo nuevo de Luís!!! para que llegue a 0
	m <- length(t) # esto es lo nuevo de Luís!!! para que llegue a 0

	soj <- rep(NA, length(t))

	p0 <- which(obj$time1 <  obj$Stime)
	n0 <- length(p0)

	if (bw == "dpik") {
	lbd2 <- dpik(x = obj[p0, covar])
				}
	else if (bw == "np") {
		options(np.messages = FALSE)
		lbd2 <- npudensbw(dat = obj[p0, covar])$bw
				  }
	else {
		lbd2 <- bw
	}
	if (!is.numeric(bw) & !(bw %in% c("dpik", "np"))) {
		stop("Argument 'bw' have to be 'dpik', 'np' or a numeric.")
	}

        be <- rep(0, n0)
        for (k in 1:n0) {
		be[k] <- Beran(obj$Stime[p0], 1 - obj$event[p0], obj[p0, covar], delta4[p0], z.value, obj$Stime[p0][k],
                         kernel = window, bw = lbd2)
        }
        ifelse(method.weights == "NW", w1 <- NWW(obj[p0, covar], z.value, kernel = window, bw = lbd2),
               w1 <- LLW(obj[p0, covar], bw = lbd2, t1 = z.value))

	for (k in 1: length(t)) {
		p <- which(obj$Stime[p0] - obj$time1[p0] <= t[k] & obj$event[p0] == 1)
		ifelse(any(be[p] == 0), soj[k] <- NA, soj[k] <- sum(w1[p]/be[p]))

		if (soj[k] > 1) soj[k] <- 1
		if (soj[k] < 0) soj[k] <- 0
                                     }

        resu <- data.frame(cbind(t, soj))
        names(resu) <- c("t", "sojourn")

	if(missing(t)) {
	ii <- duplicated(resu$sojourn)
	t <- t[!ii]
	resu <- resu[!ii,]}

	soj.ci <- matrix(NA, length(t), 2)

    if (conf == TRUE) {
        soj.ci0 <- matrix(NA, nrow= length(t), ncol = n.boot)
        n <- dim(obj)[1]

	for (j in 1:n.boot){

	xx <- sample.int(n, size = n, replace = TRUE)
        ndata <- obj[xx,]
        p0 <- which(ndata$time1 < ndata$Stime)
	n0 <- length(p0)

	be <- rep(0, n0)
        for (k in 1:n0) {
		be[k] <- Beran(ndata$Stime[p0], 1 - ndata$event[p0], ndata[p0, covar], delta4[p0], z.value, ndata$Stime[p0][k],
                         kernel = window, bw = lbd2)
        }
        ifelse(method.weights == "NW", w1 <- NWW(ndata[p0, covar], z.value, kernel = window, bw = lbd2),
               w1 <- LLW(ndata[p0, covar], bw = lbd2, t1 = z.value))

	for (k in 1: length(t)) {
		p <- which(ndata$Stime[p0] - ndata$time1[p0] <= t[k] & ndata$event[p0] == 1)
		ifelse(any(be[p] == 0), soj.ci0[k,j] <- NA, soj.ci0[k,j] <- sum(w1[p]/be[p]))

		if(!is.na(soj.ci0[k, j])){if (soj.ci0[k, j] > 1) soj.ci0[k, j] <- 1}
                if(!is.na(soj.ci0[k, j])){if (soj.ci0[k, j] < 0) soj.ci0[k, j] <- 0}
                                     }
        }

	for (k in 1: length(t)) {
		soj.ci[k,1] <- quantile(soj.ci0[k,], (1 - conf.level) / 2, na.rm=T)
		soj.ci[k,2] <- quantile(soj.ci0[k,], 1 - (1 - conf.level) / 2, na.rm=T)
				   }
        }

	callp <- paste("soj(t|", z.name, "=", z.value, ")", sep = "")

	if(conf == TRUE){
		ci <- cbind(soj.ci)
		ci <- data.frame(ci)
		names(ci) <- c("soj.li.ci", "soj.ls.ci")
	                         }

	if (conf == FALSE) {
		result <- list(est = resu, z.name = z.name,
                     z.value = z.value, t = t, bw = bw, window = window,
                     method.weights = method.weights, conf = conf, lbd = lbd2,
		     callp = callp)
	}

	if (conf == TRUE) {
		result <- list(est = resu, CI = ci, t = t, conf.level = conf.level, z.name = z.name,
                     z.value = z.value, bw = bw, window = window, method.weights = method.weights,
                     conf = conf, lbd = lbd2, callp = callp)

    }

    class(result) <- c("IPCW", "soj")
    return(invisible(result))

}

