% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sojourn.R
\name{sojourn}
\alias{sojourn}
\title{Nonparametric estimation of the Sojourn time distributions in the
recurrence state in the illness-death model}
\usage{
sojourn(
  formula,
  data,
  conf = FALSE,
  n.boot = 199,
  conf.level = 0.95,
  z.value,
  bw = "dpik",
  window = "gaussian",
  method.weights = "NW",
  method = "LM",
  presmooth = FALSE,
  cluster = FALSE,
  ncores = NULL
)
}
\arguments{
\item{formula}{A \code{formula} object, which must have a \code{survIDM}
object as the response on the left of the \code{~} operator and, if desired,
a term on the right. The term may be a qualitative or quantitative variable.
Without covariates, the right hand side should be \code{~ 1}.}

\item{data}{A data.frame including at least four columns named
\code{time1}, \code{event1}, \code{Stime} and \code{event}, which correspond
to disease free survival time, disease free survival indicator, time to death
or censoring, and death indicator, respectively.}

\item{conf}{Provides pointwise confidence bands. Defaults to \code{FALSE}.}

\item{n.boot}{The number of bootstrap replicates to compute the variance
of the non-Markovian estimator. Default is 199.}

\item{conf.level}{Level of confidence. Defaults to 0.95 (corresponding to 95\%).}

\item{z.value}{The value of the covariate on the right hand side of formula
at which the sojourn probabilities are computed. For quantitative
covariates, i.e. of class integer and numeric.}

\item{bw}{A single numeric value to compute a kernel density bandwidth.
Use \code{"dpik"} for the \pkg{KernSmooth} package based selector or \code{"np"}
for the \code{'npudensbw'} function of the \pkg{np} package.}

\item{window}{A character string specifying the desired kernel.
See details below for possible options. Defaults to \code{"gaussian"}
where the gaussian density kernel will be used.}

\item{method.weights}{A character string specifying the desired weights method.
Possible options are \code{"NW"} for the Nadaraya-Watson weights and \code{"LL"}
for local linear weights. Defaults to \code{"NW"}.}

\item{method}{The method used to compute the sojourn estimates.
Possible options are \code{"LM"} and \code{"Satten-Datta"}.
Defaults to \code{"LM"}.}

\item{presmooth}{- A logical value. If \code{TRUE}, the presmoothed landmark
estimator of the sojourn function is computed. Only valid for \code{method = "LM"}.}

\item{cluster}{A logical value. If \code{TRUE} (default), the bootstrap procedure
for the confidence intervals is parallelized. Note that there are
cases (e.g., a low number of bootstrap repetitions) that \R will gain in
performance through serial computation. \R takes time to distribute tasks
across the processors also it will need time for binding them all together
later on. Therefore, if the time for distributing and gathering pieces
together is greater than the time need for single-thread computing,
it does not worth parallelize.}

\item{ncores}{An integer value specifying the number of cores to be used in
the parallelized procedure. If \code{NULL} (default), the number of cores
to be used is equal to the number of cores of the machine - 1.

#'}
}
\value{
An object of class \code{"survIDM"} and one of the following
two classes: \code{"soj"} (Sojourn Time Distribution), and
\code{"sojIPCW"} (Inverse Probability of Censoring Weighting for the Sojourn Time Distribution). Objects are implemented as a list with elements:

\item{est}{data.frame with estimates of the sojourn probabilities.}
\item{CI}{data.frame with the confidence intervals of the sojourn
probabilities.}
\item{conf.level}{Level of confidence.}
\item{t}{The time for obtaining the estimates of sojourn probabilities.}
\item{conf}{logical; if \code{FALSE} (default) the pointwise confidence
bands are not given.}
\item{callp}{The expression of the estimated probability.}
\item{Nlevels}{The number of levels of the covariate. Provides important
information when the covariate at the right hand side of formula
is of class factor.}
\item{levels}{The levels of the qualitative covariate
(if it is of class factor) on the right hand side of formula.}
\item{formula}{A formula object.}
\item{call}{A call object.}
}
\description{
This function is used to obtain nonparametric estimates of
of the sojourn probabilities in the recurrence state in the illness-death
model.
}
\details{
Possible options for argument window are \code{"gaussian"},
\code{"epanechnikov"}, \code{"tricube"}, \code{"boxcar"},
\code{"triangular"}, \code{"quartic"} or \code{"cosine"}.
}
\examples{

res <- sojourn(survIDM(time1, event1, Stime, event) ~ 1,
data = colonIDM, conf = FALSE, conf.level = 0.95)
res
summary(res, time=365*1:6)
plot(res)

res1 <- sojourn(survIDM(time1, event1, Stime, event) ~ 1,
data = colonIDM, conf = FALSE, conf.level = 0.95, method = "LM",
presmooth = TRUE)
res1
summary(res1, time=365*1:6)
plot(res1)


# not run:
#res2 <- sojourn(survIDM(time1, event1, Stime, event) ~ 1,
#data = colonIDM, conf = FALSE, conf.level = 0.95, method = "Satten-Datta")
#res2


# with a factor
res3 <- sojourn(survIDM(time1, event1, Stime, event) ~ factor(sex),
data = colonIDM, conf = FALSE, conf.level = 0.95)
res3
summary(res3, time=365*1:6)
plot(res3)

# with a qualitative covariate
res4 <- sojourn(survIDM(time1, event1, Stime, event) ~ age, data = colonIDM,
z.value = 56, conf = FALSE, conf.level = 0.95)
res4
summary(res4, time=365*1:6)
plot(res4)

}
\references{
Satten, G.A. and Datta, S. (2002) Marginal estimation for
multi-stage models: waiting time distributions and competing risks
analyses. Statistics in Medicine, 21, 3--19.
}
\author{
Luis Meira-Machado, Marta Sestelo and Gustavo Soutinho.
}
