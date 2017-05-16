\name{Beran}
\alias{Beran}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Estimation of the conditional distribution function of the response, given the covariate under random censoring.}
\description{Computes the conditional survival probability P(T > y|Z = z)}
\usage{ Beran(time, status, covariate, delta, x, y, kernel = "gaussian", bw, lower.tail = FALSE)}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{time}{The survival time of the process.}
  \item{status}{Censoring indicator of the total time of the process; 0 if the total time is censored and 1 otherwise.}
  \item{covariate}{Covariate values for obtaining estimates for the conditional probabilities.}
  \item{delta}{Censoring indicator of the covariate.}
  \item{x}{The first time (or covariate value) for obtaining estimates for the conditional probabilities. If missing, 0 will be used.}
  \item{y}{The total time for obtaining estimates for the conditional probabilities.}
  \item{kernel}{A character string specifying the desired kernel. See details below for possible options. Defaults to "gaussian" where the gaussian density kernel will be used.}
  \item{bw}{A single numeric value to compute a kernel density bandwidth.}
  \item{lower.tail}{logical; if FALSE (default), probabilities are P(T > y|Z = z) otherwise, P(T <= y|Z = z).}
}

\details{ Possible options for argument window are "gaussian", "epanechnikov", "tricube", "boxcar", "triangular", "quartic" or "cosine".}

\references{
R. Beran. Nonparametric regression with randomly censored survival data. Technical report, University of
California, Berkeley, 1981.
}

\author{Luis Meira-Machado and Marta Sestelo}

\examples{
obj <- with(colonCS, survCS(time1, event1, Stime, event))
obj0 <- obj$data

#P(T>y|age=45)
library(KernSmooth)
h <- dpik(colonCS$age)
Beran(time = obj0$Stime, status = obj0$event, covariate = colonCS$age,
x = 45, y = 730, bw = h)

#P(T<=y|age=45)
Beran(time = obj0$Stime, status = obj0$event, covariate = colonCS$age,
x = 45, y = 730, bw = h, lower.tail = TRUE)
}

