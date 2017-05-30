\name{PKM}
\alias{PKM}
\title{Presmoothed Kaplan-Meier product-limit estimate of survival}
\description{This function provides survival estimates using the presmoothed product-limit Kaplan-Meier estimator.}
\usage{
PKM(time, status, t)
}

\arguments{
  \item{time}{Survival time of the process.}
  \item{status}{Censoring indicator of the survival time of the process; 0 if the survival time is censored and 1 otherwise. }
  \item{t}{ The time for obtaining survival estimates.}
}

\references{
R. Cao, I. Lopez-de Ullibarri, P. Janssen, and N. Veraverbeke. Presmoothed kaplan-meier and nelsonaalen
estimators. Journal of Nonparametric Statistics, 17:31-56, 2005.

G. Dikta. On semiparametric random censorship models. Journal of Statistical Planning and Inference,
66:253-279, 1998.

E. Kaplan and P. Meier. Nonparametric estimation from incomplete observations. Journal of the
American Statistical Association, 53:457-481, 1958.
}

\author{Luis Meira-Machado and Marta Sestelo}


\seealso{\code{\link{KM}}}

\examples{
obj <- with(colonIDM, survIDM(time1, event1, Stime, event))
PKM(time = obj$Stime, status = obj$event, t = 1095)
}
