\name{PKMW}
\alias{PKMW}
\title{Presmoothed Kaplan-Meier weights}
\description{This function returns a vector with the presmoothed Kaplan-Meier weights.}
\usage{
PKMW(time, status)
}

\arguments{

  \item{time}{ Survival time of the process.}
  \item{status}{Censoring indicator of the survival time of the process; 0 if the survival time is censored and 1 otherwise. }
}

\value{ Vector with presmoothed Kaplan-Meier weights. }

\references{
R. Cao, I. Lopez-de Ullibarri, P. Janssen, and N. Veraverbeke. Presmoothed kaplan-meier and nelsonaalen
estimators. Journal of Nonparametric Statistics, 17:31-56, 2005.

G. Dikta. On semiparametric random censorship models. Journal of Statistical Planning and Inference,
66:253-279, 1998.

E. Kaplan and P. Meier. Nonparametric estimation from incomplete observations. Journal of the
American Statistical Association, 53:457-481, 1958.
}

\author{Luis Meira-Machado and Marta Sestelo}


\seealso{\code{\link{KMW}}}

\examples{
obj <- with(colonIDM, survIDM(time1, event1, Stime, event))
PKMW(time = obj$Stime, status = obj$event)
}

