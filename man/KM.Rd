\name{KM}
\alias{KM}

\title{Kaplan-Meier product-limit estimate of survival}

\description{This function provides survival estimates using the product-limit Kaplan-Meier estimator.}

\usage{
KM(time, status, t)
}

\arguments{
  \item{time}{ Survival time of the process.}
  \item{status}{Censoring indicator of the survival time of the process; 0 if the survival time is censored and 1 otherwise. }
  \item{t}{ The time for obtaining survival estimates.}
}


\references{
E. Kaplan and P. Meier. Nonparametric estimation from incomplete observations. Journal of the
American Statistical Association, 53:457-481, 1958.
}

\author{Luis Meira-Machado, Marta Sestelo and Gustavo Soutinho.}


\seealso{\code{\link{PKM}}}

\examples{
require(survival)
obj <- with(colonIDM, survIDM(time1, event1, Stime, event))
KM(time = obj$Stime, status = obj$event, t = 1095)

fit <- survfit(Surv(obj$Stime, obj$event) ~ 1, data = obj)
summary(fit, time = 1095)$surv

}

