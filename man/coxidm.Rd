\name{coxidm}
\alias{coxidm}
\alias{coxidm}
\title{Fit proportional hazards regression model in each transition of
  the Illness-Death Model.}
\description{Fits a Cox proportional hazards regression model for each transition.}

\usage{coxidm(formula, data, semiMarkov = FALSE)
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
  \item{semiMarkov}{By default, \code{FALSE}.}
}

\value{A data frame or a list containing the following components:
  \item{coef}{Estimated coeficients.}
  \item{term}{Predicted values of splines and partial residual for the terms of the model.}
  \item{exp(coef)}{Exponent of the estimated coefficients.}
  \item{lower 0.95}{Lower limit of the confidence interval.}
  \item{upper 0.95}{Upper limit of the confidence interval.}
  \item{pvalue}{obtained pvalue testing that the coefficient is equals to zero.}
}

\examples{
  cmm1 <- coxidm(survIDM(time1, event1, Stime, event) ~ age, data = colonIDM)
  summary(cmm1, conf.level = 0.95)

  cmm2 <- coxidm(survIDM(time1, event1, Stime, event) ~ rx + sex + age + nodes,
                 data = colonIDM)
  summary(cmm2)

  cmm3 <- coxidm(survIDM(time1, event1, Stime, event) ~ rx + sex + age + nodes,
                 data = colonIDM, semiMarkov = TRUE)
  summary(cmm3)
}


\author{Luis Meira-Machado, Marta Sestelo and Gustavo Soutinho.}



