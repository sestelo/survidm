\name{coxidm}
\alias{coxidm}
\alias{coxidm}
\title{Fit proportional hazards regression model in each transition of
  the Illness-Death Model.}
\description{Fits a Cox proportional hazards regression model for each transition.}
\usage{coxidm(object, conf.int = 0.95)
}

\arguments{
  \item{object}{A fitted \code{cmm} object as produced by \code{coxidm()}.}
  \item{conf.int}{level for computation of the confidence intervals.
    If set to \code{FALSE} no confidence intervals are printed.}
}

\value{A data frame or a list containing the following components:
    \item{coef}{Estimated coeficients.}
  \item{exp(coef)}{Exponent of the estimated coefficients.}
  \item{lower 0.95}{Lower limit of the confidence interval.}
  \item{upper 0.95}{Upper limit of the confidence interval.}
  \item{pvalue}{obtained pvalue testing that the coefficient is equals to zero.}
}

\examples{
  cmm1 <- coxidm(survIDM(time1, event1, Stime, event) ~ age, object = colonIDM)
  summary(cmm1, conf.level = 0.95)

  cmm2 <- coxidm(survIDM(time1, event1, Stime, event) ~ rx + sex + age + nodes,
                 object = colonIDM)
  summary(cmm2)

  cmm3 <- coxidm(survIDM(time1, event1, Stime, event) ~ rx + sex + age + nodes,
                 object = colonIDM, semiMarkov = TRUE)
  summary(cmm3)
}


\author{Luis Meira-Machado and Marta Sestelo.}



