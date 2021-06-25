\name{summary.cmm}
\alias{summary.cmm}
\alias{print.cmm}
\title{Summarizing fits of \code{cmm} class}
\description{Produces a summary of a fitted coxidm model (proportional
hazards regression model in each transition of the Illness-Death Model.}
\usage{
  \method{summary}{cmm}(object, type=NULL, conf.level = 0.95, ...)
}

\arguments{
  \item{object}{A fitted \code{cmm} object as produced by \code{coxidm()}.}
  \item{type}{Type of summary of the Cox model. Default value is \code{"NULL"}
  that produces a summary of a fitted coxidm model. Possible choices are
  \code{"anova"} and \code{"ph"}. The \code{"anova"} option provides the
  Analysis of Deviance of the regression model. The \code{"ph"} option computes
  the tests of the proportional hazards assumption of the Cox Regression models.
  }
  \item{conf.level}{Level for computation of the confidence intervals. If set
  to \code{FALSE} no confidence intervals are printed.}
  \item{\dots}{Other options.}
}

\value{A data frame or a list containing the following components:
  \item{coef}{Estimated coeficients.}
  \item{exp(coef)}{Exponent of the estimated coefficients.}
  \item{lower 0.95}{Lower limit of the confidence interval.}
  \item{upper 0.95}{Upper limit of the confidence interval.}
  \item{pvalue}{obtained pvalue testing that the coefficient is equals to zero.}

  Other data frames with the Analysis of Deviance or the tests of the
  Proportional Hazards Assumption
  of the Cox Regression models.
}

\examples{
\donttest{
  cmm1 <- coxidm(survIDM(time1, event1, Stime, event) ~ age, data = colonIDM)
  summary(cmm1, conf.level = 0.95)

  cmm2 <- coxidm(survIDM(time1, event1, Stime, event) ~ rx + sex + age + nodes,
                 data = colonIDM)
  summary(cmm2)

  summary(cmm2,type = 'anova')

  cmm3 <- coxidm(survIDM(time1, event1, Stime, event) ~ rx + sex + age + nodes,
                 data = colonIDM, semiMarkov = TRUE)
  summary(cmm3)

  cmm4 <- coxidm(survIDM(time1, event1, Stime, event) ~ rx + sex + age +
                pspline(nodes) + surg + adhere, data = colonIDM)

  summary(cmm4, type = 'ph')
  }
}

\author{Luis Meira-Machado, Marta Sestelo and Gustavo Soutinho.}



