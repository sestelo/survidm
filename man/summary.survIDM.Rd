\name{summary.survIDM}
\alias{summary.survIDM}
\alias{print.survIDM}
\title{Summarizing fits of \code{survIDM} class}
\description{Returns a a data.frame or list containing the estimates of
the probabilities, its confidence limits and other information. }
\usage{
  \method{summary}{survIDM}(object, times = NULL, ...)
}

  \arguments{
\item{object}{A fitted \code{survIDM} object as produced by \code{tprob()}, \code{CIF()} or \code{sojourn()}.}
\item{times}{Vector of times; the returned data frame will contain 1 row
for each time. Missing values are not allowed.}
\item{...}{For future methods.}
}

\value{A data frame or a list containing the following components:
 \item{y}{The total time for obtaining the estimates of the probabilities.}
  \item{est}{Estimates of the probability.}
  \item{lower 95\% CI}{The lower probabilities of the interval.}
  \item{upper 95\% CI}{The upper probabilities of the interval.}
}

    \examples{
fit <- tprob(survIDM(time1,event1,Stime, event) ~ 1, s = 365,
method = "AJ", conf = TRUE, conf.level = 0.95,
conf.type = "linear", n.boot = 50, data = colonIDM)

summary(fit)
summary(fit, times = c(400, 1000, 2900))

    }


   \author{Luis Meira-Machado, Marta Sestelo and Gustavo Soutinho.}



