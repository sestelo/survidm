\name{markov.test}
\alias{markov.test}
\title{This function is used to test the markov assumption in the illness-death model.}

\description{The markov assumption may be tested including
the sojourn time in the initial state, "times1", and other covariates in the Cox model.}

\usage{markov.test(formula, s, nm.method = "LM", data)}

\arguments{
\item{formula}{A \code{formula} object, which must have a \code{survIDM}.}
\item{s}{The first time for obtaining a graphical test of markovianity by
comparison of the estimates for transition probabilities.}
\item{nm.method}{The non-markov method used to compute the transition probabilities. Defaults to \code{"LM"}.}
\item{data}{A data.frame including at least four columns named
\code{time1}, \code{event1}, \code{Stime} and \code{event}, which correspond
to disease free survival time, disease free survival indicator, time to death
or censoring, and death indicator, respectively.}
}

\details{The markov assumption may be tested including
the sojourn time in the initial state, "times1", and other covariates in the Cox model.
A graphical test for Markovianity is also available.}

\value{
\item{cox.markov.test}{An object of class \code{coxph} representing the fit. See \code{coxph.object} for details.}
\item{TPestimates}{Dataframe with estimates of the transition probabilities for Aalen-Johansen estimator (markovian)
and for non-markov estimator. Confidence intervals for the transition probability from State 1 to State 2 are also available.}
\item{nm.method}{The non-markov method used to compute the transition probabilities.}
\item{s}{The first time for obtaining a graphical test of markovianity by
comparison of the estimates for transition probabilities.}
\item{call}{A call object.}

}

\examples{

mk <- markov.test(survIDM(time1, event1, Stime, event) ~ 1, s = 365*4,
nm.method = "LM", data = colonIDM)
mk$cox.markov.test
mk$TPestimates
mk$nm.method
plot(mk)

}

\references{
L. Meira-Machado, J. de Una-Alvarez, C. Cadarso-Suarez, and P. Andersen. Multi-state models for the
analysis of time to event data. Statistical Methods in Medical Research, 18:195-222, 2009.

J. de Una-Alvarez and L. Meira-Machado. Nonparametric estimation of transition probabilities in
the non-markov illness-death model: A comparative study. Biometrics, 71(2):364-375, 2015.

L. Meira-Machado and M. Sestelo. Estimation in the progressive illness-death model: A nonexhaustive
review. Biometrical Journal, 2018.
}

\author{
Luis Meira-Machado, Marta Sestelo and Gustavo Soutinho.
}



