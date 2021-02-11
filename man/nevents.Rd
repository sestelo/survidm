\name{nevents}
\alias{nevents}
\title{Count number of observed transitions.}

\description{Given a dataset of class "survIDM", this function counts the number of observed transitions
in the multi-state model.}

\usage{nevents(dataidm, state.names=NULL)}

\arguments{
\item{dataidm}{A dataframe including at least four columns named
\code{time1}, \code{event1}, \code{Stime} and \code{event}, which correspond
to disease free survival time, disease free survival indicator, time to death
or censoring, and death indicator, respectively.}
\item{state.names}{Names for the transition states. If \code{NULL} (default),
transition states are named by \code{"healthy"}, \code{"illness"} and \code{"death"}.
}
}

\details{The colums of the dataset needs to have the format of class "survIDM",
which holds the transition matrix of the multi-state model. The arguments \code{time1}
and \code{Stime} must be numeric and nonnegative; \code{event1} and \code{event} must
be 0 or 1 if numeric and TRUE or FALSE if logical. \code{Stime} must be greater or
equal to argument arguments \code{time1}. \code{Stime} and \code{time1} must be equal
when argument \code{event1} equals 0 or FALSE. Argument \code{event} must be equal
to 0 or FALSE when argument \code{event1} equals 0 or FALSE. When arguments
\code{Stime} and \code{time1} are equal and argument \code{event1} equals 1 or
TRUE, argument \code{event} must be equal to 1 or TRUE.}

\examples{

nevents(colonIDM)

nevents(colonIDM, c('State0','State1', 'State2'))

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



