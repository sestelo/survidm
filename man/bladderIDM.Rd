\name{bladderCS}
\alias{bladderCS}
\docType{data}
\title{Bladder Cancer Recurrences.}
\description{
  bladderCS is a data frame with 8 variables and 85 observations.
  Data on recurrences of bladder cancer, used by many people to demonstrate methodology for recurrent event modelling.}

\usage{data("bladderCS")}

\format{
  A data frame with 85 observations on the following 8 variables.
  Below a brief description is given for some of these variables.
  \describe{
    \item{t1}{Time to first recurrence/censoring, whichever occurs first.}
    \item{e1}{Recurrence/censoring indicator (first recurrence=1) for the first time (t1).}
    \item{t2}{Time to second recurrence/censoring, whichever occurs first.}
    \item{e2}{Recurrence/censoring indicator (second recurrence=1) for the second time (t2)}
    \item{t3}{Time to recurrence/censoring, whichever occurs first.}
    \item{e3}{Recurrence/censoring indicator (third recurrence=1) for the third time (t3)}
    \item{t4}{Time to fourth recurrence/censoring, whichever occurs first.}
    \item{e4}{Recurrence/censoring indicator (fourth recurrence=1) for the fourth time (t4)}
    %\item{size}{Size (cm) of largest initial tumour.}
    }
}



\references{

  Byar, D. (1980) Veterans administration study of chemoprophylaxis for recurrent stage i bladder tumors:
  Comparisons of placebo, pyridoxine and topical thiotepa. Bladder Tumors and Other Topics in
  Urological Oncology, 18:363--370.
}

\examples{
  data(bladderCS)
  head(bladderCS)
}

\keyword{datasets}
