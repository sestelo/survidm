\name{gbcsCS}
\alias{gbcsCS}
\docType{data}
\title{German Breast Cancer Study Data.}
\description{
gbcsCS is a data frame with 16 variables and 686 observations. Cancer clinical trials are a rich source for examples of applications of methods for the analysis of time to event. Willi Sauerbrei and Patrick Royston have graciously provided us with data obtained from the German Breast Cancer Study Group, which they used to illustrate methods for building prognostic models (Sauerbrei and Royston, 1999). In the main study, a total of 720 patients with primary node positive breast cancer were recruited between July 1984, and December 1989, (see Schmoor, Olschweski and Schumacher M. 1996 and Schumacher et al. (1994)).
}

\usage{data("gbcsCS")}

\format{
A data frame with 686 observations on the following 16 variables.
%Below a brief description is given for some of these variables.

%\describe{
%\item{rectime}{Time to recurrence/censoring, whichever occurs first.}
%\item{censrec}{Recurrence/censoring indicator (recurrence=1, alive=0).}
%\item{survtime}{Time to censoring/death, whichever occurs first.}
%\item{censdead}{Death/censoring indicator (death=1, alive=0).}
%\item{age}{Age in years.}
%\item{size}{Tumour size.}
%}
}

\references{
Schmoor, C., Sauerbrei, W. Bastert, G., Schumacher, M. (2000). Role of Isolated Locoregional Recurrence of Breast Cancer: Results of Four Prospective Studies. Journal of Clinical Oncology, 18(8), 1696-1708.

Schumacher, M., Bastert, G., Bojar, H., Hiibner, K., Olschewski, M., Sauerbrei, W., Schmoor, C., Beyerle, C., Neumann, R.L.A. and Rauschecker, H.F. for the German Breast Cancer Study Group (GBSG) (1994). A randomized 2 x 2 trial evaluating hormonal treatment and the duration of chemotherapy in node-positive breast cancer patients. Journal of Clinical Oncology, 12, 2086-2093.

Hosmer, D.W. and Lemeshow, S. and May, S. (2008). Applied Survival Analysis: Regression Modeling of Time to Event Data: Second Edition, John Wiley and Sons Inc., New York, NY

}

\examples{
data(gbcsCS)
head(gbcsCS)
}

\keyword{datasets}
