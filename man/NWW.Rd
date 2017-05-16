\name{NWW}
\alias{NWW}

\title{Nadaraya-Watson weights.}
\description{Computes the Nadaraya-Watson weights. }
\usage{
NWW(covariate, x, kernel = "gaussian", bw)
}

\arguments{
  \item{covariate}{Covariate values for obtaining weights.}
  \item{x}{Covariate value to compute the weight at.}
  \item{kernel}{A character string specifying the desired kernel. See details below for possible options. Defaults to "gaussian" where the gaussian density kernel will be used.}
  \item{bw}{A single numeric value to compute a kernel density bandwidth.}
}
\details{Possible options for argument window are "gaussian", "epanechnikov", "tricube", "boxcar", "triangular", "quartic" or "cosine". }
\value{A vector with Nadaraya-Watson weights. }

\author{Luis Meira-Machado and Marta Sestelo}

\seealso{\code{\link{LLW}}}

\examples{
NWW(covariate = colonCS$age, x=40, kernel = "gaussian", bw = 3)
}
