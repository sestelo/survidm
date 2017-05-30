\name{LLW}
\alias{LLW}

\title{Local linear weights}
\description{Computes local linear weights based on Kernel smoothing. }
\usage{
LLW(x, kernel = "gaussian", bw, t1)
}

\arguments{
  \item{x}{Covariate values for obtaining estimates for the conditional probabilities. If missing, unconditioned probabilities will be computed.}
  \item{kernel}{A character string specifying the desired kernel. See details below for possible options. Defaults to "gaussian" where the gaussian density kernel will be used.}
  \item{bw}{A single numeric value to compute a kernel density bandwidth.}
  \item{t1}{Covariate value to compute the weight at.}
}

\details{Possible options for argument window are "gaussian", "epanechnikov", "tricube", "boxcar", "triangular", "quartic" or "cosine". }
\value{A vector with local linear weights. }

\author{Luis Meira-Machado and Marta Sestelo}


\seealso{\code{\link{NWW}}}

\examples{
LLW(x = colonIDM$age, bw = 3, t1 = 60)
}

