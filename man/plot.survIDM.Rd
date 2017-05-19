\name{plot.survIDM}
\alias{plot.survIDM}
\title{Plot for an object of class "survIDM".}
\description{It draws the estimated probabilities. }
\usage{
\method{plot}{survIDM}(x = object, y = NULL, trans = "all", conf = NULL,
type = NULL, conftype = NULL, col = 1:6, confcol = 1:6, lty = 1, conflty = 2,
xlab = "Time (years)", ylab = NULL, ylim = NULL, xlim = NULL, ...)
}

\arguments{
  \item{x}{An object of class "survIDM".}
  \item{y}{\code{NULL}}
  \item{trans}{The transition probabilities plotted. It is used only when the object is of class "AJ", "LIDA", "LDM" or "PLDM". Possible options are \code{"all"} (default), \code{"00"}, \code{"01"}, \code{"02"}, \code{"11"} or \code{"12"}.}
  \item{conf}{Draw the confidence intervals into the plot. By default it is
  \code{NULL}, they are drawn if the "surv" object contains them.}
  \item{type}{The type of plot that should be drawn. See details \code{\link{par}}
  for possible options. Defaults to "s" for the draw be stair steps.}
  \item{conftype}{The type of plot that should be drawn for confidence intervals.
  See details \code{\link{par}}
  for possible options. Defaults to "s" for the draw be stair steps.}
  \item{col}{Vector of colors. Colors are used cyclically.}
  \item{confcol}{Vector of colors for the confidence intervals.
  Colors are used cyclically.}
  \item{lty}{The line type. Line types can either be specified as an integer
  (0 = blank, 1 = solid (default), 2 = dashed, 3 = dotted, 4 = dotdash,
  5 = longdash, 6 = twodash).  See details in \code{\link{par}}.}
  \item{conflty}{The line type for confidence intervals. Line types can either
  be specified as an integer (0 = blank, 1 = solid (default), 2 = dashed,
  3 = dotted, 4 = dotdash, 5 = longdash, 6 = twodash).}
  \item{xlab}{A title for the \code{x} axis: see \code{\link{title}}.}
  \item{ylab}{A title for the \code{y} axis: see \code{\link{title}}.}
  \item{ylim}{The \code{y} limits of the plot.}
  \item{xlim}{The \code{x} limits of the plot.}
  \item{\dots}{Other options.}
}

\value{No value is returned.}

\author{Luis Meira-Machado and Marta Sestelo}


\examples{
fit <- tprob(survIDM(time1,event1,Stime, event) ~ 1, s = 365,
             method = "AJ", conf = TRUE, conf.level = 0.95,
             conf.type = "linear", n.boot = 50, data = colonCS)

plot(fit)
plot(fit, trans = "11")


fit2 <- CIF(survIDM(time1,event1,Stime, event) ~ factor(sex), data = colonCS,
            conf = FALSE, conf.level = 0.95)

plot(fit2)


}
