\name{plot.ga-method}
\docType{methods}
\alias{plot,ga-method}
\alias{plot.ga}

\title{Plot of Genetic Algorithm search path}

\description{
  The \code{plot} method for \code{\link{ga-class}} objects gives a plot
  of best and average fitness values found during the iterations of the 
  GA search.
}

\usage{
\S4method{plot}{ga}(x, y, ylim, cex.points = 0.7, 
     col = c("green3", "dodgerblue3",  adjustcolor("green3", alpha.f = 0.1)),
     pch = c(16, 1), lty = c(1,2), grid = graphics:::grid, \dots)
}
\arguments{
  \item{x}{An object of class \code{"ga"}.}
  \item{y}{Not used.}
  \item{ylim}{A vector of two values specifying the limits on the y-axis.}
  \item{cex.points}{The magnification to be used for points.}
  \item{col}{The colors to be used for best and average fitness values.}
  \item{pch}{The type of points to be used for best and average fitness values.}
  \item{lty}{The type of lines to be used for best and average fitness values.}
  \item{grid}{A function for grid drawing of NULL to avoid drawing one.}
  \item{\dots}{Further arguments, currently not used.}
}

\details{Plot best and average fitness values at each iteration of GA search.}
  
\value{
The method invisibly return a list with the following components:

\item{iter}{a vector of values for each interation.}
\item{fitnessBest}{the best value of fitness function at each iteration.}
\item{fitnessMean}{the mean value of fitness function at each iteration.}
}

\seealso{\code{\link{ga}}, \code{\link{ga-class}}.}

%\references{}

\author{Luca Scrucca}

\keyword{methods}
\keyword{hplot}
