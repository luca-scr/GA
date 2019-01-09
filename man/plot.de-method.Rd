\name{plot.de-method}
\docType{methods}
\alias{plot,de-method}

\title{Plot of Differential Evolution search path}

\description{
  The \code{plot} method for \code{\link{de-class}} objects gives a plot
  of best and average fitness values found during the iterations of the 
  DE search.
}

\usage{
\S4method{plot}{de}(x, y, ylim, cex.points = 0.7, 
     col = c("green3", "dodgerblue3",  adjustcolor("green3", alpha.f = 0.1)),
     pch = c(16, 1), lty = c(1,2), legend = TRUE, grid = graphics:::grid, \dots)
}

\arguments{
  \item{x}{An object of class \code{"ga"}.}
  \item{y}{Not used.}
  \item{ylim}{A vector of two values specifying the limits on the y-axis.}
  \item{cex.points}{The magnification to be used for points.}
  \item{col}{The colours to be used for best and average fitness values.}
  \item{pch}{The type of points to be used for best and average fitness values.}
  \item{lty}{The type of lines to be used for best and average fitness values.}
  \item{legend}{A logical specifying if a legend should be included.}
  \item{grid}{A function for grid drawing of NULL to avoid drawing one.}
  \item{\dots}{Further arguments, currently not used.}
}

\details{Plot best and average fitness values at each iteration of DE search.}
  
\value{
The method invisibly return a \code{data.frame} with the iterations and summary statistics for the fitness function evaluated at each iteration.
}

\seealso{\code{\link{de}}, \code{\link{de-class}}, \code{\link{plot,de-method}}.}

%\references{}

\author{Luca Scrucca}

\examples{
# See examples in help(de)
}

\keyword{methods}
\keyword{hplot}
