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
     pch = c(16, 1), lty = c(1,2), legend = TRUE, grid = graphics:::grid, \dots)
}
\arguments{
  \item{x}{An object of class \code{"ga"}.}
  \item{y}{Not used.}
  \item{ylim}{A vector of two values specifying the limits on the y-axis.}
  \item{cex.points}{The magnification to be used for points.}
  \item{col}{The colors to be used for best and average fitness values.}
  \item{pch}{The type of points to be used for best and average fitness values.}
  \item{lty}{The type of lines to be used for best and average fitness values.}
  \item{legend}{A logical specifying if a legend should be included.}
  \item{grid}{A function for grid drawing of NULL to avoid drawing one.}
  \item{\dots}{Further arguments, currently not used.}
}

\details{Plot best and average fitness values at each iteration of GA search.}
  
\value{
The method invisibly return a \code{data.frame} with the iterations and summary statistics for the fitness function evaluated at each iteration.}

\seealso{\code{\link{ga}}, \code{\link{ga-class}}.}

%\references{}

\author{Luca Scrucca}

\examples{
# See examples in help(ga)

# The following code shows how to obtain graphs using the 
# ggplot2 plotting system 
\dontrun{
GA <- ga(...)
out <- plot(GA)
library(reshape2)
df <- melt(out[,c(1:3,5)], id.var = "iter")
library(ggplot2)
ggplot(df, aes(x = iter, y = value, group = variable, colour = variable)) +
  xlab("Generation") +  ylab("Fitness values") +
  geom_point(aes(shape = variable)) + 
  geom_line(aes(lty = variable)) +
  scale_colour_brewer(palette = "Set1") +
  theme_bw() +
  theme(legend.title = element_blank())
}
}

\keyword{methods}
\keyword{hplot}
