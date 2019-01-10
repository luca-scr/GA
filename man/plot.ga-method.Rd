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
  \item{col}{The colours to be used for best and average fitness values.}
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
GA <- ga(type = "real-valued", 
         fitness = function(x) -(abs(x)+cos(x)), 
         lower = -20, upper = 20, 
         popSize = 20, pmutation = 0.2, maxiter = 50)
out <- plot(GA)
library(reshape2)
df <- melt(out[,c(1:3,5)], id.var = "iter")
library(ggplot2)
ggplot(out) +
  geom_ribbon(aes(x = iter, ymin = median, ymax = max, 
                  colour = "median", fill = "median")) +
  geom_line(aes(x = iter, y = max, colour = "max")) +
  geom_point(aes(x = iter, y = max, colour = "max")) +
  geom_line(aes(x = iter, y = mean, colour = "mean"), lty = 2) +
  geom_point(aes(x = iter, y = mean, colour = "mean"), pch = 1) +
  xlab("Generation") +  ylab("Fitness values") +
  scale_colour_manual(breaks = c("max", "mean", "median"),
                      values = c("green3", "dodgerblue3", adjustcolor("green3", alpha.f = 0.1))) +
  scale_fill_manual(breaks = "median",
                    values = adjustcolor("green3", alpha.f = 0.1)) +
  guides(fill = "none", 
         colour = guide_legend(override.aes = 
                                 list(fill = c(NA, NA, adjustcolor("green3", alpha.f = 0.1)),
                                       pch = c(19,1,NA)))) +
  theme_bw() +
  theme(legend.title = element_blank(), 
        legend.pos = "top",
        legend.background = element_blank())
}
}

\keyword{methods}
\keyword{hplot}
