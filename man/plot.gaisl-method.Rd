\name{plot.gaisl-method}
\docType{methods}
\alias{plot,gaisl-method}
\alias{plot.gaisl}

\title{Plot of Islands Genetic Algorithm search path}

\description{
  The \code{plot} method for \code{\link{gaisl-class}} objects gives a plot
  of best fitness values found in each island during the GA iterations.
}

\usage{
\S4method{plot}{gaisl}(x, y, \dots)
}
\arguments{
  \item{x}{An object of class \code{"gaisl"}.}
  \item{y}{Not used.}
  \item{\dots}{Further arguments passed to \code{\link{plot.default}}, such as \code{ylim}, \code{ylab}, etc., or to \code{\link{matplot}}, such as \code{col}, \code{lty}, and \code{lwd}.}
}

\details{Plot best fitness values found in each island during the GA iterations.}
  
\value{
The method invisibly return a list with the following components:

\item{iter}{a vector of values specifying the interation.}
\item{summary}{a matrix of best fitness values for each island along the columns.}
}

\seealso{\code{\link{gaisl}}, \code{\link{gaisl-class}}.}

%\references{}

\author{Luca Scrucca}

\examples{
# See examples in help(gaisl)

# The following code shows how to obtain graphs using the 
# ggplot2 plotting system 
\dontrun{
GA <- gaisl(...)
out <- plot(GA)
library(reshape2)
df <- melt(out, id.var = "iter")
library(ggplot2)
ggplot(df, aes(x = iter, y = value, group = variable, colour = variable)) +
  xlab("Generation") +  ylab("Fitness values") +
  geom_line(size = 0.8) +
  scale_colour_brewer(palette = "Paired") +
  theme_bw() +
  theme(legend.title = element_blank())
}
}

\keyword{methods}
\keyword{hplot}
