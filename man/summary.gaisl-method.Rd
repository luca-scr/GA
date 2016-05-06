\name{summary.gaisl-method}
\docType{methods}
\alias{summary,gaisl-method}
\alias{summary.gaisl}
\alias{print.summary.gaisl}

\title{Summary for Islands Genetic Algorithms}

\description{Summary method for class \code{\link{gaisl-class}}.}

\usage{
\S4method{summary}{gaisl}(object, \dots)

\S3method{print}{summary.gaisl}(x, digits = getOption("digits"), \dots)
}

\arguments{
\item{object}{an object of class \code{\link{gaisl-class}}.}

\item{x}{an object of class \code{summary.gaisl}.}

\item{digits}{number of significant digits.}

\item{\dots}{further arguments passed to or from other methods.}
}

%\details{}
  
\value{The \code{summary} function returns an object of class \code{summary.gaisl} which can be printed by the corresponding \code{print} method. The function also returns invisibly a list with the information from the islands genetic algorithm search.}

%\references{}
  
\author{Luca Scrucca}

%\note{}

\seealso{\code{\link{gaisl}}}

\examples{
f <- function(x)  abs(x)+cos(x)
GA <- gaisl(type = "real-valued", 
            fitness = function(x) -f(x), 
            min = -20, max = 20, run = 10,
            numIslands = 4, parallel = 2) # CRAN check limits ncores = 2
out <- summary(GA)
print(out)
str(out)
}

\keyword{optimize}
