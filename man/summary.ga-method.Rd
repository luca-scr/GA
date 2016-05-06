\name{summary.ga-method}
\docType{methods}
%\alias{summary.ga-class}
\alias{summary,ga-method}
\alias{summary.ga}
\alias{print.summary.ga}

\title{Summary for Genetic Algorithms}

\description{Summary method for class \code{\link{ga-class}}.}

\usage{
\S4method{summary}{ga}(object, \dots)

\S3method{print}{summary.ga}(x, digits = getOption("digits"), \dots)
}

\arguments{
\item{object}{an object of class \code{\link{ga-class}}.}

\item{x}{an object of class \code{summary.ga}.}

\item{digits}{number of significant digits.}

\item{\dots}{further arguments passed to or from other methods.}
}

%\details{}
  
\value{The \code{summary} function returns an object of class \code{summary.ga} which can be printed by the corresponding \code{print} method. The function also returns invisibly a list with the information from the genetic algorithm search.}

%\references{}
  
\author{Luca Scrucca}

%\note{}

\seealso{\code{\link{ga}}}

\examples{
f <- function(x)  abs(x)+cos(x)
GA <- ga(type = "real-valued", 
         fitness = function(x) -f(x), 
         min = -20, max = 20, run = 50)
out <- summary(GA)
print(out)
str(out)
}

\keyword{optimize}
