\name{summary.de-method}
\docType{methods}
%\alias{summary.de-class}
\alias{summary,de-method}
\alias{summary.de}
\alias{print.summary.de}

\title{Summary for Differential Evolution}

\description{Summary method for class \code{\link{de-class}}.}

\usage{
\S4method{summary}{de}(object, \dots)

\S3method{print}{summary.de}(x, digits = getOption("digits"), \dots)
}

\arguments{
\item{object}{an object of class \code{\link{de-class}}.}

\item{x}{an object of class \code{summary.de}.}

\item{digits}{number of significant digits.}

\item{\dots}{further arguments passed to or from other methods.}
}

%\details{}
  
\value{The \code{summary} function returns an object of class \code{summary.de} which can be printed by the corresponding \code{print} method. The function also returns invisibly a list with the information from the differential evolution search.}

%\references{}
  
\author{Luca Scrucca}

%\note{}

\seealso{\code{\link{de}}}

\examples{
f <- function(x) abs(x)+cos(x)
DE <- de(fitness = function(x) -f(x), 
         lower = -20, upper = 20, run = 50)
out <- summary(DE)
print(out)
str(out)
}

\keyword{optimize}
