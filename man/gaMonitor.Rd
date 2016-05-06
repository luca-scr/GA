\name{gaMonitor}
\alias{gaMonitor}
\alias{gaMonitor2}
\alias{gaislMonitor}
\alias{gaislMonitor2}

\title{Monitor genetic algorithm evolution}

\description{Functions to print summary statistics of fitness values at each iteration of a GA search.}

\usage{
gaMonitor(object, digits = getOption("digits"), \dots)

gaMonitor2(object, digits = getOption("digits"), \dots)

gaislMonitor(object, digits = getOption("digits"), \dots)

gaislMonitor2(object, digits = getOption("digits"), \dots)
}

\arguments{
  \item{object}{an object of class \code{ga-class} or \code{gaisl-class}, usually resulting from a call to function \code{\link{ga}} or \code{\link{gaisl}}, respectively.}
  \item{digits}{minimal number of significant digits.}
  \item{\dots}{further arguments passed to or from other methods.}
}

%\details{}
\value{These functions print a summary of current GA step on the console.\cr
By default, \code{gaMonitor} is called in interactive sessions by \code{\link{ga}}. The old monitoring function, used as the default until version 2.2 of \pkg{GA} package, is provided in \code{gaMonitor2}.\cr
By default, \code{gaislMonitor} is called in interactive sessions by \code{\link{gaisl}}.
}

%\references{}

\author{Luca Scrucca}

%\seealso{\code{\link{ga}}}

%\note{}

%\examples{}

%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
