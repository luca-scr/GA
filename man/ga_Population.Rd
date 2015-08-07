\name{ga_Population}
\alias{ga_Population}
%
\alias{gabin_Population}
\alias{gareal_Population}
\alias{gaperm_Population}

\title{Population initialization in genetic algorithms}

\description{Functions for creating a random initial population to be used in genetic algorithms.}

\usage{
gabin_Population(object, ...)

gareal_Population(object, ...)

gaperm_Population(object, ...) 
}

\arguments{

  \item{object}{An object of class \code{"ga"}, usually resulting from a call to function \code{\link{ga}}.}
  
  \item{\dots}{Further arguments passed to or from other methods.}
}

\details{
\code{gabin_Population} generates a random population of \code{object@nBits} binary values; 

\code{gareal_Population} generates a random (uniform) population of real values in the range [\code{object@min}, \code{object@max}]; 

\code{gaperm_Population} generates a random (uniform) population of integer values in the range [\code{object@min}, \code{object@max}].
}

\value{Return a matrix of dimension \code{object@popSize} times the number of decision variables.}

%\references{}

\author{Luca Scrucca}

%\note{}

\seealso{\code{\link{ga}}}

%\examples{}

%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
