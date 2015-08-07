\name{ga_Mutation}
\alias{ga_Mutation}
%
\alias{gabin_raMutation}
%
\alias{gareal_raMutation}
\alias{gareal_nraMutation}
\alias{gareal_rsMutation}
%
\alias{gaperm_simMutation}
\alias{gaperm_ismMutation}
\alias{gaperm_swMutation}
\alias{gaperm_dmMutation}
\alias{gaperm_scrMutation}

\title{Mutation operators in genetic algorithms}

\description{Functions implementing mutation genetic operator.}

\usage{
gabin_raMutation(object, parent, ...)

gareal_raMutation(object, parent, ...)
gareal_nraMutation(object, parent, ...)
gareal_rsMutation(object, parent, ...)

gaperm_simMutation(object, parent, ...)
gaperm_ismMutation(object, parent, ...)
gaperm_swMutation(object, parent, ...)
gaperm_dmMutation(object, parent, ...)
gaperm_scrMutation(object, parent, ...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An object of class \code{"ga"}, usually resulting from a call to function \code{\link{ga}}.}
  
  \item{parent}{A vector of values for the parent from the current population where mutation should occur.}
  
  \item{\dots}{Further arguments passed to or from other methods.}
}

%\details{}
  
\value{Return a vector of values containing the mutated string.}

%\references{}

\author{Luca Scrucca}

%\note{}

\seealso{\code{\link{ga}}}

%\examples{}

%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
