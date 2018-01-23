\name{ga_Mutation}
\alias{ga_Mutation}
%
\alias{gabin_raMutation}
\alias{gabin_raMutation_R}
\alias{gabin_raMutation_Rcpp}
%
\alias{gareal_raMutation}
\alias{gareal_raMutation_R}
\alias{gareal_raMutation_Rcpp}
\alias{gareal_nraMutation}
\alias{gareal_nraMutation_R}
\alias{gareal_nraMutation_Rcpp}
\alias{gareal_rsMutation}
\alias{gareal_rsMutation_R}
\alias{gareal_rsMutation_Rcpp}
\alias{gareal_powMutation}
\alias{gareal_powMutation_R}
\alias{gareal_powMutation_Rcpp}
%
\alias{gaperm_simMutation}
\alias{gaperm_simMutation_R}
\alias{gaperm_simMutation_Rcpp}
\alias{gaperm_ismMutation}
\alias{gaperm_ismMutation_R}
\alias{gaperm_ismMutation_Rcpp}
\alias{gaperm_swMutation}
\alias{gaperm_swMutation_R}
\alias{gaperm_swMutation_Rcpp}
\alias{gaperm_dmMutation}
\alias{gaperm_dmMutation_R}
\alias{gaperm_dmMutation_Rcpp}
\alias{gaperm_scrMutation}
\alias{gaperm_scrMutation_R}
\alias{gaperm_scrMutation_Rcpp}

\title{Mutation operators in genetic algorithms}

\description{Functions implementing mutation genetic operator.}

\usage{
gabin_raMutation(object, parent, \dots)

gareal_raMutation(object, parent, \dots)
gareal_nraMutation(object, parent, \dots)
gareal_rsMutation(object, parent, \dots)
gareal_powMutation(object, parent, pow = 10, \dots)

gaperm_simMutation(object, parent, \dots)
gaperm_ismMutation(object, parent, \dots)
gaperm_swMutation(object, parent, \dots)
gaperm_dmMutation(object, parent, \dots)
gaperm_scrMutation(object, parent, \dots)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An object of class \code{"ga"}, usually resulting from a call to function \code{\link{ga}}.}
  
  \item{parent}{A vector of values for the parent from the current population where mutation should occur.}
  
  \item{\dots}{Further arguments passed to or from other methods.}

  \item{pow}{Parameters of genetic operators.}
}

%\details{}
  
\value{Return a vector of values containing the mutated string.}

%\references{}

\author{Luca Scrucca}

%\note{}

%\seealso{\code{\link{ga}}}

%\examples{}

%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
