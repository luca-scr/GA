\name{ga_Crossover}
\alias{ga_Crossover}
\alias{ga_Crossover_R}
\alias{ga_Crossover_Rcpp}
%
\alias{ga_spCrossover}
\alias{ga_spCrossover_R}
\alias{ga_spCrossover_Rcpp}
%
\alias{gabin_spCrossover}
\alias{gabin_spCrossover_R}
\alias{gabin_spCrossover_Rcpp}
\alias{gabin_uCrossover}
\alias{gabin_uCrossover_R}
\alias{gabin_uCrossover_Rcpp}
%
\alias{gareal_spCrossover}
\alias{gareal_spCrossover_R}
\alias{gareal_spCrossover_Rcpp}
\alias{gareal_waCrossover}
\alias{gareal_waCrossover_R}
\alias{gareal_waCrossover_Rcpp}
\alias{gareal_laCrossover}
\alias{gareal_laCrossover_R}
\alias{gareal_laCrossover_Rcpp}
\alias{gareal_blxCrossover}
\alias{gareal_blxCrossover_R}
\alias{gareal_blxCrossover_Rcpp}
\alias{gareal_laplaceCrossover}
\alias{gareal_laplaceCrossover_R}
\alias{gareal_laplaceCrossover_Rcpp}
%
\alias{gaperm_cxCrossover}
\alias{gaperm_cxCrossover_R}
\alias{gaperm_cxCrossover_Rcpp}
\alias{gaperm_pmxCrossover}
\alias{gaperm_pmxCrossover_R}
\alias{gaperm_pmxCrossover_Rcpp}
\alias{gaperm_oxCrossover}
\alias{gaperm_oxCrossover_R}
\alias{gaperm_oxCrossover_Rcpp}
\alias{gaperm_pbxCrossover}
\alias{gaperm_pbxCrossover_R}
\alias{gaperm_pbxCrossover_Rcpp}

\title{Crossover operators in genetic algorithms}

\description{Functions implementing crossover genetic operator.}

\usage{
ga_spCrossover(object, parents, \dots)

gabin_spCrossover(object, parents, \dots)
gabin_uCrossover(object, parents, \dots)

gareal_spCrossover(object, parents, \dots)
gareal_waCrossover(object, parents, \dots)
gareal_laCrossover(object, parents, \dots)
gareal_blxCrossover(object, parents, a = 0.5, \dots)
gareal_laplaceCrossover(object, parents, a = 0, b = 0.15, \dots)

gaperm_cxCrossover(object, parents, \dots)
gaperm_pmxCrossover(object, parents, \dots)
gaperm_oxCrossover(object, parents, \dots)
gaperm_pbxCrossover(object, parents, \dots)

}

\arguments{
  \item{object}{An object of class \code{"ga"}, usually resulting from a call to function \code{\link{ga}}.}
  
  \item{parents}{A two-rows matrix of values indexing the parents from the current population.}
  
  \item{\dots}{Further arguments passed to or from other methods.}
  
  \item{a, b}{Parameters of genetic operators.}
}

%\details{}

\value{
Return a list with two elements:

\item{children}{a matrix of dimension 2 times the number of decision variables containing the generated offsprings;}
\item{fitness}{a vector of length 2 containing the fitness values for the offsprings. A value \code{NA} is returned if an offspring is different (which is usually the case) from the two parents.}
}

%\references{}
\author{Luca Scrucca}
%\note{}

\seealso{\code{\link{ga}}}

%\examples{}

%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
