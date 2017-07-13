\name{ga_Crossover}
\alias{ga_Crossover}
%
\alias{ga_spCrossover}
%
\alias{gabin_spCrossover}
\alias{gabin_uCrossover}
%
\alias{gareal_spCrossover}
\alias{gareal_waCrossover}
\alias{gareal_laCrossover}
\alias{gareal_blxCrossover}
\alias{gareal_laplaceCrossover}
%
\alias{gaperm_cxCrossover}
\alias{gaperm_pmxCrossover}
\alias{gaperm_oxCrossover}
\alias{gaperm_pbxCrossover}

\title{Crossover operators in genetic algorithms}

\description{Functions implementing crossover genetic operator.}

\usage{
ga_spCrossover(object, parents, \dots)

gabin_spCrossover(object, parents, \dots)
gabin_uCrossover(object, parents, \dots)

gareal_spCrossover(object, parents, \dots)
gareal_waCrossover(object, parents, \dots)
gareal_laCrossover(object, parents, \dots)
gareal_blxCrossover(object, parents, \dots)
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
