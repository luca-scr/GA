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
%
\alias{gaperm_cxCrossover}
\alias{gaperm_pmxCrossover}
\alias{gaperm_oxCrossover}
\alias{gaperm_pbxCrossover}

\title{Crossover operators in genetic algorithms}

\description{Functions implementing crossover genetic operator.}

\usage{
ga_spCrossover(object, parents, ...)

gabin_spCrossover(object, parents, ...)
gabin_uCrossover(object, parents, ...)

gareal_spCrossover(object, parents, ...)
gareal_waCrossover(object, parents, ...)
gareal_laCrossover(object, parents, ...)
gareal_blxCrossover(object, parents, ...)

gaperm_cxCrossover(object, parents, ...)
gaperm_pmxCrossover(object, parents, ...)
gaperm_oxCrossover(object, parents, ...)
gaperm_pbxCrossover(object, parents, ...)

}

\arguments{
  \item{object}{An object of class \code{"ga"}, usually resulting from a call to function \code{\link{ga}}.}
  
  \item{parents}{A two-rows matrix of values indexing the parents from the current population.}
  
  \item{\dots}{Further arguments passed to or from other methods.}
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
