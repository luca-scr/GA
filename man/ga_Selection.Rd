\name{ga_Selection}
\alias{ga_Selection}
%
\alias{ga_lrSelection}
\alias{ga_nlrSelection}
\alias{ga_rwSelection}
\alias{ga_tourSelection}
%
\alias{gabin_lrSelection}
\alias{gabin_nlrSelection}
\alias{gabin_rwSelection}
\alias{gabin_tourSelection}
%
\alias{gareal_lrSelection}
\alias{gareal_nlrSelection}
\alias{gareal_rwSelection}
\alias{gareal_tourSelection}
\alias{gareal_lsSelection}
\alias{gareal_sigmaSelection}
%
\alias{gaperm_lrSelection}
\alias{gaperm_nlrSelection}
\alias{gaperm_rwSelection}
\alias{gaperm_tourSelection}

\title{Selection operators in genetic algorithms}

\description{Functions implementing selection genetic operator.}

\usage{
ga_lrSelection(object, r = 2/(object@popSize * (object@popSize - 1)),
                       q = 2/object@popSize, ...)
ga_nlrSelection(object, q = 0.25, ...)
ga_rwSelection(object, ...)
ga_tourSelection(object, k = 3, ...)

gabin_lrSelection(object, r = 2/(object@popSize * (object@popSize - 1)), 
                          q = 2/object@popSize, ...)
gabin_nlrSelection(object, q = 0.25, ...)
gabin_rwSelection(object, ...)
gabin_tourSelection(object, k = 3, ...)

gareal_lrSelection(object, r = 2/(object@popSize * (object@popSize - 1)), 
                           q = 2/object@popSize, ...)
gareal_nlrSelection(object, q = 0.25, ...)
gareal_rwSelection(object, ...)
gareal_tourSelection(object, k = 3, ...)
gareal_lsSelection(object, ...)
gareal_sigmaSelection(object, ...)

gaperm_lrSelection(object, r = 2/(object@popSize * (object@popSize - 1)), 
                           q = 2/object@popSize, ...)
gaperm_nlrSelection(object, q = 0.25, ...)
gaperm_rwSelection(object, ...)
gaperm_tourSelection(object, k = 3, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An object of class \code{"ga"}, usually resulting from a call to function \code{\link{ga}}.}

  \item{r}{A tuning parameter for the specific selection operator.}
    
  \item{q}{A tuning parameter for the specific selection operator.}
    
  \item{k}{A tuning parameter for the specific selection operator.}

  \item{\dots}{Further arguments passed to or from other methods.}

}

%\details{}

\value{
Return a list with two elements:

\item{population}{a matrix of dimension \code{object@popSize} times the number of decision variables containing the selected individuals or strings;}

\item{fitness}{a vector of length \code{object@popSize} containing the fitness values for the selected individuals.}
}

%\references{}
\author{Luca Scrucca}
%\note{}

\seealso{\code{\link{ga}}}

%\examples{}

%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
