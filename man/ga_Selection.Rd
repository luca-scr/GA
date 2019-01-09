\name{ga_Selection}
\alias{ga_Selection}
%
\alias{ga_lrSelection}
\alias{ga_lrSelection_R}
\alias{ga_lrSelection_Rcpp}
\alias{ga_nlrSelection}
\alias{ga_nlrSelection_R}
\alias{ga_nlrSelection_Rcpp}
\alias{ga_rwSelection}
\alias{ga_rwSelection_R}
\alias{ga_rwSelection_Rcpp}
\alias{ga_tourSelection}
\alias{ga_tourSelection_R}
\alias{ga_tourSelection_Rcpp}
%
\alias{gabin_lrSelection}
\alias{gabin_lrSelection_R}
\alias{gabin_lrSelection_Rcpp}
\alias{gabin_nlrSelection}
\alias{gabin_nlrSelection_R}
\alias{gabin_nlrSelection_Rcpp}
\alias{gabin_rwSelection}
\alias{gabin_rwSelection_R}
\alias{gabin_rwSelection_Rcpp}
\alias{gabin_tourSelection}
\alias{gabin_tourSelection_R}
\alias{gabin_tourSelection_Rcpp}
%
\alias{gareal_lrSelection}
\alias{gareal_lrSelection_R}
\alias{gareal_lrSelection_Rcpp}
\alias{gareal_nlrSelection}
\alias{gareal_nlrSelection_R}
\alias{gareal_nlrSelection_Rcpp}
\alias{gareal_rwSelection}
\alias{gareal_rwSelection_R}
\alias{gareal_rwSelection_Rcpp}
\alias{gareal_tourSelection}
\alias{gareal_tourSelection_R}
\alias{gareal_tourSelection_Rcpp}
\alias{gareal_lsSelection}
\alias{gareal_lsSelection_R}
\alias{gareal_lsSelection_Rcpp}
\alias{gareal_sigmaSelection}
\alias{gareal_sigmaSelection_R}
\alias{gareal_sigmaSelection_Rcpp}
%
\alias{gaperm_lrSelection}
\alias{gaperm_lrSelection_R}
\alias{gaperm_lrSelection_Rcpp}
\alias{gaperm_nlrSelection}
\alias{gaperm_nlrSelection_R}
\alias{gaperm_nlrSelection_Rcpp}
\alias{gaperm_rwSelection}
\alias{gaperm_rwSelection_R}
\alias{gaperm_rwSelection_Rcpp}
\alias{gaperm_tourSelection}
\alias{gaperm_tourSelection_R}
\alias{gaperm_tourSelection_Rcpp}
%
\alias{gareal_de}
\alias{gareal_de_R}
\alias{gareal_de_Rcpp}

\title{Selection operators in genetic algorithms}

\description{Functions implementing selection genetic operator.}

\usage{
ga_lrSelection(object, r = 2/(object@popSize * (object@popSize - 1)),
                       q = 2/object@popSize, \dots)
ga_nlrSelection(object, q = 0.25, \dots)
ga_rwSelection(object, \dots)
ga_tourSelection(object, k = 3, \dots)

gabin_lrSelection(object, r = 2/(object@popSize * (object@popSize - 1)), 
                          q = 2/object@popSize, \dots)
gabin_nlrSelection(object, q = 0.25, \dots)
gabin_rwSelection(object, \dots)
gabin_tourSelection(object, k = 3, \dots)

gareal_lrSelection(object, r = 2/(object@popSize * (object@popSize - 1)), 
                           q = 2/object@popSize, \dots)
gareal_nlrSelection(object, q = 0.25, \dots)
gareal_rwSelection(object, \dots)
gareal_tourSelection(object, k = 3, \dots)
gareal_lsSelection(object, \dots)
gareal_sigmaSelection(object, \dots)

gaperm_lrSelection(object, r = 2/(object@popSize * (object@popSize - 1)), 
                           q = 2/object@popSize, \dots)
gaperm_nlrSelection(object, q = 0.25, \dots)
gaperm_rwSelection(object, \dots)
gaperm_tourSelection(object, k = 3, \dots)

gareal_de(object, F = 0.8, p = 0.5, \dots) 
}

\arguments{
  \item{object}{An object of class \code{"ga"}, usually resulting from a call to function \code{\link{ga}}.}

  \item{r}{A tuning parameter for the GA selection operator.}
    
  \item{q}{A tuning parameter for the GA selection operator.}
    
  \item{k}{A tuning parameter for the GA selection operator.}
  
  \item{F, p}{Tuning parameters for the DE selection operator.}

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

\seealso{\code{\link{ga}}, \code{\link{de}}.}

%\examples{}
