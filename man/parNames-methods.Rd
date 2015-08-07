\name{parNames-methods}
\docType{methods}
\alias{parNames}
\alias{parNames,ga-method}

\title{Parameters or decision variables names from an object of class \code{\link{ga-class}}.}

\description{
  A method for obtaining the names of parameters or decision variables 
  from an object of class \code{\link{ga-class}}.
}

\usage{
parNames(object, \dots)
\S4method{parNames}{ga}(object, \dots)
}

\arguments{
  \item{object}{An object of class \code{"ga"}, usually resulting from a call to function \code{\link{ga}}.}
  \item{\dots}{Further arguments, currently not used.}
}

%\details{}
  
\value{A list of character values providing the names of parameters or decision variables.}

%\references{}

\author{Luca Scrucca}

\seealso{\code{\link{ga}}}

%\examples{}
\keyword{methods}
