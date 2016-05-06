\name{gaSummary}
\alias{gaSummary}

\title{Summarize genetic algorithm evolution}

\description{A function which returns fitness summary statistics at each iteration of GA search.}

\usage{
gaSummary(x, \dots)
}

\arguments{
  \item{x}{a vector of fitness values for which summary statistics should be computed.}
  \item{\dots}{further arguments passed to or from other methods.}
}

\details{This function computes summary statistics for a vector of fitness values at current iteration of GA search.}

\value{A vector with the following values: \code{(max, mean, median, min)}
}

%\references{}

\author{Luca Scrucca}

\seealso{\code{\link{ga}}}

%\note{}

%\examples{}

%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
