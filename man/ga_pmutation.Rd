\name{ga_pmutation}
\alias{ga_pmutation}

\title{Variable mutation probability in genetic algorithms}

\description{A function which calculates the mutation probability for the current iteration. This enables to use GAs with variable mutation rate (see examples).}

\usage{
ga_pmutation(object, p0 = 0.5, p = 0.01, T = round(object@maxiter/2), \dots)
}

\arguments{
  \item{object}{An object of class \code{"ga"}, usually resulting from a call to function \code{\link{ga}}.}
  
  \item{p0}{initial probability of mutation.}
  \item{p}{limiting probability of mutation.}
  \item{T}{maximum iteration after which it should converges to \code{p}.}  
  \item{\dots}{Further arguments passed to or from other methods.}
}

%\details{}
  
\value{Return a numeric value in the range (0,1).}

%\references{}

\author{Luca Scrucca}

%\note{}

\seealso{\code{\link{ga}}, \code{\link{ga_Mutation}}}

\examples{
\dontrun{
Rastrigin <- function(x1, x2)
{
  20 + x1^2 + x2^2 - 10*(cos(2*pi*x1) + cos(2*pi*x2))
}

GA <- ga(type = "real-valued", 
         fitness =  function(x) -Rastrigin(x[1], x[2]),
         min = c(-5.12, -5.12), max = c(5.12, 5.12), 
         popSize = 50, maxiter = 500, run = 100,
         pmutation = ga_pmutation)
plot(GA)

GA <- ga(type = "real-valued", 
         fitness =  function(x) -Rastrigin(x[1], x[2]),
         min = c(-5.12, -5.12), max = c(5.12, 5.12), 
         popSize = 50, maxiter = 500, run = 100,
         pmutation = function(...) ga_pmutation(..., p0 = 0.1))
plot(GA)
}
}

%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
