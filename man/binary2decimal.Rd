\name{binary2decimal}
\alias{binary2decimal}
\alias{decimal2binary}

\title{Binary encoding of decimal numbers and viceversa.}

\description{Functions for computing binary to decimal conversion of numbers and viceversa.}

\usage{
decimal2binary(x, length)
binary2decimal(x)
}

\arguments{
  \item{x}{input value.}
  \item{length}{an optional value giving the length of binary string to return.}
}

\details{
\code{decimal2binary} converts a numerical value (which is forced to be an integer) to a binary representation, i.e. a vector of 0s and 1s. For real numerical values see the example below.

\code{binary2binary} converts a binary value, i.e. a vector of 0s and 1s, to a decimal representation.
}

%\value{}
%\references{}
\author{Luca Scrucca}

\seealso{\code{\link{binary2gray}}}

\examples{
# for integer values
dval <- 12
(bval <- decimal2binary(dval))
binary2decimal(bval)

# for real values
dval <- 12.456
# use
(bval <- decimal2binary(dval*1000))
binary2decimal(bval)/1000
}

%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
