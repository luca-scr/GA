\name{binary2gray}
\alias{binary2gray}
\alias{gray2binary}

\title{Gray encoding for binary strings}

\description{Functions for computing Gray encoding from/to binary strings.}

\usage{
binary2gray(x)
gray2binary(x)
}

\arguments{
  \item{x}{the string to be evaluated}
}

\details{
Gray encoding allows to obtain binary strings not affected by the well-known Hamming cliff problem. With Gray encoding the number of bit differences between any two consecutive values is one, whereas in binary strings this is not always true.
}

%\value{}
%\references{}

\author{Luca Scrucca}

\seealso{\code{\link{binary2decimal}}}

\examples{
# Consider a five-bit encoding of values 15 and 16  using the standard 
# binary coding
decimal2binary(15, 5)
decimal2binary(16, 5)
# Moving from 15 to 16 (or viceversa) all five bits need to be changed,
# but using Gray encoding the two binary strings differ by one bit.
binary2gray(decimal2binary(15, 5))
binary2gray(decimal2binary(16, 5))
}
