\name{jet.colors}
\alias{jet.colors}

\title{Jet Colors Palette}

\description{Create a vector of n colors beginning with dark blue, ranging through shades of blue, cyan, green, yellow and red, and ending with dark red.
}

\usage{
jet.colors(n)
}

\arguments{
  \item{n}{a numerical value specifying the number of colors in the palette.}
}

\details{This function creates a palette of colors beginning with dark blue, ranging through shades of blue, cyan, green, yellow and red, and ending with dark red.}

\value{Returns vector of \code{n} color names.} 

%\references{}
%\author{}
%\note{}

\seealso{\code{\link{colors}}.}

\examples{
jet.colors(5)

palette(jet.colors(21))
pie(rep(1,21), col = 1:21)

}

\keyword{hplot}
