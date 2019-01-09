\name{palettes}
\alias{jet.colors}
\alias{spectral.colors}
\alias{bl2gr.colors}

\title{Colours palettes}

\description{Functions for creating a vector of colours from pre-specified palettes.}

\usage{
jet.colors(n)

spectral.colors(n)

bl2gr.colors(n)
}

\arguments{
  \item{n}{a numerical value specifying the number of colours in the palette.}
}

\details{\code{jet.colors()} creates a palette of colours which tend to have high brightness and not uniform luminance. Furthermore, the brightest colours, yellow and cyan, are used for intermediate data values, and this has the effect of emphasizing uninteresting (and arbitrary) values while de-emphasizing the extremes. For these reasons this popular palette is not recommended.\cr
\code{spectral.colors()} creates a palette based on ColorBrewer \url{http://colorbrewer2.org}, so the resulting colours have a much uniform luminance.\cr
The \code{bl2gr.colors()} palette returns a palette of colours from blue to green.
}

\value{Returns a character vector of colours encoded in hexadecimal values.} 

%\references{}
%\author{}
%\note{}

\seealso{\code{\link{colors}}, \code{\link{rgb}}.}

\examples{
jet.colors(9)
spectral.colors(9)
bl2gr.colors(9)

par(mfrow = c(3,1), mar = c(1,1,1,1))
n = 21
image(1:21, 1, as.matrix(1:21), col = jet.colors(21),
      ylab = "", xlab = "", xaxt = "n", yaxt = "n", bty = "n")
image(1:21, 1, as.matrix(1:21), col = spectral.colors(21),
      ylab = "", xlab = "", xaxt = "n", yaxt = "n", bty = "n")
image(1:21, 1, as.matrix(1:21), col = bl2gr.colors(21),
      ylab = "", xlab = "", xaxt = "n", yaxt = "n", bty = "n")
}

\keyword{hplot}
