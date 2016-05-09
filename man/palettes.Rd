\name{palettes}
\alias{jet.colors}
\alias{spectral.colors}
\alias{bl2gr.colors}

\title{Colors palettes}

\description{Functions for creating a vector of colors from pre-speficied palettes.}

\usage{
jet.colors(n)

spectral.colors(n)

bl2gr.colors(n)
}

\arguments{
  \item{n}{a numerical value specifying the number of colors in the palette.}
}

\details{\code{jet.colors()} creates a palette of colors which tend to have high brightness and not uniform luminance. Furthermore, the brightest colors, yellow and cyan, are used for intermediate data values, and this has the effect of emphasizing uninteresting (and arbitrary) values while demphasizing the extremes. For these reasons this popular palette is not recommended.\cr
\code{spectral.colors()} creates a palette based on ColorBrewer \url{http://colorbrewer2.org}, so the resulting colors have a much uniform luminance.\cr
The \code{bl2gr.colors()} palette returns a palette of colors from blue to green.
}

\value{Returns a character vector of colors encoded in hexadecimal values.} 

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
