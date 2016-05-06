\name{palettes}
\alias{jet.colors}
\alias{spectral.colors}
\alias{bl2gr.colors}

\title{Diverging colors palettes}

\description{Functions for creating a vector of diverging colors.}

\usage{
jet.colors(n)

spectral.colors(n)

bl2gr.colors(n)
}

\arguments{
  \item{n}{a numerical value specifying the number of colors in the palette.}
}

\details{Jet colors tend to have high brightness and not uniform luminance. Furthermore, the brightest colors, yellow and cyan, are used for intermediate data values, and this has the effect of emphasizing uninteresting (and arbitrary) values while demphasizing the extremes. For these reasons this popular palette is not recommended.\cr
The spectral palette is based on ColorBrewer \url{http://colorbrewer2.org} and the resulting colors have a much uniform luminance.\cr
The `bl2gr.colors` palette returns colors from blue to green. 
}

\value{Returns vector of \code{n} color names.} 

%\references{}
%\author{}
%\note{}

\seealso{\code{\link{colors}}.}

\examples{
jet.colors(9)
spectral.colors(9)

par(mfrow = c(2,1), mar = c(1,1,1,1))
n = 21
image(1:21, 1, as.matrix(1:21), col = jet.colors(21),
      ylab = "", xlab = "", xaxt = "n", yaxt = "n", bty = "n")
image(1:21, 1, as.matrix(1:21), col = spectral.colors(21),
      ylab = "", xlab = "", xaxt = "n", yaxt = "n", bty = "n")
}

\keyword{hplot}
