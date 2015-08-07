\name{persp3D}
\alias{persp3D}

\title{Perspective plot with colour levels}

\description{This function draws a perspective plot of a surface with different levels in different colors.}

\usage{
persp3D(x, y, z, theta = 30, phi = 20, d = 5, expand = 2/3,
        xlim = range(x, finite = TRUE), ylim = range(y, finite = TRUE), 
        zlim = range(z, finite = TRUE), levels = pretty(zlim, nlevels),
        nlevels = 20, color.palette = jet.colors, border = NA, 
        ticktype = "detailed", xlab = NULL, ylab = NULL, zlab = NULL, 
        \dots)
}

\arguments{
  \item{x, y}{locations of grid lines at which the values in \code{z} are measured. These must be in ascending order. By default, equally spaced values from 0 to 1 are used. If \code{x} is a list, its components \code{x$x} and \code{x$y} are used for \code{x} and \code{y}, respectively.}

  \item{z}{a matrix containing the values to be plotted (NAs are allowed).} 

  \item{theta, phi}{angles defining the viewing direction. \code{theta} gives the azimuthal direction and \code{phi} the colatitude.}

  \item{d}{a value which can be used to vary the strength of the perspective transformation.}
  
  \item{expand}{a expansion factor applied to the \code{z} coordinates.}

  \item{xlim, ylim, zlim}{x-, y- and z-limits for the axes.}
  
  \item{levels}{a vector of values specifying the levels to be used for plotting the surface with different colors.}

  \item{nlevels}{a value specifying the numbe of levels to be used for plotting. This value is used if \code{levels} argument is not specified.}
  
  \item{color.palette}{the color palette used for plotting.}
  
  \item{border}{the color of the line drawn around the surface facets. By default is set to \code{NA} so no borders are drawn.}
  
  \item{ticktype}{a character specifying the type of axes tickmarks. By default \code{"detailed"} ticks are drawn.}
  
  \item{xlab, ylab, zlab}{character strings specifying the titles for the axes.}

  \item{\dots}{Further arguments passed to the function \code{\link{persp}}.}
}

\details{This function enhances the default perspective plot for drawing 3-dimensional surfaces.}

\value{Return a list with the following elements:
\item{persp}{the viewing transformation matrix (see \code{link{persp}});}
\item{levels}{a vector of values giving the levels used for plotting the surface;}
\item{colors}{a vector of strings giving the color used for plotting the surface.}
}

%\references{}

\author{Luca Scrucca}

%\note{}

\seealso{\code{link{persp}}}

\examples{
y <- x <- seq(-10, 10, length=60)
f <- function(x,y) { r <- sqrt(x^2+y^2); 10 * sin(r)/r }
z <- outer(x, y, f)
persp3D(x, y, z, theta = 30, phi = 30, expand = 0.5)
persp3D(x, y, z, color.palette = heat.colors, phi = 30, theta = 225, 
        box = TRUE, border = NA, shade = .4)

x1 = seq(-3,3,length=50)
x2 = seq(-3,3,length=50)
y = function(x1, x2) sin(x1)+cos(x2)
persp3D(x1, x2, outer(x1,x2,y), zlab="y", theta = 150, phi = 20, expand = 0.6)
}

\keyword{hplot}

