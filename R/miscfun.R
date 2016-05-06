##
## Conversion from/to decimal-binary
##

decimal2binary <- function(x, length)
{
  x <- as.integer(x)
  b <- if(missing(length)) NULL else rep(0, length)
  i <- 0
  while(x >= 1)
       { i <- i + 1
         b[i] <- x %% 2
         x <- x %/% 2 }
  return(rev(b))
}

binary2decimal <- function(x)
{
  sum(x * 2^(rev(seq(along=x)) - 1))  
}

## old versions
# decimal2binary <- function(x)
# {
#   x <- as.integer(x)
#   b <- NULL
#   while(x > 1)
#        { b <- c(x %% 2, b)
#          x <- x %/% 2 }
#   b <- c(x %% 2, b)
#   return(b)
# }
#
# binary2decimal <- function(x)
# {
#   l <- length(x)
#   i <- seq(l-1, 0)
#   d <- sum(x*2^i)
#   return(d)
# }



##
## Gray coding for binary genetic algorithm
## Based on algorithm on Eiben Smith (2003) Introduction to Evolutionary Computing
##

binary2gray <- function(x)
{
  x <- as.logical(x)
  n <- length(x)
  g <- vector(mode = "logical", length = n)
  g[1] <- x[1]
  if(n > 1)
    for(i in 2:n)
       { g[i] <- xor(x[i-1], x[i]) }
  g <- as.numeric(g)
  return(g)
}

gray2binary <- function(x)
{
  x <- as.logical(x)
  n <- length(x)
  b <- vector(mode = "logical", length = n)
  b[1] <- value <- x[1]
  if(n > 1)
    for(i in 2:n)
     { if(x[i]) value <- !value
       b[i] <- value }
  b <- as.numeric(b)
  return(b)
}

#############################################################################

clearConsoleLine <- function()
{
  cat(paste0(rep("\b", getOption("width")), collapse = ""))
  flush.console()
}

#############################################################################

repairSolution <- function(x, lo, up) 
{
# Repair solutions to the specified bound of decision variables
#
# Gilli, M., Maringer, D. & Schumann, E. (2011) Numerical Methods and 
#   Optimization in Finance, Academic Press, p. 551, algorithm 66
#
# x = n-length vector of solutions for n decision variables
# lo, up = n-length vector of upper and lower boundaries for each decision
#          variable
#
# Example:  
# set.seed(1)
# lo <- c(0, 0, 0, 0)
# up <- c(1, 1, 2, 1)
# x <- rnorm(4)
# x
# repairSolution(x, lo, up)

  xl <- lo - x
  xl <- xl + abs(xl)
  xu <- x - up
  xu <- xu + abs(xu)
  x <- x - (xu - xl)/2
  return(x)
}

#############################################################################

reflectSolution <- function(x, lo, up, tol = sqrt(.Machine$double.eps)) 
{
# Reflects solution values that are too large or too small around the boundary. 
# It restricts the change in a variable x[i] to the range up[i] - lo[i].
#  
# x = n-length vector of solutions for n decision variables
# lo, up = n-length vector of upper and lower boundaries for each decision
#          variable
#
# Example:
# set.seed(1)
# lo <- rep(0,4)
# up <- rep(1,4)
# x <- rnorm(4)
# x
# reflectSolution(x, lo, up)

  done <- TRUE
  e <- sum(x - up + abs(x - up) + lo - x + abs(lo - x)) 
  if(e > tol) done <- FALSE
  r <- up - lo
  while(!done) 
  { xu <- x - up
    xu <- xu + abs(xu)
    xu <- xu + r - abs(xu - r)
    xl <- lo - x
    xl <- xl + abs(xl)
    xl <- xl + r - abs(xl - r)
    x <- x - (xu - xl)/2
    e <- sum(x - up + abs(x - up) + lo - x + abs(lo - x)) 
    if(e < tol) done <- TRUE
  }
  return(x)  
}

#############################################################################

jet.colors <- function(n)
{
# Creates a palette of n colors beginning with dark blue, ranging through
# shades of blue, cyan, green, yellow and red, and ending with dark red. 
# This is inspired by the colormap 'jet' available in Matlab.
  palette <- colorRampPalette(c("#00007F", "blue", "#007FFF", 
                                "cyan", "#7FFF7F", "yellow", 
                                "#FF7F00", "red", "#7F0000"))
  palette(n)
}

spectral.colors <- function (n) 
{
  col <- c("#2B83BA", "#ABDDA4", "#FFFFBF", "#FDAE61", "#D7191C")
  # colors obtained as rev(brewer.pal(5, "Spectral"))
  palette <- grDevices::colorRampPalette(col)
  palette(n)
}

bl2gr.colors <- function (n) 
{
  palette <- grDevices::colorRampPalette(c("#084081", "#0868AC", "#2B8CBE", 
                                "#4EB3D3", "#7BCCC4", "#A8DDB5", 
                                "#CCEBC5", "#E0F3DB"), 
                              space = "Lab")
  palette(n)
}

persp3D <- function(x, y, z, theta = 30, phi = 20, d = 5, expand = 2/3, xlim = range(x, finite = TRUE), ylim = range(y, finite = TRUE), zlim = range(z, finite = TRUE), levels = pretty(zlim, nlevels), nlevels = 20, color.palette = jet.colors, border = NA, ticktype = "detailed", xlab = NULL, ylab = NULL, zlab = NULL, ...)
{
#----------------------------------------------------------------------------#  
# 3D plot, i.e. perspective plot, with different levels in different colors
#
# Example
# y <- x <- seq(-10, 10, length=60)
# f <- function(x,y) { r <- sqrt(x^2+y^2); 10 * sin(r)/r }
# z <- outer(x, y, f)
# persp3D(x, y, z, theta = 30, phi = 30, expand = 0.5)
# persp3D(x, y, z, color.palette = heat.colors, phi = 30, theta = 225, box = TRUE, border = NA, shade = .4)
# persp3D(x, y, z, color.palette = terrain.colors, phi = 30, theta = 225, box = FALSE, border = NA, shade = .4)
#
# x1 = seq(-3,3,length=50)
# x2 = seq(-3,3,length=50)
# y = function(x1, x2) sin(x1)+cos(x2)
# persp3D(x1, x2, outer(x1,x2,y), zlab="y", theta = 150, phi = 20, expand = 0.6)
#
#----------------------------------------------------------------------------#

  if(is.null(xlab)) 
     xlab <- if(!missing(x)) 
                deparse(substitute(x))
             else "X"
  if(is.null(ylab)) 
     ylab <- if(!missing(y)) 
                deparse(substitute(y))
             else "Y"
   if(is.null(zlab)) 
      zlab <- if(!missing(z)) 
                 deparse(substitute(z))
              else "Z"
  if(missing(z))
    { if(!missing(x)) 
        { if(is.list(x)) 
            { z <- x$z
              y <- x$y
              x <- x$x }
          else 
            { z <- x
              x <- seq.int(0, 1, length.out = nrow(z)) }
         }
      else stop("no 'z' matrix specified")
    }
  else if(is.list(x))
         { y <- x$y
           x <- x$x }
  if(any(diff(x) <= 0) || any(diff(y) <= 0)) 
     stop("increasing 'x' and 'y' values expected")

  # getting the value of the midpoint
  zz <- (z[-1,-1] + z[-1,-ncol(z)] + z[-nrow(z),-1] + z[-nrow(z),-ncol(z)])/4
  # set colors for levels
  cols <- color.palette(length(levels)-1)
  zzz <- cut(zz, breaks = levels, labels = cols)
  # plot
  out <- persp(x, y, z, theta = theta, phi = phi, d = d, expand = expand,
               col = as.character(zzz),
               xlim = xlim, ylim = ylim, zlim = zlim,
               border = border, ticktype = ticktype, 
               xlab = xlab, ylab = ylab, zlab = zlab, ...)
  # add breaks and colors for a legend
  out <- list(persp = out, levels = levels, colors = cols)
  invisible(out)
}

#----------------------------------------------------------------------------#
# print a short version of a matrix by allowing to select the number of 
# head/tail rows and columns to display

.printShortMatrix <- function(x, head = 2, tail = 1, chead = 5, ctail = 1, ...)
{ 
  x <- as.matrix(x)
  nr <- nrow(x)
  nc <- ncol(x)
  if(is.na(head <- as.numeric(head))) head <- 2
  if(is.na(tail <- as.numeric(tail))) tail <- 1
  if(is.na(chead <- as.numeric(chead))) chead <- 5
  if(is.na(ctail <- as.numeric(ctail))) ctail <- 1
  
  if(nr > (head + tail + 1))
    { rnames <- rownames(x)
      if(is.null(rnames)) 
        rnames <- paste("[", 1:nr, ",]", sep ="")
      x <- rbind(x[1:head,,drop=FALSE], 
                 rep(NA, nc), 
                 x[(nr-tail+1):nr,,drop=FALSE])
      rownames(x) <- c(rnames[1:head], "...", rnames[(nr-tail+1):nr])
  }
  if(nc > (chead + ctail + 1))
    { cnames <- colnames(x)
      if(is.null(cnames)) 
        cnames <- paste("[,", 1:nc, "]", sep ="")
      x <- cbind(x[,1:chead,drop=FALSE], 
                 rep(NA, nrow(x)), 
                 x[,(nc-ctail+1):nc,drop=FALSE])
      colnames(x) <- c(cnames[1:chead], "...", cnames[(nc-ctail+1):nc])
  }
          
  print(x, na.print = "", ...)
}
