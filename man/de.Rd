\name{de}
\alias{de}
\alias{show,de-method}
\alias{print,de-method}

\title{Differential Evolution via Genetic Algorithms}

\description{Maximization of a fitness function using Differential Evolution (DE).
DE is a population-based evolutionary algorithm for optimisation of fitness functions defined over a continuous parameter space.}

\usage{
de(fitness,
   lower, upper,
   popSize = 10*d,
   stepsize = 0.8,
   pcrossover = 0.5,
   \dots)
}

\arguments{
\item{fitness}{the fitness function, any allowable R function which takes as input a vector of values representing a potential solution, and returns a numerical value describing its ``fitness''.}

\item{lower}{a vector of length equal to the decision variables providing the lower bounds of the search space.}

\item{upper}{a vector of length equal to the decision variables providing the upper bounds of the search space.}

\item{popSize}{the population size. By default is set at 10 times the number of decision variables.}

\item{pcrossover}{the probability of crossover, by default set to 0.5.}

\item{stepsize}{the stepsize or weighting factor. A value in the interval [0,2], by default set to 0.8. If set at \code{NA} a random value is selected in the interval [0.5, 1.0] (so called dithering).}

\item{\dots}{additional arguments to be passed to the fitness function. This allows to write fitness functions that keep some variables fixed during the search.}
}

\details{
Differential Evolution (DE) is a stochastic evolutionary algorithm that optimises multidimensional real-valued fitness functions without requiring the optimisation problem to be differentiable.

This implimentation follows the description in Simon (2013; Sec. 12.4, and Fig. 12.12) and uses the functionalities available in the \code{\link{ga}} function for Genetic Algorithms. 

The DE selection operator is defined by \code{\link{gareal_de}} with parameters \code{p = pcrossover} and \code{F = stepsize}.
}
  
\value{Returns an object of class \code{de-class}. See \code{\link{de-class}} for a description of available slots information.}

\references{

Scrucca L. (2013). GA: A Package for Genetic Algorithms in R.
\emph{Journal of Statistical Software}, 53(4), 1-37, \url{http://www.jstatsoft.org/v53/i04/}.

Simon D. (2013) \emph{Evolutionary Optimization Algorithms}. John Wiley & Sons.

Price K., Storn R.M., Lampinen J.A. (2005) \emph{Differential Evolution: A Practical Approach to Global Optimization}. Springer. 
}
  
\author{Luca Scrucca \email{luca.scrucca@unipg.it}}

%\note{}

\seealso{
  \code{\link{summary,de-method}}, 
  \code{\link{plot,de-method}}, 
  \code{\link{de-class}}
}

\examples{
# 1) one-dimensional function
f <- function(x)  abs(x)+cos(x)
curve(f, -20, 20)

DE <- de(fitness = function(x) -f(x), lower = -20, upper = 20)
plot(DE)
summary(DE)

curve(f, -20, 20, n = 1000)
abline(v = DE@solution, lty = 3)

# 2) "Wild" function, global minimum at about -15.81515

wild <- function(x) 10*sin(0.3*x)*sin(1.3*x^2) + 0.00001*x^4 + 0.2*x + 80
plot(wild, -50, 50, n = 1000)

# from help("optim")
SANN <- optim(50, fn = wild, method = "SANN",
              control = list(maxit = 20000, temp = 20, parscale = 20))
unlist(SANN[1:2])

DE <- de(fitness = function(...) -wild(...), lower = -50, upper = 50)
plot(DE)
summary(DE)

# 3) two-dimensional Rastrigin function

Rastrigin <- function(x1, x2)
{
  20 + x1^2 + x2^2 - 10*(cos(2*pi*x1) + cos(2*pi*x2))
}

x1 <- x2 <- seq(-5.12, 5.12, by = 0.1)
f <- outer(x1, x2, Rastrigin)
persp3D(x1, x2, f, theta = 50, phi = 20, col.palette = bl2gr.colors)

DE <- de(fitness = function(x) -Rastrigin(x[1], x[2]),
         lower = c(-5.12, -5.12), upper = c(5.12, 5.12),
         popSize = 50)
plot(DE)
summary(DE)

filled.contour(x1, x2, f, color.palette = bl2gr.colors,
               plot.axes = { axis(1); axis(2); 
                             points(DE@solution, 
                                    col = "yellow", pch = 3, lwd = 2) })

# 4) two-dimensional Ackley function

Ackley <- function(x1, x2)
{
  -20*exp(-0.2*sqrt(0.5*(x1^2 + x2^2))) - 
  exp(0.5*(cos(2*pi*x1) + cos(2*pi*x2))) + exp(1) + 20
}

x1 <- x2 <- seq(-3, 3, by = 0.1)
f <- outer(x1, x2, Ackley)
persp3D(x1, x2, f, theta = 50, phi = 20, col.palette = bl2gr.colors)

DE <- de(fitness = function(x) -Ackley(x[1], x[2]),
         lower = c(-3, -3), upper = c(3, 3),
         stepsize = NA)
plot(DE)
summary(DE)

filled.contour(x1, x2, f, color.palette = bl2gr.colors,
               plot.axes = { axis(1); axis(2); 
                             points(DE@solution, 
                                    col = "yellow", pch = 3, lwd = 2) })
}

\keyword{optimize}
