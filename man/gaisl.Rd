\name{gaisl}
\alias{gaisl}
\alias{show,gaisl-method}
\alias{print,gaisl-method}

\title{Islands Genetic Algorithms}

\description{Maximization of a fitness function using islands genetic algorithms (ISLGAs). This is a distributed multiple-population GA, where the population is partitioned into several subpopulations and assigned to separated islands. Independent GAs are executed in each island, and only occasionally sparse exchanges of individuals are performed among the islands. In principle islands can evolve sequentially, but increased computational efficiency is obtained by running GAs in parallel on each island. The latter is called island parallel GAs (ISLPGAs) and it is used by default.}

\usage{
gaisl(type = c("binary", "real-valued", "permutation"), 
      fitness, \dots,
      min, max, nBits,
      population = gaControl(type)$population,
      selection = gaControl(type)$selection,
      crossover = gaControl(type)$crossover, 
      mutation = gaControl(type)$mutation,
      popSize = 100, 
      numIslands = 4, 
      migrationRate = 0.1, 
      migrationInterval = 10,
      pcrossover = 0.8, 
      pmutation = 0.1, 
      elitism = base::max(1, round(popSize/numIslands*0.05)), 
      updatePop = FALSE,
      postFitness = NULL,
      maxiter = 1000,
      run = maxiter,
      maxFitness = Inf,
      names = NULL,
      suggestions = NULL, 
      optim = FALSE,
      optimArgs = list(method = "L-BFGS-B", 
                       poptim = 0.05,
                       pressel = 0.5,
                       control = list(fnscale = -1, maxit = 100)),
      parallel = TRUE,
      monitor = if(interactive()) 
                  { if(is.RStudio()) gaislMonitor else gaislMonitor2 } 
                else FALSE,
      seed = NULL)
}

\arguments{
\item{type}{the type of genetic algorithm to be run depending on the nature of decision variables. Possible values are: 

\tabular{ll}{
\code{"binary"} \tab for binary representations of decision variables. \cr
\code{"real-valued"} \tab for optimization problems where the decision variables are floating-point representations of real numbers. \cr
\code{"permutation"} \tab for problems that involves reordering of a list.
}
}

\item{fitness}{the fitness function, any allowable R function which takes as input an individual \code{string} representing a potential solution, and returns a numerical value describing its ``fitness''.}

\item{\dots}{additional arguments to be passed to the fitness function. This allows to write fitness functions that keep some variables fixed during the search.}

\item{min}{a vector of length equal to the decision variables providing the minimum of the search space in case of real-valued or permutation encoded optimizations.}

\item{max}{a vector of length equal to the decision variables providing the maximum of the search space in case of real-valued or permutation encoded optimizations.}

\item{nBits}{a value specifying the number of bits to be used in binary encoded optimizations.}

\item{population}{an R function for randomly generating an initial population. See \code{\link{ga_Population}} for available functions.}

\item{numIslands}{an integer value specifying the number of islands to be used in a \emph{ring topology}, in which each island is connected unidirectionally with another island, hence forming a single continuous pathway.}
 
\item{migrationRate}{a value in the range $[0,1]$ providing the proportion of individuals that should migrate between the islands.}

\item{migrationInterval}{an integer value specifying the number of iterations at which exchange of individuals takes place.}

\item{selection}{an R function performing selection, i.e. a function which generates a new population of individuals from the current population probabilistically according to individual fitness. See \code{\link{ga_Selection}} for available functions.}

\item{crossover}{an R function performing crossover, i.e. a function which forms offsprings by combining part of the genetic information from their parents. See \code{\link{ga_Crossover}} for available functions.}

\item{mutation}{an R function performing mutation, i.e. a function which randomly alters the values of some genes in a parent chromosome. See \code{\link{ga_Mutation}} for available functions.}

\item{popSize}{the population size.}

\item{updatePop}{a logical defaulting to \code{FALSE}. If set at \code{TRUE} the first attribute attached to the value returned by the user-defined fitness function is used to update the population. \cr
\emph{Be careful though, this is an experimental feature!}}

\item{postFitness}{a user-defined function which, if provided, receives the current \code{ga-class} object as input, performs post fitness-evaluation steps, then returns an updated version of the object which is used to update the GA search.\cr
\emph{Be careful though, this is an experimental feature!}}

\item{pcrossover}{the probability of crossover between pairs of chromosomes. Typically this is a large value and by default is set to 0.8.}

\item{pmutation}{the probability of mutation in a parent chromosome. Usually mutation occurs with a small probability, and by default is set to 0.1.}

\item{elitism}{the number of best fitness individuals to survive at each generation. By default the top 5\% individuals in each island will survive at each iteration.}

\item{maxiter}{the maximum number of iterations to run before the GA search is halted.}

\item{run}{the number of consecutive generations without any improvement in the best fitness value before the GA is stopped.}

\item{maxFitness}{the upper bound on the fitness function after that the GA search is interrupted.}

\item{names}{a vector of character strings providing the names of decision variables.}

\item{suggestions}{a matrix of solutions strings to be included in the initial population. If provided the number of columns must match the number of decision variables.}

\item{optim}{a logical defaulting to \code{FALSE} determining whether or not a local search using general-purpose optimisation algorithms should be used. See argument \code{optimArgs} for further details and finer control.}

\item{optimArgs}{a list controlling the local search algorithm with the following components:

\describe{
\item{\code{method}}{a string specifying the general-purpose optimisation method to be used, by default is set to \code{"L-BFGS-B"}. Other possible methods are those reported in \code{\link[stats]{optim}}.}
\item{\code{poptim}}{a value in the range [0,1] specifying the probability of performing a local search at each iteration of GA (default 0.1).}
\item{\code{pressel}}{a value in the range [0,1] specifying the pressure selection (default 0.5). The local search is started from a random solution selected with probability proportional to fitness. High values of \code{pressel} tend to select the solutions with the largest fitness, whereas low values of \code{pressel} assign quasi-uniform probabilities to any solution.}
\item{\code{control}}{a list of control parameters. See 'Details' section in \code{\link[stats]{optim}}.}
}
}

\item{parallel}{An optional argument which allows to specify if the Islands Genetic Algorithm should be run sequentially or in parallel. 

For a single machine with multiple cores, possible values are: 
\itemize{
  \item a logical value specifying if parallel computing should be used (\code{TRUE}) or not (\code{FALSE}, default) for running GAs on each island;
  \item a numerical value which gives the number of cores to employ. By default, this is obtained from the function \code{\link[parallel]{detectCores}}; 
  \item a character string specifying the type of parallelisation to use. This depends on system OS: on Windows OS only \code{"snow"} type functionality is available, while on Unix/Linux/Mac OSX both \code{"snow"} and \code{"multicore"} (default) functionalities are available.
}
In all the cases described above, at the end of the search the cluster is automatically stopped by shutting down the workers.

If a cluster of multiple machines is available, evolution of GAs on each island can be executed in parallel using all, or a subset of, the cores available to the machines belonging to the cluster. However, this option requires more work from the user, who needs to set up and register a parallel back end. 
In this case the cluster must be explicitely stopped with \code{\link[parallel]{stopCluster}}.
}

\item{monitor}{a logical or an R function which takes as input the current state of the \code{gaisl-class} object and show the evolution of the search in different epochs. By default, for interactive sessions, the function \code{\link{gaislMonitor}} or \code{\link{gaislMonitor2}}, depending on whether or not is an RStudio session, prints the average and best fitness values at each epoch for each island. In non interactive sessions, by default \code{monitor = FALSE} so any output is suppressed.}

\item{seed}{an integer value containing the random number generator state. This argument can be used to replicate the results of a ISLGA search. Note that if parallel computing is required, the \pkg{doRNG} package must be installed.}
}

\details{
Genetic algorithms (GAs) are stochastic search algorithms inspired by the basic principles of biological evolution and natural selection. GAs simulate the evolution of living organisms, where the fittest individuals dominate over the weaker ones, by mimicking the biological mechanisms of evolution, such as selection, crossover and mutation.

The \code{gaisl} function implements the islands GAs approach, where the population is partitioned into several subpopulations and assigned to separated islands. Independent GAs are executed in each island, and only occasionally sparse exchanges of individuals are performed among the islands. The algorithm can be run in parallel or sequentially. 
For more information on GAs see \code{\link{ga}}. 
}
  
\value{Returns an object of class \code{gaisl-class}. See \code{\link{gaisl-class}} for a description of available slots information.}

\references{

Luque G., Alba E. (2011) \emph{Parallel Genetic Algorithms: Theory and Real World Applications}. Springer.

Luke S. (2013) \emph{Essentials of Metaheuristics}, 2nd edition. Lulu. Freely available at \url{http://cs.gmu.edu/~sean/book/metaheuristics/}.

Scrucca, L. (2017) On some extensions to GA package: hybrid optimisation, parallelisation and islands evolution. \emph{The R Journal}, 9/1, 187–206. \url{https://journal.r-project.org/archive/2017/RJ-2017-008}
}
  
\author{Luca Scrucca \email{luca.scrucca@unipg.it}}

%\note{}

\seealso{
  \code{\link{summary,gaisl-method}}, 
  \code{\link{plot,gaisl-method}}, 
  \code{\link{gaisl-class}},
  \code{\link{ga}}
}

\examples{
\dontrun{
# two-dimensional Rastrigin function
Rastrigin <- function(x1, x2)
{
  20 + x1^2 + x2^2 - 10*(cos(2*pi*x1) + cos(2*pi*x2))
}

x1 <- x2 <- seq(-5.12, 5.12, by = 0.1)
f <- outer(x1, x2, Rastrigin)
persp3D(x1, x2, f, theta = 50, phi = 20)
filled.contour(x1, x2, f, color.palette = jet.colors)

GA <- gaisl(type = "real-valued", 
            fitness =  function(x) -Rastrigin(x[1], x[2]),
            min = c(-5.12, -5.12), max = c(5.12, 5.12), 
            popSize = 80, maxiter = 500, 
            numIslands = 4, migrationInterval = 50)
summary(GA)
plot(GA)
}
}

\keyword{optimize}
