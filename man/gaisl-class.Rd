\name{gaisl-class}
\docType{class}
\alias{gaisl-class}

\title{Class "gaisl"}

\description{An S4 class for islands genetic algorithms (ISLGAs)}

\section{Objects from the Class}{
Objects can be created by calls to the \code{\link{gaisl}} function.
}

\section{Slots}{
  \describe{
    \item{\code{call}}{an object of class \code{"call"} representing the matched call;}
    \item{\code{type}}{a character string specifying the type of genetic algorithm used;}
    \item{\code{min}}{a vector providing for each decision variable the minimum of the search space in case of real-valued or permutation encoded optimizations;}
    \item{\code{max}}{a vector providing for each decision variable the maximum of the search space in case of real-valued or permutation encoded optimizations;}
    \item{\code{nBits}}{a value specifying the number of bits to be used in binary encoded optimizations;}
    \item{\code{names}}{a vector of character strings providing the names of decision variables (optional);}
    \item{\code{popSize}}{the population size;}
    \item{\code{numIslands}}{the number of islands;}
    \item{\code{migrationRate}}{the migration rate;}
    \item{\code{migrationInterval}}{the migration interval;}
    \item{\code{maxiter}}{the maximum number of ISLGA iterations before the search is halted;}
    \item{\code{run}}{the number of consecutive generations without any improvement in the best fitness value before the ISLGA is stopped;}
    \item{\code{maxiter}}{the maximum number of iterations to run before the GA search is halted;}
    \item{\code{suggestions}}{a matrix of user provided solutions and included in the initial population;}
    \item{\code{elitism}}{the number of best fitness individuals to survive at each generation;}
    \item{\code{pcrossover}}{the crossover probability;}
    \item{\code{pmutation}}{the mutation probability;}
    \item{\code{islands}}{a list containing the "ga" objects corresponding to each island GA evolution;}
    \item{\code{summary}}{a list of matrices of summary statistics for fitness values at each iteration (along the rows). Each element of the list corresponds to the evolution of an island;}
    \item{\code{fitnessValues}}{a list of best fitness values found in each island at the final iteration;}
    \item{\code{solutions}}{a list of matrices, one for each island, containing the values of the decision variables giving the best fitness at the final iteration;}
    \item{\code{fitnessValue}}{the best fitness value at the final iteration;}
    \item{\code{solution}}{a matrix containing the values of the decision variables giving the best fitness at the final iteration.}
}}

\author{Luca Scrucca}

\keyword{classes}

\seealso{For examples of usage see \code{\link{gaisl}}.}

%\examples{}

\keyword{classes}
\keyword{optimize}


