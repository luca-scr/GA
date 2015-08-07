\name{ga-class}
\docType{class}
\alias{ga-class}

\title{Class "ga"}

\description{An S4 class for genetic algorithms}

\section{Objects from the Class}{
Objects can be created by calls to the \code{\link{ga}} function.
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
    \item{\code{iter}}{the actual (or final) iteration of GA search;}
    \item{\code{run}}{the number of consecutive generations without any improvement in the best fitness value before the GA is stopped;}
    \item{\code{maxiter}}{the maximum number of iterations to run before the GA search is halted;}
    \item{\code{suggestions}}{a matrix of user provided solutions and included in the initial population;}
    \item{\code{population}}{the current (or final) population;}
    \item{\code{elitism}}{the number of best fitness individuals to survive at each generation;}
    \item{\code{pcrossover}}{the crossover probability;}
    \item{\code{pmutation}}{the mutation probability;}
    \item{\code{fitness}}{the values of fitness function for the current (or final) population;}
    \item{\code{summary}}{a matrix of summary statistics for fitness values at each iteration (along the rows);}
    \item{\code{bestSol}}{if \code{keepBest = TRUE}, the best solutions at each iteration;}
    \item{\code{fitnessValue}}{the best fitness value at the final iteration;}
    \item{\code{solution}}{the value(s) of the decision variables giving the best fitness at the final iteration.}
  }
}

\author{Luca Scrucca}

\keyword{classes}

\seealso{For examples of usage see \code{\link{ga}}.}

%\examples{}

\keyword{classes}
\keyword{optimize}


