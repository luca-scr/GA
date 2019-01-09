\name{de-class}
\docType{class}
\alias{de-class}

\title{Class "de"}

\description{An S4 class for differential evolution algorithm}

\section{Objects from the Class}{
Objects can be created by calls to the \code{\link{de}} function.
}

\section{Slots}{
  \describe{
    \item{\code{call}}{an object of class \code{"call"} representing the matched call;}
    \item{\code{type}}{a character string specifying the type of genetic algorithm used;}
    \item{\code{lower}}{a vector providing for each decision variable the lower bounds of the search space in case of real-valued or permutation encoded optimisations. Formerly this slot was named \code{min};}
    \item{\code{upper}}{a vector providing for each decision variable the upper bounds of the search space in case of real-valued or permutation encoded optimizations. Formerly this slot was named \code{max};}
    \item{\code{names}}{a vector of character strings providing the names of decision variables (optional);}
    \item{\code{popSize}}{the population size;}
    \item{\code{iter}}{the actual (or final) iteration of DE search;}
    \item{\code{run}}{the number of consecutive generations without any improvement in the best fitness value before the DE is stopped;}
    \item{\code{maxiter}}{the maximum number of iterations to run before the DE search is halted;}
    \item{\code{suggestions}}{a matrix of user provided solutions and included in the initial population;}
    \item{\code{population}}{the current (or final) population;}
    \item{\code{elitism}}{the number of best fitness individuals to survive at each generation;}
    \item{\code{stepsize}}{the stepsize or weighting factor;}
    \item{\code{pcrossover}}{the crossover probability;}
    \item{\code{pmutation}}{the mutation probability;}
    \item{\code{optim}}{a logical specifying whether or not a local search using general-purpose optimisation algorithms should be used;}
    \item{\code{fitness}}{the values of fitness function for the current (or final) population;}
    \item{\code{summary}}{a matrix of summary statistics for fitness values at each iteration (along the rows);}
    \item{\code{bestSol}}{if \code{keepBest = TRUE}, the best solutions at each iteration;}
    \item{\code{fitnessValue}}{the best fitness value at the final iteration;}
    \item{\code{solution}}{the value(s) of the decision variables giving the best fitness at the final iteration.}
  }
}

\author{Luca Scrucca}

\keyword{classes}

\seealso{For examples of usage see \code{\link{de}}.}

%\examples{}

\keyword{classes}
\keyword{optimize}


