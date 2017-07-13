\name{GA-package}
\alias{GA-package}
%\alias{ga-package}
\alias{GA}
\docType{package}

\title{Genetic Algorithms}

\description{
An R package for optimization using genetic algorithms. The package provides a flexible general-purpose set of tools for implementing genetic algorithms search in both the continuous and discrete case, whether constrained or not. Users can easily define their own objective function depending on the problem at hand. Several genetic operators are available and can be combined to explore the best settings for the current task. Furthermore, users can define new genetic operators and easily evaluate their performances. Local search using general-purpose optimisation algorithms can be applied stochastically to exploit interesting regions. GAs can be run sequentially or in parallel, using an explicit master-slave parallelisation or a coarse-grain islands approach.
}

\details{
For a quick intro to GA package see the \href{../doc/index.html}{vignette} accompanying the package. Further details are provided in the papers referenced below.
%GA: A Package for Genetic Algorithms in R \url{../doc/GA.pdf}
%Using parallel computing in GA package  \url{../doc/GAparallel.pdf}
}

\author{Luca Scrucca \email{luca.scrucca@unipg.it}}

\section{References}{
Scrucca, L. (2013) GA: A Package for Genetic Algorithms in R. \emph{Journal of Statistical Software}, 53(4), 1-37. \url{http://www.jstatsoft.org/v53/i04/}

Scrucca, L. (2017) On some extensions to GA package: hybrid optimisation, parallelisation and islands evolution. \emph{The R Journal}, 9/1, 187–206. \url{https://journal.r-project.org/archive/2017/RJ-2017-008}
}

% \examples{}

\keyword{package}
