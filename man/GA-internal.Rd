\name{GA-internal}
\alias{.printShortMatrix}
\alias{startParallel}
\alias{stopParallel}
\alias{garun}
\alias{optimProbsel}
\alias{reflectSolution}
\alias{repairSolution}
\alias{GAStartupMessage}

\title{Internal GA functions}

\description{
Internal functions not intended to be called directly by users.
}

\usage{
.printShortMatrix(x, head = 2, tail = 1, chead = 5, ctail = 1, \dots)
startParallel(parallel = TRUE, \dots) 
stopParallel(cluster, \dots)
garun(x)
optimProbsel(x, q = 0.25, \dots)
reflectSolution(x, lo, up, tol = sqrt(.Machine$double.eps))
repairSolution(x, lo, up) 
GAStartupMessage()
}

\keyword{internal}
