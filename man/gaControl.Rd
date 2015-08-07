\name{gaControl}
\alias{gaControl}

\title{A function for setting or retrieving defaults genetic operators}

\description{Default settings for genetic operators used in the GA package.}

\usage{
gaControl(\dots)
}

\arguments{
  \item{\dots}{
  no arguments, a single character vector, or a named list with components.
}}

\details{
If the function is called with no arguments returns the current default settings, i.e., a list with the following default components:
\itemize{
\item{\code{"binary"}}{
  \itemize{
  \item{\code{population}} = \code{"gabin_Population"}
  \item{\code{selection}} = \code{"gabin_lrSelection"}
  \item{\code{crossover}} = \code{"gabin_spCrossover"}
  \item{\code{mutation}} = \code{"gabin_raMutation"}
  }
}
\item{\code{"real-valued"}}{
  \itemize{
  \item{\code{population}} = \code{"gareal_Population"}
  \item{\code{selection}} = \code{"gareal_lsSelection"}
  \item{\code{crossover}} = \code{"gareal_laCrossover"}
  \item{\code{mutation}} = \code{"gareal_raMutation"}
  }
}
\item{\code{"permutation"}}{
  \itemize{
  \item{\code{population}} = \code{"gaperm_Population"}
  \item{\code{selection}} = \code{"gaperm_lrSelection"}
  \item{\code{crossover}} = \code{"gaperm_oxCrossover"}
  \item{\code{mutation}} = \code{"gaperm_simMutation"}
  }
}

\item{\code{"eps"} = }{the tolerance value used by the package functions. By default set at \code{sqrt(.Machine$double.eps)}.}
}

The function may be called with a single string specifying the name of the component. In this case the function returns the current default settings.

To change the default values, a named component must be followed by a single value (in case of \code{"eps"}) or a list of component(s) specifying the name of the function for a genetic operator. See the Examples section.
}

\value{
If the argument list is empty the function returns the current list of values.  
If the argument list is not empty, the returned list is invisible.
}

\note{
The parameter values set via a call to this function will remain in effect for the rest of the session, affecting the subsequent behaviour of the functions for which the given parameters are relevant.
}

\author{Luca Scrucca}

\seealso{\code{\link{ga}}}

\examples{
# get and save defaults
defaultControl <- gaControl()
print(defaultControl)
# get current defaults only for binary search
gaControl("binary")
# set defaults for selection operator of binary search
gaControl("binary" = list(selection = "gabin_tourSelection"))
gaControl("binary")
# set defaults for selection and crossover operators of binary search
gaControl("binary" = list(selection = "ga_rwSelection", 
                          crossover = "gabin_uCrossover"))
gaControl("binary")
# restore defaults
gaControl(defaultControl)
gaControl()
}
