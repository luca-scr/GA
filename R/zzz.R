.onAttach <- function(lib, pkg)
{
  unlockBinding(".ga.default", asNamespace("GA")) 
  version <- read.dcf(file.path(lib, pkg, "DESCRIPTION"), "Version")
  packageStartupMessage("Package 'GA' version ", version)
  packageStartupMessage("Type 'citation(\"GA\")' for citing this R package in publications.")
  invisible()
}



  