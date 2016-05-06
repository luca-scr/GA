.onAttach <- function(lib, pkg)
{
  unlockBinding(".ga.default", asNamespace("GA")) 
  version <- read.dcf(file.path(lib, pkg, "DESCRIPTION"), "Version")
  
  if(interactive())
    { # > figlet GA
      packageStartupMessage(
"  ____    _    
 / ___|  / \\     Genetic 
| |  _  / _ \\    Algorithms
| |_| |/ ___ \\   
 \\____/_/   \\_\\  version ", version)
}
else
  { packageStartupMessage("Package 'GA' version ", version) } 

  packageStartupMessage("Type 'citation(\"GA\")' for citing this R package in publications.")
  invisible()
}
  