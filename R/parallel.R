startParallel <- function(parallel = TRUE, ...)
{
# Start parallel computing for GA package
  
  # check availability of parallel and doParallel (their dependencies, i.e. 
  # foreach and iterators, are specified as Depends on package DESCRIPTION file)
  if(!all(requireNamespace("parallel", quietly = TRUE),
        requireNamespace("doParallel", quietly = TRUE)))     
     stop("packages 'parallel' and 'doParallel' required for parallelization!")

  # if a cluster is provided as input argument use that cluster and exit
  if(any(class(parallel) == "cluster"))
    { cl <- parallel
      parallel <- TRUE
      attr(parallel, "type") <- getDoParName()
      attr(parallel, "cores") <- getDoParWorkers()
      attr(parallel, "cluster") <- cl
      return(parallel)
  }
    
  # set default parallel functionality depending on system OS:
  # - snow functionality on Windows OS
  # - multicore functionality on Unix-like systems (Unix/Linux & Mac OSX)
  parallelType <- if(.Platform$OS.type == "windows") 
                    "snow" else "multicore"

  # get the current number of cores available
  numCores <- parallel::detectCores()

  # set parameters for parallelization
  if(is.logical(parallel))
    { NULL }
  else if(is.numeric(parallel))
    { numCores <- as.integer(parallel)
      parallel <- TRUE }
  else if(is.character(parallel))
    { parallelType <- parallel
      parallel <- TRUE 
    }
  else parallel <- FALSE
  
  attr(parallel, "type") <- parallelType
  attr(parallel, "cores") <- numCores

  # start "parallel backend" if needed
  if(parallel)
  { 
    if(parallelType == "snow")
      { 
        # snow functionality on Unix-like systems & Windows
        cl <- parallel::makeCluster(numCores, type = "PSOCK")
        attr(parallel, "cluster") <- cl
        # export parent environment
        varlist <- ls(envir = parent.frame(), all.names = TRUE)
        varlist <- varlist[varlist != "..."]
        parallel::clusterExport(cl, varlist = varlist,
                                # envir = parent.env(environment())
                                envir = parent.frame() )
        # export global environment (workspace)
        parallel::clusterExport(cl, 
                                varlist = ls(envir = globalenv(), 
                                             all.names = TRUE),
                                envir = globalenv())
        # load current packages in workers
        pkgs <- .packages()
        lapply(pkgs, function(pkg) 
               parallel::clusterCall(cl, library, package = pkg, 
                                     character.only = TRUE))
        #
        doParallel::registerDoParallel(cl, cores = numCores)
      }
      else if(parallelType == "multicore")
        { # multicore functionality on Unix-like systems
          cl <- parallel::makeCluster(numCores, type = "FORK")
          doParallel::registerDoParallel(cl, cores = numCores) 
          attr(parallel, "cluster") <- cl
        }
      else 
        { stop("Only 'snow' and 'multicore' clusters allowed!") }
  }

  return(parallel)
}

stopParallel <- function(cluster, ...)
{ 
# Stop parallel computing for GA package
  parallel::stopCluster(cluster)
  foreach::registerDoSEQ()
  invisible()
}