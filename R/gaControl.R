#-------------------------------------------------------------------#
# Function for options retrieval and setting
#-------------------------------------------------------------------#

gaControl <- function(...)
{
  current <- .ga.default
  if(nargs() == 0) return(current)
  args <- list(...)
  if(length(args) == 1 && is.null(names(args))) 
    { arg <- args[[1]]
      switch(mode(arg),
             list = args <- arg,
             character = return(.ga.default[[arg]]),
             stop("invalid argument: ", dQuote(arg)))
    }

  if(length(args) == 0) return(current)
  nargs <- names(args)
  if(is.null(nargs)) stop("options must be given by name")

  if(is.list(args))
    { changed <- current[nargs]
      for(i in 1:length(nargs))
         { if(is.list(args[[i]]))
             { what <- names(args[[i]])
               changed[[i]][what] <- args[[i]][what] }
            else
              { changed[i] <- args[[i]] }
         }
      current[nargs] <- changed
    }
  else
    { changed <- current[nargs]
      current[nargs] <- args
    }

  if(sys.parent() == 0) env <- asNamespace("GA") else env <- parent.frame()
  assign(".ga.default", current, envir = env)
  invisible(current)
}
    
#  call <- match.call(expand.dots = TRUE)
#  # use version in the globalenv, if exists (needed due to Namespace)
#  if(exists(".ga.default", where = .GlobalEnv)) 
#    { .ga.default <- get(".ga.default", pos = .GlobalEnv) }
#
#  current <- .ga.default
#  # current <- getFromNamespace(".ga.default", "GA")
#  args <- list(...)
#  if(length(args) == 0) return(current)
#
#  if(length(args) == 1 && is.null(names(args))) 
#    { arg <- match.arg(args[[1]], names(.ga.default))
#      switch(mode(arg),
#             list = temp <- arg,
#             character = return(.ga.default[[arg]]),
#             stop("invalid argument: ", sQuote(arg)))
#    }
#
#  n <- names(args)
#  if(is.null(n)) stop("options must be given by name")
#  if(is.list(args[[1]]))
#    { changed <- current[[names(args)]]
#      what <- sapply(args, names) 
#      changed[what] <- args[[1]][what]
#      current[[names(args)]] <- changed
#    }
#  else
#    { changed <- current[names(args)]
#      what <- names(args) 
#      changed[what] <- args[what]
#      current[names(args)] <- changed
#    }
#
#  ## This assigns .ga.default in the global environment.  That way one
#  ## can get back to the `factory defaults' by removing the variable from
#  ## the global environment.  It also means that options are remembered
#  ## between sessions (if the environment is saved).  
#  assign(".ga.default", current, pos = .GlobalEnv)
#  invisible(current)


.ga.default <- list("binary" = list(population = "gabin_Population",
                                    selection  = "gabin_lrSelection",
                                    crossover  = "gabin_spCrossover",
                                    mutation   = "gabin_raMutation"),
                    "real-valued" = list(population = "gareal_Population",
                                         selection  = "gareal_lsSelection",
                                         crossover  = "gareal_laCrossover",
                                         mutation   = "gareal_raMutation"),                                    
                    "permutation" = list(population = "gaperm_Population",
                                         selection  = "gaperm_lrSelection",
                                         crossover  = "gaperm_oxCrossover",
                                         mutation   = "gaperm_simMutation"),
                    "eps" = sqrt(.Machine$double.eps)                     
                   )
