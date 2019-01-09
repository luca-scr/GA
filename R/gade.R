##############################################################################
#                                                                            #
#         DIFFERENTIAL EVOLUTION via GENETIC ALGORITHMS in R                 #
#                                                                            #
##############################################################################

# Note: this file must be named to alphabetically follow ga.R

de <- function(fitness,
               lower, upper,
               popSize = 10*d,
               stepsize = 0.8,
               pcrossover = 0.5,
               ...) 
{
  call <- match.call()
  args <- list(...)
  args$type <- "real-valued"
  args$nBits <- NULL
  # DE selection including crossover
  args$selection <- function(...) 
    gareal_de(..., F = stepsize, p = pcrossover)
  args$pcrossover <- 0 # skip GA crossover
  if(is.null(args$elitism))
    args$elitism <- 0
  if(is.null(args$pmutation))
    args$pmutation <- 0 
  if(is.null(args$optim))
    args$optim <- FALSE
  if(is.null(args$monitor) & interactive())
    args$monitor <- deMonitor 
  
  lower <- as.vector(lower)
  upper <- as.vector(upper)
  stopifnot(length(lower) == length(upper))
  d <- length(lower)
  popSize <- as.numeric(popSize)
  
  object <- do.call("ga", c(args, 
                            list(fitness = fitness, 
                                 lower = lower, 
                                 upper = upper,
                                 popSize = popSize)))
  object <- as(object, "de")
  object@call <- call
  object@pcrossover <- pcrossover
  object@stepsize <- stepsize
  
  return(object)
}

setClass(Class = "de", 
         representation(call = "language",
                        type = "character",
                        lower = "numericOrNA", 
                        upper = "numericOrNA", 
                        names = "character",
                        popSize = "numeric",
                        iter = "numeric", 
                        run = "numeric", 
                        maxiter = "numeric",
                        suggestions = "matrix",
                        population = "matrix",
                        elitism = "numeric", 
                        stepsize = "numericOrNA",
                        pcrossover = "numeric",
                        pmutation = "numeric", 
                        optim = "logical",
                        fitness = "numericOrNA",
                        summary = "matrix",
                        bestSol = "list",
                        fitnessValue = "numeric",
                        solution = "matrix"
                      ),
         package = "GA" 
) 

# register conversion function
setAs("ga", "de", function(from, to) 
{
  to <- new(to)
  for (n in slotNames(from)) 
  {
    if(.hasSlot(to, n))
      slot(to, n) <- slot(from, n)
  }
  to
})

setMethod("print", "de", function(x, ...) str(x))

setMethod("show", "de",
function(object)
 { cat("An object of class \"de\"\n")
   cat("\nCall:\n", deparse(object@call), "\n\n",sep="")
   cat("Available slots:\n")
   print(slotNames(object))
}) 

summary.de <- function(object, ...)
{
  nvars <- ncol(object@population)
  varnames <- parNames(object)
  domain <- NULL
  domain <- rbind(object@lower, object@upper)
  rownames(domain) <- c("lower", "upper")
  if(ncol(domain) == nvars) 
    colnames(domain) <- varnames
  suggestions <- NULL
  if(nrow(object@suggestions) > 0) 
    { suggestions <- object@suggestions
      dimnames(suggestions) <- list(1:nrow(suggestions), varnames) 
    }
  
  out <- list(type = object@type,
              popSize = object@popSize,
              maxiter = object@maxiter,
              elitism = object@elitism,
              stepsize = if(is.na(object@stepsize)) 
                            "runif(0.5,1.0)" else object@stepsize,
              pcrossover = object@pcrossover,
              pmutation = object@pmutation,
              domain = domain,
              suggestions = suggestions,
              iter = object@iter,
              fitness = object@fitnessValue,
              solution = object@solution)
  class(out) <- "summary.de"
  return(out)
}

setMethod("summary", "de", summary.de)

print.summary.de <- function(x, digits = getOption("digits"), ...)
{
  dotargs <- list(...)
  if(is.null(dotargs$head)) dotargs$head <- 10
  if(is.null(dotargs$tail)) dotargs$tail <- 2
  if(is.null(dotargs$chead)) dotargs$chead <- 10
  if(is.null(dotargs$ctail)) dotargs$ctail <- 2
  
  cat(cli::rule(left = crayon::bold("Differential Evolution"), 
                width = min(getOption("width"),40)), "\n\n")
  cat("DE settings: \n")
  cat(paste("Type                  = ", x$type, "\n"))
  cat(paste("Population size       = ", x$popSize, "\n"))
  cat(paste("Number of generations = ", x$maxiter, "\n"))
  cat(paste("Elitism               = ", x$elitism, "\n"))
  cat(paste("Stepsize              = ", format(x$stepsize, digits = digits), "\n"))
  cat(paste("Crossover probability = ", format(x$pcrossover, digits = digits), "\n"))
  cat(paste("Mutation probability  = ", format(x$pmutation, digits = digits), "\n"))
  #
  cat(paste("Search domain = \n"))
  do.call(".printShortMatrix", 
          c(list(x$domain, digits = digits), 
            dotargs[c("head", "tail", "chead", "ctail")]))
  #
  if(!is.null(x$suggestions))
    { cat(paste("Suggestions =", "\n"))
      do.call(".printShortMatrix", 
              c(list(x$suggestions, digits = digits), 
                dotargs[c("head", "tail", "chead", "ctail")]))
    }
  #
  cat("\nDE results: \n")
  cat(paste("Iterations             =", format(x$iter, digits = digits), "\n"))
  cat(paste("Fitness function value =", format(x$fitness, digits = digits), "\n"))
  if(nrow(x$solution) > 1) 
    { cat(paste("Solutions = \n")) }
  else
    { cat(paste("Solution = \n")) }
  do.call(".printShortMatrix", 
          c(list(x$solution, digits = digits), 
            dotargs[c("head", "tail", "chead", "ctail")]))
  #
  invisible()
}

setMethod("plot", "de", plot.ga)

setMethod("parNames", "de",
function(object, ...)
{ 
  names <- object@names
  nvars <- ncol(object@population)
  if(length(names) == 0)
    { names <- paste("x", 1:nvars, sep = "") }
  return(names)
})


gareal_de <- function(object, F = 0.8, p = 0.5, ...)
{
  if(gaControl("useRcpp"))
    gareal_de_Rcpp(object, fitness = object@call$fitness, F, p)
  else
    gareal_de_R(object, fitness = object@call$fitness, F, p)
}

gareal_de_R <- function(object, fitness,
                        F = 0.8, p = 0.5)
{
# Differential Evolution operator based on the description in Simon (2013)
# Evolutionary Optimization Algorithms, Sec. 12.4, Fig. 12.12
# See also http://mfouesneau.github.io/docs/de/
#  
# object = 'ga' object
# fitness = the fitness function
# p = probability of exchange on [0,1] (crossover probability in DE literature)
# F = stepsize from the interval [0,2]; if NA a random value is selected in
#     [0.5, 1.0] (dithering)

  p <- max(0, min(p, 1))
  F <- max(0, min(F, 2))
  
  pop <- object@population
  f   <- object@fitness
  # fitness <- eval(object@call$fitness) # extract the fitness function
  popSize <- object@popSize
  popseq <- seq_len(popSize)
  n <- ncol(pop)
  nseq <- seq_len(n)
  lb <- object@lower
  ub <- object@upper
  
  for(i in popseq)
  {
    r <- sample(popseq, size = 3, replace = FALSE)
    Fi <- if(is.na(F)) runif(1, 0.5, 1) else F
    v <- pop[r[1],] + Fi*(pop[r[2],] - pop[r[3],])
    J <- sample(nseq, size = 1)
    x <- pop[i,]
    for(j in nseq)
    {
      if(runif(1) < p | j == J) x[j] <- v[j]
      # reset to random amount if outside the bounds
      if(x[j] < lb[j])
        x[j] <- lb[j] + runif(1)*(ub[j] - lb[j])
      if(x[j] > ub[j])
        x[j] <- ub[j] - runif(1)*(ub[j] - lb[j])
    }
    fx <- fitness(x)
    if(fx > f[i])
    {
      f[i] <- fx
      pop[i,] <- x
    }
  }
  out <- list(population = pop, fitness = f)
  return(out)
}

