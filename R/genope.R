##############################################################################
#                                                                            #
#                          GENETIC OPERATORS                                 #
#                                                                            #
##############################################################################

## 
## Generic GA operators
##

ga_lrSelection <- function(object, 
                            r = 2/(object@popSize*(object@popSize-1)), 
                            q = 2/object@popSize, ...)
{
# Linear-rank selection
# Michalewicz (1996) Genetic Algorithms + Data Structures = Evolution Programs. p. 60
  rank <- (object@popSize+1) - rank(object@fitness, ties.method = "min")
  prob <- q - (rank-1)*r
  prob <- pmin(pmax(0, prob/sum(prob)), 1, na.rm = TRUE)
  sel <- sample(1:object@popSize, size = object@popSize, 
                prob = prob, replace = TRUE)
  out <- list(population = object@population[sel,,drop=FALSE],
              fitness = object@fitness[sel])
  return(out)
}

ga_nlrSelection <- function(object, q = 0.25, ...)
{
# Nonlinear-rank selection
# Michalewicz (1996) Genetic Algorithms + Data Structures = Evolution Programs. p. 60
  rank <- (object@popSize + 1) - rank(object@fitness, ties.method = "random")
  prob <- q*(1-q)^(rank-1)
  prob <- pmin(pmax(0, prob/sum(prob)), 1, na.rm = TRUE)
  sel <- sample(1:object@popSize, size = object@popSize, 
                prob = prob, replace = TRUE)
  out <- list(population = object@population[sel,,drop=FALSE],
              fitness = object@fitness[sel])
  return(out)
}

ga_rwSelection <- function(object, ...)
{
# Proportional (roulette wheel) selection
  prob <- abs(object@fitness)/sum(abs(object@fitness))
  prob <- pmin(pmax(0, prob/sum(prob)), 1, na.rm = TRUE)
  sel <- sample(1:object@popSize, size = object@popSize, 
                prob = prob, replace = TRUE)
  out <- list(population = object@population[sel,,drop=FALSE],
              fitness = object@fitness[sel])
  return(out)
}

ga_tourSelection <- function(object, k = 3, ...)
{
# (unbiased) Tournament selection 
  sel <- rep(NA, object@popSize)
  for(i in 1:object@popSize)
     { s <- sample(1:object@popSize, size = k)
       sel[i] <- s[which.max(object@fitness[s])]
     }
  out <- list(population = object@population[sel,,drop=FALSE],
              fitness = object@fitness[sel])
  return(out)
}

ga_spCrossover <- function(object, parents, ...)
{
# Single-point crossover
  fitness <- object@fitness[parents]
  parents <- object@population[parents,,drop = FALSE]
  n <- ncol(parents)
  children <- matrix(as.double(NA), nrow = 2, ncol = n)
  fitnessChildren <- rep(NA, 2)
  crossOverPoint <- sample(0:n, size = 1)
  if(crossOverPoint == 0)
    { children[1:2,] <- parents[2:1,]
      fitnessChildren[1:2] <- fitness[2:1] }
  else if(crossOverPoint == n) 
         { children <- parents
           fitnessChildren <- fitness }
       else 
         { children[1,] <- c(parents[1,1:crossOverPoint],
                             parents[2,(crossOverPoint+1):n])
           children[2,] <- c(parents[2,1:crossOverPoint],
                             parents[1,(crossOverPoint+1):n])
         }
  out <- list(children = children, fitness = fitnessChildren)
  return(out)
}

## 
## Binary GA operators 
##

gabin_Population <- function(object, ...)
{
# Generate a random population of nBits 0/1 values of size popSize
  population <- matrix(as.double(NA), 
                       nrow = object@popSize, 
                       ncol = object@nBits)
  for(j in 1:object@nBits) 
     { population[,j] <- round(runif(object@popSize)) }
  return(population)
}

gabin_lrSelection <- ga_lrSelection

gabin_nlrSelection <- ga_nlrSelection

gabin_rwSelection <- ga_rwSelection

gabin_tourSelection <- ga_tourSelection

gabin_spCrossover <- ga_spCrossover

gabin_uCrossover <- function(object, parents, ...)
{
# Uniform crossover
  parents <- object@population[parents,,drop = FALSE]
  n <- ncol(parents)
  u <- runif(n)
  children <- parents
  children[1:2, u > 0.5] <- children[2:1, u > 0.5]
  out <- list(children = children, fitness = rep(NA,2))  
  return(out)
}

gabin_raMutation <- function(object, parent, ...)
{
# Uniform random mutation
  mutate <- parent <- as.vector(object@population[parent,])
  n <- length(parent)
  j <- sample(1:n, size = 1)
  mutate[j] <- abs(mutate[j]-1)
  return(mutate)
}


## 
## Real-value GA operators 
##

gareal_Population <- function(object, ...)
{
# Generate a random population of size popSize in the range [min, max]  
  min <- object@min
  max <- object@max
  nvars <- length(min)
  population <- matrix(as.double(NA), nrow = object@popSize, ncol = nvars)
  for(j in 1:nvars) 
     { population[,j] <- runif(object@popSize, min[j], max[j]) }
  return(population)
}

gareal_lrSelection <- ga_lrSelection

gareal_nlrSelection <- ga_nlrSelection

gareal_rwSelection <- ga_rwSelection

gareal_tourSelection <- ga_tourSelection

gareal_lsSelection <- function(object, ...)
{
# Fitness proportional selection with fitness linear scaling  
  popSize <- object@popSize
  f <- object@fitness
  fmin <- min(f, na.rm = TRUE)
  if(fmin < 0) 
    { f <- f - fmin
      fmin <- min(f, na.rm = TRUE) }
  fave <- mean(f, na.rm = TRUE)
  fmax <- max(f, na.rm = TRUE)
  sfactor <- 2 # scaling factor
  # transform f -> f' = a*f + b such that
  if(fmin > (sfactor*fave - fmax)/(sfactor-1))
    { # ave(f) = ave(f')
      # 2*ave(f') = max(f')
      delta <- fmax - fave
      a <- (sfactor - 1.0)*fave/delta
      b <- fave * (fmax - sfactor*fave)/delta 
    }
  else
    { # ave(f) = ave(f')
      # min(f') = 0
      delta <- fave - fmin
      a <- fave/delta
      b <- -1*fmin*fave/delta 
    }
  fscaled <- a*f + b
  prob <- abs(fscaled)/sum(abs(fscaled), na.rm = TRUE)
  prob <- pmin(pmax(0, prob/sum(prob)), 1, na.rm = TRUE)
  sel <- sample(1:object@popSize, size = object@popSize, 
                prob = prob, replace = TRUE)
  out <- list(population = object@population[sel,,drop=FALSE],
              fitness = object@fitness[sel])
  return(out)
}

gareal_sigmaSelection <- function(object, ...)
{
# Fitness proportional selection with Goldberg's Sigma Truncation Scaling
  popSize <- object@popSize
  mf <- mean(object@fitness, na.rm = TRUE)
  sf <- sd(object@fitness, na.rm = TRUE)
  fscaled <- pmax(object@fitness - (mf - 2*sf), 0, na.rm = TRUE)
  prob <- abs(fscaled)/sum(abs(fscaled))
  prob <- pmin(pmax(0, prob/sum(prob)), 1, na.rm = TRUE)
  sel <- sample(1:object@popSize, size = object@popSize,
                prob = prob, replace = TRUE)
  out <- list(population = object@population[sel,,drop=FALSE],
              fitness = object@fitness[sel])
  return(out)
}

gareal_spCrossover <- ga_spCrossover

gareal_waCrossover <- function(object, parents, ...)
{
# Whole arithmetic crossover
  parents <- object@population[parents,,drop = FALSE]
  n <- ncol(parents)
  children <- matrix(as.double(NA), nrow = 2, ncol = n)
  a <- runif(1)
  children[1,] <- a*parents[1,] + (1-a)*parents[2,]
  children[2,] <- a*parents[2,] + (1-a)*parents[1,]
  out <- list(children = children, fitness = rep(NA,2))
  return(out)
}

gareal_laCrossover <- function(object, parents, ...)
{
# Local arithmetic crossover
  parents <- object@population[parents,,drop = FALSE]
  n <- ncol(parents)
  children <- matrix(as.double(NA), nrow = 2, ncol = n)
  a <- runif(n)
  children[1,] <- a*parents[1,] + (1-a)*parents[2,]
  children[2,] <- a*parents[2,] + (1-a)*parents[1,]
  out <- list(children = children, fitness = rep(as.double(NA),2))
  return(out)
}

# gareal_blxCrossover <- function(object, parents, ...)
# {
# # Blend crossover
#   parents <- object@population[parents,,drop = FALSE]
#   n <- ncol(parents)
#   a <- 0.5
#   # a <- exp(-pi*iter/max(iter)) # annealing factor
#   children <- matrix(as.double(NA), nrow = 2, ncol = n)
#   for(i in 1:n)
#      { x <- sort(parents[,i])
#        xl <- max(x[1] - a*(x[2]-x[1]), object@min[i])
#        xu <- min(x[2] + a*(x[2]-x[1]), object@max[i])
#        children[,i] <- runif(2, xl, xu) 
#      }
#   out <- list(children = children, fitness = rep(NA,2))
#   return(out)
# }

# versione ottimizzata 
gareal_blxCrossover <- function(object, parents, ...)
{
# Blend crossover
  parents <- object@population[parents,,drop = FALSE]
  n <- ncol(parents)
  a <- 0.5
  x <- apply(parents, 2, range)
  xl <- pmax(x[1,] - a*(x[2,]-x[1,]), object@min)
  xu <- pmin(x[2,] + a*(x[2,]-x[1,]), object@max)
  children <- matrix(as.double(NA), nrow = 2, ncol = n)
  for(i in 1:n)
    children[,i] <- runif(2, xl[i], xu[i]) 
  out <- list(children = children, fitness = rep(NA,2))
  return(out)
}


gareal_laplaceCrossover <- function (object, parents, a = 0, b = 0.15, ...) 
{
# Laplace crossover(a, b)
#
# a is the location parameter and b > 0 is the scaling parameter of a Laplace
# distribution, which is generated as described in 
# Krishnamoorthy K. (2006) Handbook of Statistical Distributions with 
#   Applications, Chapman & Hall/CRC.
#
# For smaller values of b offsprings are likely to be produced nearer to 
# parents, and for larger values of b offsprings are expected to be produced
# far from parents.
# Deep et al. (2009) suggests to use a = 0, b = 0.15 for real-valued 
# variables, and b = 0.35 for integer variables.
#
# References
#
# Deep K., Thakur M. (2007) A new crossover operator for real coded genetic
#   algorithms, Applied Mathematics and Computation, 188, 895–912.
# Deep K., Singh K.P., Kansal M.L., Mohan C. (2009) A real coded genetic
#   algorithm for solving integer and mixed integer optimization problems.
#   Applied Mathematics and Computation, 212(2), pp. 505-518.
  
  parents <- object@population[parents, , drop = FALSE]
  n <- ncol(parents)
  children <- matrix(as.double(NA), nrow = 2, ncol = n)
  r <- runif(n)
  u <- runif(n)
  beta <- a + ifelse(r > 0.5, b*log(u), -b*log(u))
  bpar <- beta*abs(parents[1,] - parents[2,])
  children[1,] <- pmin(pmax(parents[1,] + bpar, object@min), object@max)
  children[2,] <- pmin(pmax(parents[2,] + bpar, object@min), object@max)
  out <- list(children = children, fitness = rep(NA, 2))
  return(out)
}


gareal_raMutation <- function(object, parent, ...)
{
# Uniform random mutation
  mutate <- parent <- as.vector(object@population[parent,])
  n <- length(parent)
  j <- sample(1:n, size = 1)
  mutate[j] <- runif(1, object@min[j], object@max[j])
  return(mutate)
}

gareal_nraMutation <- function(object, parent, ...)
{
# Non uniform random mutation
  mutate <- parent <- as.vector(object@population[parent,])
  n <- length(parent)
  g <- 1 - object@iter/object@maxiter # dempening factor
  sa <- function(x) x*(1-runif(1)^g)
  j <- sample(1:n, 1)
  u <- runif(1)
  if(u < 0.5)
    { mutate[j] <- parent[j] - sa(parent[j] - object@max[j]) }
  else
    { mutate[j] <- parent[j] + sa(object@max[j] - parent[j]) }
  return(mutate)
}

gareal_rsMutation <- function(object, parent, ...)
{
# Random mutation around the solution
  mutate <- parent <- as.vector(object@population[parent,])
  dempeningFactor <- 1 - object@iter/object@maxiter
  direction <- sample(c(-1,1),1)
  value <- (object@max - object@min)*0.67
  mutate <- parent + direction*value*dempeningFactor
  outside <- (mutate < object@min | mutate > object@max)
  for(j in which(outside))
     { mutate[j] <- runif(1, object@min[j], object@max[j]) }
  return(mutate)
}

gareal_powMutation <- function(object, parent, pow = 10, ...)
{
# Power mutation(pow)
#
# Deep et al. (2009) suggests to use pow = 10 for real-valued variables, and
# pow = 4 for integer variables.
#
# References
#
# Deep K., Singh K.P., Kansal M.L., Mohan C. (2009) A real coded genetic
#   algorithm for solving integer and mixed integer optimization problems.
#   Applied Mathematics and Computation, 212(2), pp. 505-518.
# Deep K., Thakur M. (2007) A new mutation operator for real coded genetic
#  algorithms, Applied Mathematics and Computation, 193, pp. 211–230.

  mutate <- parent <- as.vector(object@population[parent,])
  n <- length(parent)
  s <- runif(1)^pow
  t <- (parent - object@min)/(object@max - parent)
  r <- runif(n)
  mutate <- parent + ifelse(r < t, 
                            -s*(parent - object@min),
                            +s*(object@max - parent))
  return(mutate)
}

## 
## Permutation GA operators 
##

gaperm_Population <- function(object, ...)
{
# Generate popSize random permutations in the range [min, max]
  int <- seq.int(object@min, object@max)
  n <- length(int)
  population <- matrix(NA, nrow = object@popSize, ncol = n)
  for(i in 1:object@popSize)
     population[i,] <- sample(int, replace = FALSE)
  return(population)
}

gaperm_lrSelection <- ga_lrSelection

gaperm_nlrSelection <- ga_nlrSelection

gaperm_rwSelection <- ga_rwSelection

gaperm_tourSelection <- ga_tourSelection

gaperm_cxCrossover <- function(object, parents, ...)
{
# Cycle crossover (CX)
  parents <- object@population[parents,,drop = FALSE]
  n <- ncol(parents)
  cxPoint <- 1 # sample(1:n, size = 1)
  children <- parents
  children[1:2,cxPoint] <- parents[2:1,cxPoint]
  while( length(dup <- which(duplicated(children[1,], fromLast = TRUE))) > 0 )
       { children[1:2,dup] <- children[2:1,dup] }
  out <- list(children = children, fitness = rep(NA,2))
  return(out)
}

gaperm_pmxCrossover <- function(object, parents, ...)
{
# Partially matched crossover (PMX)
  parents <- object@population[parents,,drop = FALSE]
  n <- ncol(parents)
  cxPoints <- sample(1:n, size = 2)
  cxPoints <- seq(min(cxPoints), max(cxPoints))
  children <- matrix(as.double(NA), nrow = 2, ncol = n)
  children[,cxPoints] <- parents[,cxPoints]
  for(i in setdiff(1:n, cxPoints))
     { if(!any(parents[2,i] == children[1,cxPoints]))
         { children[1,i] <- parents[2,i] }
       if(!any(parents[1,i] == children[2,cxPoints]))
         { children[2,i] <- parents[1,i] }
     }
  children[1,is.na(children[1,])] <- setdiff(parents[2,], children[1,])
  children[2,is.na(children[2,])] <- setdiff(parents[1,], children[2,])
  out <- list(children = children, fitness = rep(NA,2))
  return(out)
}

gaperm_oxCrossover <- function(object, parents, ...)
{
# Order Crossover (OX)
  parents <- object@population[parents,,drop = FALSE]
  n <- ncol(parents)
  #
  cxPoints <- sample(seq(2,n-1), size = 2)
  cxPoints <- seq(min(cxPoints), max(cxPoints))
  children <- matrix(as.double(NA), nrow = 2, ncol = n)
  children[,cxPoints] <- parents[,cxPoints]
  #
  for(j in 1:2)
     { pos <- c((max(cxPoints)+1):n, 1:(max(cxPoints)))
       val <- setdiff(parents[-j,pos], children[j,cxPoints])
       i <- intersect(pos, which(is.na(children[j,])))
       children[j,i] <- val
     }
  #
  out <- list(children = children, fitness = rep(NA,2))
  return(out)
}

gaperm_pbxCrossover <- function(object, parents, ...)
{
# Position-based crossover (PBX)
  parents <- object@population[parents,,drop = FALSE]
  n <- ncol(parents)
  #
  cxPoints <- unique(sample(1:n, size = n, replace = TRUE))
  children <- matrix(as.double(NA), nrow = 2, ncol = n)
  children[1,cxPoints] <- parents[2,cxPoints]
  children[2,cxPoints] <- parents[1,cxPoints]
  #
  for(j in 1:2)
     { pos <- which(is.na(children[j,]))
       val <- setdiff(parents[-j,], children[j,cxPoints])
       children[j,pos] <- val
     }
  #
  out <- list(children = children, fitness = rep(NA,2))
  return(out)
}

gaperm_simMutation <- function(object, parent, ...)
{
# Simple inversion mutation
  parent <- as.vector(object@population[parent,])
  n <- length(parent)
  m <- sort(sample(1:n, size = 2))
  m <- seq(m[1], m[2], by = 1)
  if(min(m)==1 & max(m)==n)
    i <- rev(m)
  else if(min(m)==1) 
       i <- c(rev(m), seq(max(m)+1, n, by = 1))
       else if(max(m)==n) 
            i <- c(seq(1, min(m)-1, by = 1), rev(m))
            else i <- c(seq(1, min(m)-1, by = 1), rev(m), seq(max(m)+1, n, by = 1))
  mutate <- parent[i]
  return(mutate)
} 

gaperm_ismMutation <- function(object, parent, ...)
{
# Insertion mutation
  parent <- as.vector(object@population[parent,])
  n <- length(parent)
  m <- sample(1:n, size = 1)
  pos <- sample(1:(n-1), size = 1)
  i <- c(setdiff(1:pos,m), m, setdiff((pos+1):n,m))
  mutate <- parent[i]
  return(mutate)
}

gaperm_swMutation <- function(object, parent, ...)
{
# Exchange mutation or swap mutation
  mutate <- parent <- as.vector(object@population[parent,])
  n <- length(parent)
  m <- sample(1:n, size = 2)
  mutate[m[1]] <- parent[m[2]]
  mutate[m[2]] <- parent[m[1]]
  return(mutate)
}

gaperm_dmMutation <- function(object, parent, ...)
{
# Displacement mutation
  parent <- as.vector(object@population[parent,])
  n <- length(parent)
  m <- sort(sample(1:n, size = 2))
  m <- seq(m[1], m[2], by = 1)
  l <- max(m)-min(m)+1
  pos <- sample(1:max(1,(n-l)), size = 1)
  i <- c(setdiff(1:n,m)[1:pos], m, setdiff(1:n,m)[-(1:pos)])
  mutate <- parent[na.omit(i)]
  return(mutate)
} 

gaperm_scrMutation <- function(object, parent, ...)
{
# Scramble mutation
  parent <- as.vector(object@population[parent,])
  n <- length(parent)
  m <- sort(sample(1:n, size = 2))
  m <- seq(min(m), max(m), by = 1)
  m <- sample(m, replace = FALSE)
  i <- c(setdiff(1:min(m),m), m, setdiff(max(m):n,m))
  mutate <- parent[i]
  return(mutate)
} 

ga_pmutation <- function(object, p0 = 0.5, p = 0.01, 
                         T = round(object@maxiter/2), ...)
{
# variable probability of mutation
# p0 = initial pmutation
# p = final pmutation
# T = maximum iteration after which converges to p
#
# Example:
# p0 = 0.5; p = 0.01; 
# maxiter = 1000; T = round(maxiter/2); t = seq(maxiter)
# pm1 = ifelse(t > T, p, p0 - (p0-p)/T * (t-1)) # linear decay
# pm2 = (p0 - p)*exp(-2*(t-1)/T) + p # exponential decay
# plot(t, pm1, type = "l")
# lines(t, pm2, col = 2)

  t <- object@iter
  # linear decay
  # pm <- if(t > T) p else p0 - (p0-p)/T * (t-1)
  # exponential decay
  pm = (p0 - p)*exp(-2*(t-1)/T) + p
  #
  return(pm)
}

# Probability of selection based on fitness values in vector x with
# selection pressure given by q in (0,1).
# This is used in optim() local search to select which solution should
# be used for starting the algorithm.
optimProbsel <- function(x, q = 0.25)
{
  x <- as.vector(x)
  n <- length(x)
  # selection pressure parameter
  q <- min(max(sqrt(.Machine$double.eps), q), 
           1 - sqrt(.Machine$double.eps))
  rank <- (n + 1) - rank(x, ties.method = "first", na.last = FALSE)
  # prob <- q*(1-q)^(rank-1) * 1/(1-(1-q)^n)
  prob <- q*(1-q)^(rank-1)
  prob[is.na(x)] <- 0
  prob <- pmin(pmax(0, prob/sum(prob)), 1, na.rm = TRUE)
  return(prob)
}
  