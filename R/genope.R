##--------------------------------------------------------------------##
##                                                                    ##
##                        GENETIC OPERATORS                           ##
##                                                                    ##
##--------------------------------------------------------------------##

##
## Generic GA operators ----
##

# Linear-rank selection ----
# Michalewicz (1996) Genetic Algorithms + Data Structures = Evolution Programs. p. 60

ga_lrSelection <- function(object, 
                           r = 2/(object@popSize*(object@popSize-1)), 
                           q = 2/object@popSize, ...)
{
  if(gaControl("useRcpp"))
    ga_lrSelection_Rcpp(object, r, q)
  else
    ga_lrSelection_R(object, r, q)
}

ga_lrSelection_R <- function(object, r, q)
{
  if(missing(r)) r <- 2/(object@popSize * (object@popSize - 1))
  if(missing(q)) q <- 2/object@popSize
  rank <- (object@popSize+1) - rank(object@fitness, ties.method = "min")
  prob <- 1 + q - (rank-1)*r
  prob <- pmin(pmax(0, prob/sum(prob)), 1, na.rm = TRUE)
  sel <- sample(1:object@popSize, size = object@popSize, 
                prob = prob, replace = TRUE)
  out <- list(population = object@population[sel,,drop=FALSE],
              fitness = object@fitness[sel])
  return(out)
}

# Nonlinear-rank selection ----
# Michalewicz (1996) Genetic Algorithms + Data Structures = Evolution Programs. p. 60

ga_nlrSelection <- function(object, q = 0.25, ...)
{
  if(gaControl("useRcpp"))
    ga_nlrSelection_Rcpp(object, q)
  else
    ga_nlrSelection_R(object, q)
}

ga_nlrSelection_R <- function(object, q)
{
  if(missing(q)) q <- 0.25
  rank <- (object@popSize + 1) - rank(object@fitness, ties.method = "min")
  prob <- q*(1-q)^(rank-1)
  prob <- pmin(pmax(0, prob/sum(prob)), 1, na.rm = TRUE)
  sel <- sample(1:object@popSize, size = object@popSize, 
                prob = prob, replace = TRUE)
  out <- list(population = object@population[sel,,drop=FALSE],
              fitness = object@fitness[sel])
  return(out)
}

# Proportional (roulette wheel) selection ----

ga_rwSelection <- function(object, ...)
{
  if(gaControl("useRcpp"))
    ga_rwSelection_Rcpp(object)
  else
    ga_rwSelection_R(object)
}

ga_rwSelection_R <- function(object, ...)
{
  prob <- abs(object@fitness)/sum(abs(object@fitness))
  prob <- pmin(pmax(0, prob/sum(prob)), 1, na.rm = TRUE)
  sel <- sample(1:object@popSize, size = object@popSize, 
                prob = prob, replace = TRUE)
  out <- list(population = object@population[sel,,drop=FALSE],
              fitness = object@fitness[sel])
  return(out)
}

# (unbiased) Tournament selection ----

ga_tourSelection <- function(object, k = 3, ...)
{
  if(gaControl("useRcpp"))
    ga_tourSelection_Rcpp(object, k)
  else
    ga_tourSelection_R(object, k)
}

ga_tourSelection_R <- function(object, k = 3, ...)
{
  sel <- rep(NA, object@popSize)
  for(i in 1:object@popSize)
     { s <- sample(1:object@popSize, size = k)
       sel[i] <- s[which.max(object@fitness[s])]
     }
  out <- list(population = object@population[sel,,drop=FALSE],
              fitness = object@fitness[sel])
  return(out)
}

# Single-point crossover ----

ga_spCrossover <- function(object, parents, ...)
{
  if(gaControl("useRcpp"))
    ga_spCrossover_Rcpp(object, parents)
  else
    ga_spCrossover_R(object, parents)
}

ga_spCrossover_R <- function(object, parents)
{
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
## Binary GA operators  ----
##

# Generate a binary random population ----

gabin_Population <- function(object, ...)
{
  if(gaControl("useRcpp"))
    gabin_Population_Rcpp(object)
  else
    gabin_Population_R(object)
}

gabin_Population_R <- function(object)
{
  population <- matrix(as.double(NA), 
                       nrow = object@popSize, 
                       ncol = object@nBits)
  for(j in 1:object@nBits) 
     { population[,j] <- round(runif(object@popSize)) }
  storage.mode(population) <- "integer"
  return(population)
}

gabin_lrSelection      <- ga_lrSelection
gabin_lrSelection_R    <- ga_lrSelection_R
gabin_lrSelection_Rcpp <- ga_lrSelection_Rcpp

gabin_nlrSelection      <- ga_nlrSelection
gabin_nlrSelection_R    <- ga_nlrSelection_R
gabin_nlrSelection_Rcpp <- ga_nlrSelection_Rcpp

gabin_rwSelection      <- ga_rwSelection
gabin_rwSelection_R    <- ga_rwSelection_R
gabin_rwSelection_Rcpp <- ga_rwSelection_Rcpp

gabin_tourSelection      <- ga_tourSelection
gabin_tourSelection_R    <- ga_tourSelection_R
gabin_tourSelection_Rcpp <- ga_tourSelection_Rcpp

gabin_spCrossover      <- ga_spCrossover
gabin_spCrossover_R    <- ga_spCrossover_R
gabin_spCrossover_Rcpp <- ga_spCrossover_Rcpp

# Uniform crossover ----

gabin_uCrossover <- function(object, parents, ...)
{
  if(gaControl("useRcpp"))
    gabin_uCrossover_Rcpp(object, parents)
  else
    gabin_uCrossover_R(object, parents)
}

gabin_uCrossover_R <- function(object, parents)
{
  parents <- object@population[parents,,drop = FALSE]
  n <- ncol(parents)
  u <- runif(n)
  children <- parents
  children[1:2, u > 0.5] <- children[2:1, u > 0.5]
  out <- list(children = children, fitness = rep(NA,2))  
  return(out)
}

# Uniform random mutation ----

gabin_raMutation <- function(object, parent, ...)
{
  if(gaControl("useRcpp"))
    gabin_raMutation_Rcpp(object, parent)
  else
    gabin_raMutation_R(object, parent)
}

gabin_raMutation_R <- function(object, parent)
{
  mutate <- parent <- as.vector(object@population[parent,])
  n <- length(parent)
  j <- sample(1:n, size = 1)
  mutate[j] <- abs(mutate[j]-1)
  return(mutate)
}

## 
## Real-value GA operators ----
##

# Generate a random population ----

gareal_Population <- function(object, ...)
{
  if(gaControl("useRcpp"))
    gareal_Population_Rcpp(object)
  else
    gareal_Population_R(object)
}

gareal_Population_R <- function(object)
{
  lower <- object@lower
  upper <- object@upper
  nvars <- length(lower)
  population <- matrix(as.double(NA), nrow = object@popSize, ncol = nvars)
  for(j in 1:nvars) 
     { population[,j] <- runif(object@popSize, lower[j], upper[j]) }
  return(population)
}

gareal_lrSelection      <- ga_lrSelection
gareal_lrSelection_R    <- ga_lrSelection_R
gareal_lrSelection_Rcpp <- ga_lrSelection_Rcpp

gareal_nlrSelection      <- ga_nlrSelection
gareal_nlrSelection_R    <- ga_nlrSelection_R
gareal_nlrSelection_Rcpp <- ga_nlrSelection_Rcpp

gareal_rwSelection      <- ga_rwSelection
gareal_rwSelection_R    <- ga_rwSelection_R
gareal_rwSelection_Rcpp <- ga_rwSelection_Rcpp

gareal_tourSelection      <- ga_tourSelection
gareal_tourSelection_R    <- ga_tourSelection_R
gareal_tourSelection_Rcpp <- ga_tourSelection_Rcpp

# Fitness proportional selection with fitness linear scaling ----

gareal_lsSelection <- function(object, ...)
{
  if(gaControl("useRcpp"))
    gareal_lsSelection_Rcpp(object)
  else
    gareal_lsSelection_R(object)
}

gareal_lsSelection_R <- function(object)
{
  f <- object@fitness
  fmin <- min(f, na.rm = TRUE)
  if(fmin < 0) 
    { f <- f - fmin
      fmin <- min(f, na.rm = TRUE) }
  fave <- mean(f, na.rm = TRUE)
  fmax <- max(f, na.rm = TRUE)
  sfactor <- 2 # scaling factor
  eps <- sqrt(.Machine$double.eps)
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
  prob[is.na(prob)] <- eps
  prob <- pmin(pmax(0.0, prob/sum(prob)), 1.0)
  sel <- sample(1:object@popSize, size = object@popSize, 
                prob = prob, replace = TRUE)
  out <- list(population = object@population[sel,,drop=FALSE],
              fitness = object@fitness[sel])
  return(out)
}


# Fitness proportional selection with Goldberg's Sigma Truncation Scaling ----

gareal_sigmaSelection <- function(object, ...)
{
  if(gaControl("useRcpp"))
    gareal_sigmaSelection_Rcpp(object)
  else
    gareal_sigmaSelection_R(object)
}

gareal_sigmaSelection_R <- function(object)
{
  popSize <- object@popSize
  mf <- mean(object@fitness, na.rm = TRUE)
  sf <- sd(object@fitness, na.rm = TRUE)
  fscaled <- pmax(object@fitness - (mf - 2*sf), 0, na.rm = TRUE)
  prob <- abs(fscaled)/sum(abs(fscaled))
  prob <- pmin(pmax(0, prob, na.rm = TRUE), 1, na.rm = TRUE)
  sel <- sample(1:popSize, size = popSize, prob = prob, replace = TRUE)
  out <- list(population = object@population[sel,,drop=FALSE],
              fitness = object@fitness[sel])
  return(out)
}

gareal_spCrossover      <- ga_spCrossover
gareal_spCrossover_R    <- ga_spCrossover_R
gareal_spCrossover_Rcpp <- ga_spCrossover_Rcpp

# Whole arithmetic crossover ----

gareal_waCrossover <- function(object, parents, ...)
{
  if(gaControl("useRcpp"))
    gareal_waCrossover_Rcpp(object, parents)
  else
    gareal_waCrossover_R(object, parents)
}

gareal_waCrossover_R <- function(object, parents)
{
  parents <- object@population[parents,,drop = FALSE]
  n <- ncol(parents)
  children <- matrix(as.double(NA), nrow = 2, ncol = n)
  a <- runif(1)
  children[1,] <- a*parents[1,] + (1-a)*parents[2,]
  children[2,] <- a*parents[2,] + (1-a)*parents[1,]
  out <- list(children = children, fitness = rep(NA,2))
  return(out)
}

# Local arithmetic crossover

gareal_laCrossover <- function(object, parents, ...)
{
  if(gaControl("useRcpp"))
    gareal_laCrossover_Rcpp(object, parents)
  else
    gareal_laCrossover_R(object, parents)
}

gareal_laCrossover_R <- function(object, parents)
{
  parents <- object@population[parents,,drop = FALSE]
  n <- ncol(parents)
  children <- matrix(as.double(NA), nrow = 2, ncol = n)
  a <- runif(n)
  children[1,] <- a*parents[1,] + (1-a)*parents[2,]
  children[2,] <- a*parents[2,] + (1-a)*parents[1,]
  out <- list(children = children, fitness = rep(as.double(NA),2))
  return(out)
}

# Blend crossover ----

gareal_blxCrossover <- function(object, parents, a = 0.5, ...)
{
  if(gaControl("useRcpp"))
    gareal_blxCrossover_Rcpp(object, parents, a)
  else
    gareal_blxCrossover_R(object, parents, a)
}

gareal_blxCrossover_R <- function(object, parents, a)
{
  if(missing(a)) a <- 0.5
  parents <- object@population[parents,,drop = FALSE]
  n <- ncol(parents)
  x <- apply(parents, 2, range)
  xl <- pmax(x[1,] - a*(x[2,]-x[1,]), object@lower)
  xu <- pmin(x[2,] + a*(x[2,]-x[1,]), object@upper)
  children <- matrix(as.double(NA), nrow = 2, ncol = n)
  for(i in 1:n)
    children[,i] <- runif(2, xl[i], xu[i]) 
  out <- list(children = children, fitness = rep(NA,2))
  return(out)
}

# Laplace crossover ----
#
# Laplace crossover(a, b), where a is the location parameter and b > 0 is 
# the scaling parameter of a Laplace distribution, which is generated as 
# described in 
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
#   algorithms, Applied Mathematics and Computation, 188, 895-912.
# Deep K., Singh K.P., Kansal M.L., Mohan C. (2009) A real coded genetic
#   algorithm for solving integer and mixed integer optimization problems.
#   Applied Mathematics and Computation, 212(2), pp. 505-518.

gareal_laplaceCrossover <- function(object, parents, a = 0, b = 0.15, ...) 
{
  if(gaControl("useRcpp"))
    gareal_laplaceCrossover_Rcpp(object, parents, a, b)
  else
    gareal_laplaceCrossover_R(object, parents, a, b)
}

gareal_laplaceCrossover_R <- function(object, parents, a, b) 
{
  if(missing(a)) a <- 0.00
  if(missing(b)) b <- 0.15
  parents <- object@population[parents,,drop = FALSE]
  n <- ncol(parents)
  if(length(a) == 1) a <- rep(a, n)
  if(length(b) == 1) b <- rep(b, n)
  children <- matrix(as.double(NA), nrow = 2, ncol = n)
  r <- runif(n)
  u <- runif(n)
  beta <- a + ifelse(r > 0.5, b*log(u), -b*log(u))
  bpar <- beta*abs(parents[1,] - parents[2,])
  children[1,] <- pmin(pmax(parents[1,] + bpar, object@lower), object@upper)
  children[2,] <- pmin(pmax(parents[2,] + bpar, object@lower), object@upper)
  out <- list(children = children, fitness = rep(NA, 2))
  return(out)
}

# Uniform random mutation ----

gareal_raMutation <- function(object, parent, ...)
{
  if(gaControl("useRcpp"))
    gareal_raMutation_Rcpp(object, parent)
  else
    gareal_raMutation_R(object, parent)
}

gareal_raMutation_R <- function(object, parent)
{
  mutate <- parent <- as.vector(object@population[parent,])
  n <- length(parent)
  j <- sample(1:n, size = 1)
  mutate[j] <- runif(1, object@lower[j], object@upper[j])
  return(mutate)
}

# Non uniform random mutation ----

gareal_nraMutation <- function(object, parent, ...)
{
  if(gaControl("useRcpp"))
    gareal_nraMutation_Rcpp(object, parent)
  else
    gareal_nraMutation_R(object, parent)
}

gareal_nraMutation_R <- function(object, parent, ...)
{
  mutate <- parent <- as.vector(object@population[parent,])
  n <- length(parent)
  g <- 1 - object@iter/object@maxiter # dempening factor
  sa <- function(x) x*(1-runif(1)^g)
  j <- sample(1:n, 1)
  u <- runif(1)
  if(u < 0.5)
    { mutate[j] <- parent[j] - sa(parent[j] - object@lower[j]) }
  else
    { mutate[j] <- parent[j] + sa(object@upper[j] - parent[j]) }
  return(mutate)
}

# Random mutation around the solution ----

gareal_rsMutation <- function(object, parent, ...)
{
  if(gaControl("useRcpp"))
    gareal_rsMutation_Rcpp(object, parent)
  else
    gareal_rsMutation_R(object, parent)
}

gareal_rsMutation_R <- function(object, parent)
{
  mutate <- parent <- as.vector(object@population[parent,])
  dempeningFactor <- 1 - object@iter/object@maxiter
  direction <- sample(c(-1,1),1)
  value <- (object@upper - object@lower)*0.67
  mutate <- parent + direction*value*dempeningFactor
  outside <- (mutate < object@lower | mutate > object@upper)
  for(j in which(outside))
     { mutate[j] <- runif(1, object@lower[j], object@upper[j]) }
  return(mutate)
}

# Power mutation ----
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
#  algorithms, Applied Mathematics and Computation, 193, pp. 211-230.

gareal_powMutation <- function(object, parent, pow = 10, ...)
{
  if(gaControl("useRcpp"))
    gareal_powMutation_Rcpp(object, parent, pow)
  else
    gareal_powMutation_R(object, parent, pow)
}

gareal_powMutation_R <- function(object, parent, pow)
{
  if(missing(pow)) pow <- 10
  mutate <- parent <- as.vector(object@population[parent,])
  n <- length(parent)
  if(length(pow) == 1) pow <- rep(pow, n)
  s <- runif(1)^pow
  t <- (parent - object@lower)/(object@upper - parent)
  r <- runif(n)
  mutate <- parent + ifelse(r < t, 
                            -s*(parent - object@lower),
                            +s*(object@upper - parent))
  return(mutate)
}

## 
## Permutation GA operators ----
##

# Generate a permutation random population ----

gaperm_Population <- function(object, ...)
{
  if(gaControl("useRcpp"))
    gaperm_Population_Rcpp(object)
  else
    gaperm_Population_R(object)
}

gaperm_Population_R <- function(object)
{
  int <- seq.int(object@lower, object@upper)
  n <- length(int)
  population <- matrix(NA, nrow = object@popSize, ncol = n)
  for(i in 1:object@popSize)
     population[i,] <- sample(int, replace = FALSE)
  return(population)
}

gaperm_lrSelection      <- ga_lrSelection
gaperm_lrSelection_R    <- ga_lrSelection_R
gaperm_lrSelection_Rcpp <- ga_lrSelection_Rcpp

gaperm_nlrSelection      <- ga_nlrSelection
gaperm_nlrSelection_R    <- ga_nlrSelection_R
gaperm_nlrSelection_Rcpp <- ga_nlrSelection_Rcpp

gaperm_rwSelection      <- ga_rwSelection
gaperm_rwSelection_R    <- ga_rwSelection_R
gaperm_rwSelection_Rcpp <- ga_rwSelection_Rcpp

gaperm_tourSelection      <- ga_tourSelection
gaperm_tourSelection_R    <- ga_tourSelection_R
gaperm_tourSelection_Rcpp <- ga_tourSelection_Rcpp

# Cycle crossover (CX) ----

gaperm_cxCrossover <- function(object, parents, ...)
{
  if(gaControl("useRcpp"))
    gaperm_cxCrossover_Rcpp(object, parents)
  else
    gaperm_cxCrossover_R(object, parents)
}

gaperm_cxCrossover_R <- function(object, parents)
{
  parents <- object@population[parents,,drop = FALSE]
  n <- ncol(parents)
  children <- matrix(NA_integer_, nrow = 2, ncol = n) 
  k <- 1 # cx point 
  ALL <- 1:n
  while(length(ALL) > 0)
  {
    i <- ALL[1]
    # perform a cycle
    base <- parents[1,i]
    vi <- parents[2,i]
    I <- i
    while(vi != base)
    {
      i <- which(parents[1,] == parents[2,i])
      vi <- parents[2,i]
      I <- c(I,i)
    }
    ALL = setdiff(ALL,I)
    if(k %%2 == 1) 
      children[,I] <- parents[,I] 
    else 
      children[,I] <- parents[2:1,I]
    k <- k+1 
  } 
  out <- list(children = children, fitness = rep(NA,2))
  return(out)
}

# Partially matched crossover (PMX) ----

gaperm_pmxCrossover <- function(object, parents, ...)
{
  if(gaControl("useRcpp"))
    gaperm_pmxCrossover_Rcpp(object, parents)
  else
    gaperm_pmxCrossover_R(object, parents)
}

gaperm_pmxCrossover_R <- function(object, parents)
{
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

# Order Crossover (OX) ----

gaperm_oxCrossover <- function(object, parents, ...)
{
  if(gaControl("useRcpp"))
    gaperm_oxCrossover_Rcpp(object, parents)
  else
    gaperm_oxCrossover_R(object, parents)
}

gaperm_oxCrossover_R <- function(object, parents)
{
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
       ival <- intersect(pos, which(is.na(children[j,])))
       children[j,ival] <- val
     }
  #
  out <- list(children = children, fitness = rep(NA,2))
  return(out)
}

# Position-based crossover (PBX) ----

gaperm_pbxCrossover <- function(object, parents, ...)
{
  if(gaControl("useRcpp"))
    gaperm_pbxCrossover_Rcpp(object, parents)
  else
    gaperm_pbxCrossover_R(object, parents)
}

gaperm_pbxCrossover_R <- function(object, parents)
{
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
       val <- setdiff(parents[j,], children[j,cxPoints])
       children[j,pos] <- val
     }
  #
  out <- list(children = children, fitness = rep(NA,2))
  return(out)
}

# Simple inversion mutation ----

gaperm_simMutation <- function(object, parent, ...)
{
  if(gaControl("useRcpp"))
    gaperm_simMutation_Rcpp(object, parent)
  else
    gaperm_simMutation_R(object, parent)
}

gaperm_simMutation_R <- function(object, parent)
{
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

# Insertion mutation ----

gaperm_ismMutation <- function(object, parent, ...)
{
  if(gaControl("useRcpp"))
    gaperm_ismMutation_Rcpp(object, parent)
  else
    gaperm_ismMutation_R(object, parent)
}

gaperm_ismMutation_R <- function(object, parent)
{
  parent <- as.vector(object@population[parent,])
  n <- length(parent)
  m <- sample(1:n, size = 1)
  pos <- sample(1:(n-1), size = 1)
  i <- c(setdiff(1:pos,m), m, setdiff((pos+1):n,m))
  mutate <- parent[i]
  return(mutate)
}

# Exchange mutation or swap mutation ----

gaperm_swMutation <- function(object, parent, ...)
{
  if(gaControl("useRcpp"))
    gaperm_swMutation_Rcpp(object, parent)
  else
    gaperm_swMutation_R(object, parent)
}

gaperm_swMutation_R <- function(object, parent)
{
  mutate <- parent <- as.vector(object@population[parent,])
  n <- length(parent)
  m <- sample(1:n, size = 2)
  mutate[m[1]] <- parent[m[2]]
  mutate[m[2]] <- parent[m[1]]
  return(mutate)
}

# Displacement mutation ----

gaperm_dmMutation <- function(object, parent, ...)
{
  if(gaControl("useRcpp"))
    gaperm_dmMutation_Rcpp(object, parent)
  else
    gaperm_dmMutation_R(object, parent)
}

gaperm_dmMutation_R <- function(object, parent)
{
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

# Scramble mutation ----

gaperm_scrMutation <- function(object, parent, ...)
{
  if(gaControl("useRcpp"))
    gaperm_scrMutation_Rcpp(object, parent)
  else
    gaperm_scrMutation_R(object, parent)
}

gaperm_scrMutation_R <- function(object, parent)
{
  parent <- as.vector(object@population[parent,])
  n <- length(parent)
  m <- sort(sample(1:n, size = 2))
  m <- seq(min(m), max(m), by = 1)
  m <- sample(m, replace = FALSE)
  i <- c(setdiff(1:min(m),m), m, setdiff(max(m):n,m))
  mutate <- parent[i]
  return(mutate)
} 

# Variable probability of mutation ----
#
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

ga_pmutation <- function(object, p0 = 0.5, p = 0.01, 
                               T = round(object@maxiter/2), ...)
{
  if(gaControl("useRcpp"))
    ga_pmutation_Rcpp(object, p0, p, T)
  else
    ga_pmutation_R(object, p0, p, T)
}

ga_pmutation_R <- function(object, p0, p, T)
{
  if(missing(p0)) p0 <- 0.5
  if(missing(p))  p <- 0.01
  if(missing(T))  T <- round(object@maxiter/2)
  t <- object@iter
  # linear decay
  # pm <- if(t > T) p else p0 - (p0-p)/T * (t-1)
  # exponential decay
  pm <- (p0 - p)*exp(-2*(t-1)/T) + p
  #
  return(pm)
}

# Probability of selection ----
#
# Probability of selection based on fitness values in vector x with
# selection pressure given by q in (0,1).
# This is used in optim() local search to select which solution should
# be used for starting the algorithm.

optimProbsel <- function(x, q = 0.25, ...)
{
  if(gaControl("useRcpp"))
    optimProbsel_Rcpp(x, q)
  else
    optimProbsel_R(x, q)
}

optimProbsel_R <- function(x, q)
{
  if(missing(q)) q <- 0.25
  x <- as.vector(x)
  n <- length(x)
  eps <- sqrt(.Machine$double.eps)
  # selection pressure parameter
  q <- min(max(eps, q), 1 - eps)
  r <- (n + 1) - rank(x, ties.method = "first", na.last = FALSE)
  # prob <- q*(1-q)^(r-1) * 1/(1-(1-q)^n)
  prob <- exp(log(q) + (r-1.0)*log(1.0-q));
  prob[is.na(x)] <- 0
  prob <- pmin(pmax(0, prob/sum(prob)), 1, na.rm = TRUE)
  return(prob)
}
  