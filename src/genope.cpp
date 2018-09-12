//////////////////////////////////////////////////////////////////////////////
//                                                                          //
//                         GA Rcpp functions                                //
//                                                                          //
//////////////////////////////////////////////////////////////////////////////

#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;

//
// Miscellaneous functions
//

// [[Rcpp::export]]
IntegerVector rank_asR(NumericVector x, bool decreasing = false)
{  
  IntegerVector rank = match(x, clone(x).sort());
  if(decreasing) rank = rank.length()+1 - rank; 
  return rank;
}

/***
x = round(runif(10), 2)
identical(rank(x, ties.method = "min"), rank_asR(x))
identical(rank(x, ties.method = "random"), rank_asR(x))
# with no ties are identical
 
rank(x, ties.method = "min") - rank_asR(x)
length(x)+1-rank(x, ties.method = "min") - rank_asR(x, decreasing = TRUE)

library(microbenchmark)
microbenchmark(rank(x, ties.method = "min"),
               rank_asR(x),
               unit = "relative")
*/

// [[Rcpp::export]]
IntegerVector which_asR(LogicalVector x)
{
  IntegerVector i = Rcpp::seq(0, x.length()-1);
  i = i[x];
  return(i);
}

// [[Rcpp::export]]
IntegerVector c_int(IntegerVector x, IntegerVector y)
{
// combine two integer vectors
  std::vector<int> xy;
  xy.reserve( x.size() + y.size() ); // preallocate memory
  xy.insert( xy.end(), x.begin(), x.end() );
  xy.insert( xy.end(), y.begin(), y.end() );
  return(wrap(xy));
}
  
// [[Rcpp::export]]
NumericVector c_double(NumericVector x, NumericVector y)
{
// combine two numerical vectors
  std::vector<double> xy;
  xy.reserve( x.size() + y.size() ); // preallocate memory
  xy.insert( xy.end(), x.begin(), x.end() );
  xy.insert( xy.end(), y.begin(), y.end() );
  return(wrap(xy));
}

// [[Rcpp::export]]
IntegerVector setdiff_asR(IntegerVector x, IntegerVector y)
{
  // difference of sets x & y (without reordering)
  x = x[duplicated(x) == 0]; x = na_omit(x);
  y = y[duplicated(y) == 0]; y = na_omit(y);
  IntegerVector out(0, NA_INTEGER);
  for(int i=0; i < x.length(); i++)
  {
    if(is_false(any(x[i] == y)))
      out.push_back(x[i]);
  }
  return(out);
}

// [[Rcpp::export]]
IntegerVector intersect_asR(IntegerVector x, IntegerVector y)
{
  // intersection between sets x & y (without reordering)
  x = x[duplicated(x) == 0]; x = na_omit(x);
  y = y[duplicated(y) == 0]; y = na_omit(y);
  IntegerVector out(0, NA_INTEGER);
  for(int i=0; i < x.length(); i++)
  {
    if(is_true(any(x[i] == y)))
      out.push_back(x[i]);
  }
  return(out);
}

/***
x = 1:5
foo(x)
setdiff(x, c(3,5))
setdiff_asR(x, c(3,5))
setdiff(rev(x), c(3,5))
setdiff_asR(rev(x), c(3,5))
*/

// [[Rcpp::export]]
double round_double(double x, int n)
{
  // round function for double
  // for vectors Rcpp sugar already has round 
	int d = 0;
	if( ((x * pow(10.0, n)) - floor(x * pow(10.0, n))) >= 0.5) d = 1;
	x = floor(x * pow(10.0, n) + d) / pow(10.0, n);
	return x;
}

// 
// Generic GA operators
//

// [[Rcpp::export]]
List ga_lrSelection_Rcpp(RObject object, double r = NA_REAL, double q = NA_REAL)
{
  NumericVector fitness = object.slot("fitness");
  NumericMatrix pop = object.slot("population");
  int           popSize = pop.nrow();
  int           n       = pop.ncol();
  NumericMatrix newpop(popSize, n);
  if(std::isnan(r))  r = 2.0/(popSize*(popSize-1));
  if(std::isnan(q))  q = 2.0/popSize;
  double eps = std::numeric_limits<double>::epsilon();
  NumericVector rank = as<NumericVector>(rank_asR(fitness, true));
  NumericVector prob = 1 + q - (rank-1)*r;
                prob[is_na(prob)] = eps;
                prob = pmin(pmax(0.0, prob/sum(prob)), 1.0);
  IntegerVector seq = Rcpp::seq(0, popSize-1);
  IntegerVector sel = Rcpp::sample(seq, popSize, true, prob);
  for(int i=0; i < sel.length(); i++)
  {
    newpop(i,_) = pop(sel[i],_);
  }
  fitness = fitness[sel];
  
  List out = List::create(Rcpp::Named("population") = newpop,
                          Rcpp::Named("fitness") = fitness);
  return out;  
}

/***
library(GA)
library(microbenchmark)
f <- function(x)  abs(x)+cos(x)
fitness <- function(x) -f(x)
GA <- ga(type = "real-valued", fitness = fitness, 
         popSize = 10, lower = -20, upper = 20, maxiter = 10, seed = 1)

set.seed(123)
out1 = ga_lrSelection_R(GA)
set.seed(123)
out2 = ga_lrSelection_Rcpp(GA)
identical(out1, out2)

microbenchmark(ga_lrSelection_R(GA),
               ga_lrSelection_Rcpp(GA),
               unit = "relative")
*/


// [[Rcpp::export]]
List ga_nlrSelection_Rcpp(RObject object, double q = NA_REAL)
{
  if(std::isnan(q))  q = 0.25;
  NumericVector fitness = object.slot("fitness");
  NumericMatrix pop = object.slot("population");
  int           popSize = pop.nrow();
  int           n       = pop.ncol();
  NumericMatrix newpop(popSize, n);
  NumericVector rank = as<NumericVector>(rank_asR(fitness, true));
  double        eps = std::numeric_limits<double>::epsilon();
  NumericVector prob = exp(log(q)+(rank-1.0)*log(1-q));
                prob[is_na(prob)] = eps;
                prob = pmin(pmax(0.0, prob/sum(prob)), 1.0);
  IntegerVector seq = Rcpp::seq(0, popSize-1);
  IntegerVector sel = Rcpp::sample(seq, popSize, true, prob);
  for(int i=0; i < sel.length(); i++) 
  {
    newpop(i,_) = pop(sel[i],_);
  }
  fitness = fitness[sel];
  
  List out = List::create(Rcpp::Named("population") = newpop,
                          Rcpp::Named("fitness") = fitness);
  return out;  
}

/***
library(GA)
library(microbenchmark)
f <- function(x)  abs(x)+cos(x)
fitness <- function(x) -f(x)
GA <- ga(type = "real-valued", fitness = fitness, seed = 1, popSize = 10, lower = -20, upper = 20, maxiter = 1)

set.seed(12)
out1 = ga_nlrSelection(GA)
set.seed(12)
out2 = ga_nlrSelection_Rcpp(GA)
identical(out1, out2)
 
GA <- ga(type = "real-valued", fitness = fitness, popSize = 100,
         lower = -20, upper = 20, run = 1, monitor = FALSE)
microbenchmark(ga_nlrSelection_R(GA),
               ga_nlrSelection_Rcpp(GA),
               unit = "relative")
*/

// [[Rcpp::export]]
List ga_rwSelection_Rcpp(RObject object)
{
  NumericVector fitness = object.slot("fitness");
  NumericMatrix pop = object.slot("population");
  int           popSize = pop.nrow();
  int           n       = pop.ncol();
  NumericMatrix newpop(popSize, n);
  double        eps = std::numeric_limits<double>::epsilon();
  NumericVector prob = abs(fitness);
                prob[is_na(prob)] = eps;
                prob = pmin(pmax(0.0, prob/sum(prob)), 1.0);
  IntegerVector seq = Rcpp::seq(0, popSize-1);
  IntegerVector sel = Rcpp::sample(seq, popSize, true, prob);
  
  for(int i=0; i < sel.length(); i++) 
  {
    newpop(i,_) = pop(sel[i],_);
  }
  fitness = fitness[sel];
  
  List out = List::create(Rcpp::Named("population") = newpop,
                          Rcpp::Named("fitness") = fitness);
  return out;  
}

/***
library(GA)
library(microbenchmark)

f <- function(x)  abs(x)+cos(x)
GA <- ga(type = "real-valued", fitness = function(x) -f(x), 
         popSize = 50, lower = -20, upper = 20, maxiter = 10, seed = 1)

set.seed(123)
out1 = ga_rwSelection(GA)
set.seed(123)
out2 = ga_rwSelection_Rcpp(GA)
identical(out1, out2)

microbenchmark(ga_rwSelection_R(GA),
               ga_rwSelection_Rcpp(GA),
               unit = "relative")
 
f <- function(x1, x2) 
  { 20 + x1^2 + x2^2 - 10*(cos(2*pi*x1) + cos(2*pi*x2)) }
GA <- ga(type = "real-valued", fitness = function(x) -f(x[1],x[2]),  
         popSize = 50, lower = c(-5.12,-5.12), upper = c(5.12,5.12),
         maxiter = 10, seed = 1)
set.seed(123)
out1 = ga_rwSelection(GA)
set.seed(123)
out2 = ga_rwSelection_Rcpp(GA)
identical(out1, out2)
 
microbenchmark(ga_rwSelection_R(GA),
               ga_rwSelection_Rcpp(GA),
               unit = "relative")
*/

// [[Rcpp::export]]
List ga_tourSelection_Rcpp(RObject object, double k = NA_REAL)
{
  if(std::isnan(k))  k = 3;
  NumericVector fitness = object.slot("fitness");
  NumericMatrix pop = object.slot("population");
  int           popSize = pop.nrow();
  int           n       = pop.ncol();
  NumericMatrix newpop(popSize, n);
  IntegerVector seq = Rcpp::seq(0, popSize-1);
  IntegerVector sel  (popSize, NA_INTEGER);
  for(int i=0; i < sel.length(); i++) 
  {  
    IntegerVector s = Rcpp::sample(seq, k, false);
    int j = Rcpp::which_max(as<NumericVector>(fitness[s]));
    sel[i] = s[j];
  }
  for(int i=0; i < sel.length(); i++) 
  {
    newpop(i,_) = pop(sel[i],_);
  }
  fitness = fitness[sel];
  
  List out = List::create(Rcpp::Named("population") = newpop,
                          Rcpp::Named("fitness") = fitness);
  return out;  
}

/***
library(GA)
library(microbenchmark)

f <- function(x)  abs(x)+cos(x)
fitness <- function(x) -f(x)
GA <- ga(type = "real-valued", fitness = fitness, 
         popSize = 20, lower = -20, upper = 20, maxiter = 10, seed = 1)

set.seed(123)
out1 = ga_tourSelection_R(GA)
set.seed(123)
out2 = ga_tourSelection_Rcpp(GA)
identical(out1, out2)
 
microbenchmark(ga_tourSelection_R(GA),
               ga_tourSelection_Rcpp(GA),
               unit = "relative")
*/

// [[Rcpp::export]]
List ga_spCrossover_Rcpp(RObject object, IntegerVector parents)
{
  NumericMatrix pop = object.slot("population");
  int           n = pop.ncol();
  arma::rowvec  pop0 = pop(parents[0]-1,_);
  arma::rowvec  pop1 = pop(parents[1]-1,_);
  NumericVector fitness = object.slot("fitness");
  arma::mat     children(2, n, arma::fill::zeros);
  NumericVector fitnessChildren(2, NA_REAL);
  IntegerVector seq = Rcpp::seq(0, n);
  int crossOverPoint = Rcpp::sample(seq, 1, false)[0];
  if(crossOverPoint == 0)
  { 
    children = arma::join_cols(pop1, pop0);
    fitnessChildren = fitness[rev(parents-1)];
  }
  else if(crossOverPoint == n)
  { 
    children = arma::join_cols(pop0, pop1);
    fitnessChildren = fitness[parents-1];
  }
  else
  { 
    arma::rowvec x1_1 = pop0(arma::span(0,crossOverPoint-1));
    arma::rowvec x1_2 = pop1(arma::span(crossOverPoint,n-1));
    arma::rowvec x2_1 = pop1(arma::span(0,crossOverPoint-1));
    arma::rowvec x2_2 = pop0(arma::span(crossOverPoint,n-1));
    arma::rowvec x1   = arma::join_rows(x1_1, x1_2);
    arma::rowvec x2   = arma::join_rows(x2_1, x2_2);
    children = arma::join_cols(x1, x2);
  }

  List out = List::create(Rcpp::Named("children") = children,
                          Rcpp::Named("fitness") = fitnessChildren);
  return out;  
} 
 
/***
library(GA)
library(microbenchmark)
f <- function(x) { 20 + x[1]^2 + x[2]^2 - 10*(cos(2*pi*x[1]) + cos(2*pi*x[2])) }
GA <- ga(type = "real-valued", fitness = f, lower = c(-5.12,-5.12), upper = c(5.12,5.12), maxiter = 10, seed = 1)

i = c(1,2); GA@population[i,]
set.seed(1)
(out1 = ga_spCrossover_R(GA, i))
set.seed(1)
(out2 = ga_spCrossover_Rcpp(GA, i))
identical(out1$children, out2$children)

microbenchmark(ga_spCrossover_R(GA, 1:2),
               ga_spCrossover_Rcpp(GA, 1:2),
               unit = "relative")
*/

// 
// "binary" GA genetic operators
// 

// [[Rcpp::export]]
NumericMatrix gabin_Population_Rcpp(RObject object)
{
  int popSize = object.slot("popSize");
  int nBits   = object.slot("nBits");
  NumericMatrix pop(popSize, nBits);
  std::fill( pop.begin(), pop.end(), NumericVector::get_na() );

  for(int j=0; j < nBits; j++) 
  {
    pop(_,j) = round(Rcpp::runif(popSize, 0, 1), 0);
  }
    
  return pop;
}

/***
library(GA)
library(microbenchmark)

data(fat, package = "UsingR")
mod <- lm(body.fat.siri ~ age + weight + height + neck + chest + abdomen +
          hip + thigh + knee + ankle + bicep + forearm + wrist, data = fat)
x <- model.matrix(mod)[,-1]
y <- model.response(mod$model)
fitness <- function(string)
{ 
  inc <- which(string==1)
  X <- cbind(1, x[,inc])
  mod <- lm.fit(X, y)
  class(mod) <- "lm"
  -AIC(mod)
}

GA <- ga("binary", fitness = fitness, nBits = ncol(x), maxiter = 10, seed = 1)

set.seed(1)
out1 <- gabin_Population_R(GA)
set.seed(1)
out2 <- gabin_Population_Rcpp(GA)
identical(out1, out2)

microbenchmark(gabin_Population_R(GA),
               gabin_Population_Rcpp(GA),
               unit = "relative")
*/

// [[Rcpp::export]]
List gabin_uCrossover_Rcpp(RObject object, IntegerVector parents)
{
  NumericMatrix pop = object.slot("population");
  int           n = pop.ncol();
  NumericVector fitness(2, NA_REAL);
  NumericMatrix children(2, n);
  NumericVector u = Rcpp::runif(n, 0, 1);
  for(int j=0; j < n; j++) 
  {
    if(u[j] > 0.5) 
      { children(0,j) = pop(parents[1]-1,j);
        children(1,j) = pop(parents[0]-1,j);
      }
    else
      { children(0,j) = pop(parents[0]-1,j);
        children(1,j) = pop(parents[1]-1,j);
      }
  }

  List out = List::create(Rcpp::Named("children") = children,
                          Rcpp::Named("fitness")  = fitness);
  return out;  
} 
 
/***
library(GA)
library(microbenchmark)
data(fat, package = "UsingR")
mod <- lm(body.fat.siri ~ age + weight + height + neck + chest + abdomen +
          hip + thigh + knee + ankle + bicep + forearm + wrist, data = fat)
x <- model.matrix(mod)[,-1]
y <- model.response(mod$model)
fitness <- function(string)
{ 
  inc <- which(string==1)
  X <- cbind(1, x[,inc])
  mod <- lm.fit(X, y)
  class(mod) <- "lm"
  -AIC(mod)
}

GA = ga("binary", fitness = fitness, nBits = ncol(x), maxiter = 10, seed = 1)

i <- c(1,2); GA@population[i,]
set.seed(1)
(out1 <- gabin_uCrossover_R(GA, i))
set.seed(1)
(out2 <- gabin_uCrossover_Rcpp(GA, i))
out1$children - out2$children
 
microbenchmark(gabin_uCrossover_R(GA, i),
               gabin_uCrossover_Rcpp(GA, i),
               unit = "relative")
*/

// [[Rcpp::export]]
NumericVector gabin_raMutation_Rcpp(RObject object, int parent)
{
  NumericMatrix pop = object.slot("population");
  int           nBits = object.slot("nBits");
  NumericVector mutate = pop(parent-1,_);
  IntegerVector seq = Rcpp::seq(0, nBits-1);
  IntegerVector j = Rcpp::sample(seq, 1, true);
  mutate[j] = std::abs(as<double>(mutate[j]) - 1.0);
  return mutate;
}

/***
library(GA)
library(microbenchmark)
data(fat, package = "UsingR")
mod <- lm(body.fat.siri ~ age + weight + height + neck + chest + abdomen +
          hip + thigh + knee + ankle + bicep + forearm + wrist, data = fat)
x <- model.matrix(mod)[,-1]
y <- model.response(mod$model)
fitness <- function(string)
{ 
  inc <- which(string==1)
  X <- cbind(1, x[,inc])
  mod <- lm.fit(X, y)
  class(mod) <- "lm"
  -AIC(mod)
}

GA = ga("binary", fitness = fitness, nBits = ncol(x), maxiter = 10, seed = 1)

i <- 14; GA@population[i,]
set.seed(3)
(out1 <- gabin_raMutation_R(GA, i))
set.seed(3)
(out2 <- gabin_raMutation_Rcpp(GA, i))
out1-out2

microbenchmark(gabin_raMutation_R(GA, i),
               gabin_raMutation_Rcpp(GA, i),
               unit = "relative")
*/


// 
// "real-valued" GA genetic operators
// 

// [[Rcpp::export]]
NumericMatrix gareal_Population_Rcpp(RObject object)
{
  int           popSize = object.slot("popSize");
  NumericVector lower   = object.slot("lower");
  NumericVector upper   = object.slot("upper");
  int nvars = lower.length();
  NumericMatrix pop(popSize, nvars);
  std::fill( pop.begin(), pop.end(), NumericVector::get_na() );

  for(int j=0; j < nvars; j++) 
  {
    pop(_,j) = Rcpp::runif(popSize, lower[j], upper[j]);
  }
    
  return pop;
}

/***
library(GA)
library(microbenchmark)
f <- function(x) -(abs(x)+cos(x))
GA <- ga(type = "real-valued", fitness = f, lower = -20, upper = 20, maxiter = 10, seed = 1)
f <- function(x) { 20 + x[1]^2 + x[2]^2 - 10*(cos(2*pi*x[1]) + cos(2*pi*x[2])) }
GA <- ga(type = "real-valued", fitness = f, lower = c(-5.12,-5.12), upper = c(5.12,5.12), maxiter = 10, seed = 1)

set.seed(123)
out1 <- gareal_Population_R(GA)
set.seed(123)
out2 <- gareal_Population_Rcpp(GA)
identical(out1, out2)

microbenchmark(gareal_Population_R(GA),
               gareal_Population_Rcpp(GA),
               unit = "relative")
*/

// [[Rcpp::export]]
List gareal_lsSelection_Rcpp(RObject object)
{
  NumericVector fitness = object.slot("fitness");
  NumericMatrix pop = object.slot("population");
  int           popSize = pop.nrow();
  int           n = pop.ncol();
  double        eps = std::numeric_limits<double>::epsilon();
  NumericMatrix newpop(popSize, n);
  NumericVector f = clone(fitness);
  NumericVector fo = na_omit(f);
  double fmin = min(fo);
  if(fmin < 0) 
  { 
    f = f - fmin;
    fo = na_omit(f);
    fmin = min(fo); 
  }
  double fave = mean(fo);
  double fmax = max(fo);
  double sfactor = 2.0; // scaling factor
  // transform f -> f' = a*f + b such that
  double delta, a, b;
  if(fmin > (sfactor*fave - fmax)/(sfactor-1.0))
    { // ave(f) = ave(f')
      // 2*ave(f') = max(f')
      delta = fmax - fave;
      a = (sfactor - 1.0)*fave/delta;
      b = fave * (fmax - sfactor*fave)/delta;
    }
  else
    { // ave(f) = ave(f')
      // min(f') = 0
      delta = fave - fmin;
      a = fave/delta;
      b = -1.0*fmin*fave/delta;
    }
  NumericVector fscaled = a*f + b;
  NumericVector prob = abs(fscaled); 
                prob[is_na(prob)] = eps;
                prob[is_infinite(prob)] = eps;
                prob = pmin(pmax(0.0, prob/sum(prob)), 1.0);
  IntegerVector seq = Rcpp::seq(0, popSize-1);
  IntegerVector sel = Rcpp::sample(seq, popSize, true, prob);
  for(int i=0; i < sel.length(); i++)
  {
    newpop(i,_) = pop(sel[i],_);
  }
  fitness = fitness[sel];
  
  List out = List::create(Rcpp::Named("population") = newpop,
                          Rcpp::Named("fitness") = fitness);
  return out;  
}

/***
library(GA)
library(microbenchmark)
f <- function(x) -(abs(x)+cos(x))
GA <- ga(type = "real-valued", fitness = f, lower = -20, upper = 20, maxiter = 10, seed = 1)
f <- function(x) { 20 + x[1]^2 + x[2]^2 - 10*(cos(2*pi*x[1]) + cos(2*pi*x[2])) }
GA <- ga(type = "real-valued", fitness = f, lower = c(-5.12,-5.12), upper = c(5.12,5.12), maxiter = 10, seed = 1)

set.seed(123)
(out1 <- GA:::gareal_lsSelection_R(GA))
set.seed(123)
(out2 <- GA:::gareal_lsSelection_Rcpp(GA))
identical(out1, out2)

microbenchmark(GA:::gareal_lsSelection_R(GA),
               GA:::gareal_lsSelection_Rcpp(GA),
               unit = "relative")
*/

// [[Rcpp::export]]
List gareal_sigmaSelection_Rcpp(RObject object)
{
  NumericVector fitness = object.slot("fitness");
  NumericMatrix pop = object.slot("population");
  int           popSize = pop.nrow();
  int           n = pop.ncol();
  NumericMatrix newpop(popSize, n);
  // NumericVector f = clone(fitness);
  // double        mf = mean(na_omit(f)); 
  // double        sf = sd(na_omit(f)); 
  double        mf = mean(na_omit(fitness)); 
  double        sf = sd(na_omit(fitness)); 
  double        eps = std::numeric_limits<double>::epsilon();
  NumericVector fscaled = pmax(fitness - (mf - 2*sf), 0);
  NumericVector prob = abs(fscaled);
                prob[is_na(prob)] = eps;
                prob = pmin(pmax(0.0, prob/sum(prob)), 1.0);
  IntegerVector seq = Rcpp::seq(0, popSize-1);
  IntegerVector sel = Rcpp::sample(seq, popSize, true, prob);
  for(int i=0; i < sel.length(); i++)
  {
    newpop(i,_) = pop(sel[i],_);
  }
  fitness = fitness[sel];
  
  List out = List::create(Rcpp::Named("population") = newpop,
                          Rcpp::Named("fitness") = fitness);
  return out;  
}

/***
library(GA)
library(microbenchmark)
f <- function(x) -(abs(x)+cos(x))
GA <- ga(type = "real-valued", fitness = f, lower = -20, upper = 20, maxiter = 10, seed = 1)
f <- function(x) { 20 + x[1]^2 + x[2]^2 - 10*(cos(2*pi*x[1]) + cos(2*pi*x[2])) }
GA <- ga(type = "real-valued", fitness = f, lower = c(-5.12,-5.12), upper = c(5.12,5.12), maxiter = 10, seed = 1)

set.seed(1)
(out1 <- gareal_sigmaSelection_R(GA))
set.seed(1)
(out2 <- gareal_sigmaSelection_Rcpp(GA))
identical(out1, out2)

microbenchmark(gareal_sigmaSelection_R(GA),
               gareal_sigmaSelection_Rcpp(GA),
               unit = "relative")
*/

// [[Rcpp::export]]
List gareal_waCrossover_Rcpp(RObject object, IntegerVector parents)
{
  NumericMatrix pop = object.slot("population");
  int           n = pop.ncol();
  double        a = R::runif(0,1);
  NumericVector fitness(2, NA_REAL);
  NumericMatrix children(2, n);
  children(0,_) = a*pop(parents[0]-1,_) + (1-a)*pop(parents[1]-1,_);
  children(1,_) = a*pop(parents[1]-1,_) + (1-a)*pop(parents[0]-1,_);
  
  List out = List::create(Rcpp::Named("children") = children,
                          Rcpp::Named("fitness")  = fitness);
  return out;  
}  

/***
library(GA)
library(microbenchmark)
f <- function(x)  abs(x)+cos(x)
fitness <- function(x) -f(x)
GA <- ga(type = "real-valued", fitness = fitness, seed = 1, lower = -20, upper = 20, maxiter = 10)

set.seed(123)
(out1 <- gareal_waCrossover_R(GA, c(11,33)))
set.seed(123)
(out2 <- gareal_waCrossover_Rcpp(GA, c(11,33)))
identical(out1$children, out2$children)

microbenchmark(gareal_waCrossover_R(GA, 1:2),
               gareal_waCrossover_Rcpp(GA, 1:2),
               unit = "relative")
*/

// [[Rcpp::export]]
List gareal_laCrossover_Rcpp(RObject object, IntegerVector parents)
{
  NumericMatrix pop = object.slot("population");
  int           n = pop.ncol();
  NumericVector a = Rcpp::runif(n);
  NumericVector fitness(2, NA_REAL);
  NumericMatrix children(2, n);
  children(0,_) = a*pop(parents[0]-1,_) + (1-a)*pop(parents[1]-1,_);
  children(1,_) = a*pop(parents[1]-1,_) + (1-a)*pop(parents[0]-1,_);

  List out = List::create(Rcpp::Named("children") = children,
                          Rcpp::Named("fitness") = fitness);
  return out;  
}  

/***
library(GA)
library(microbenchmark)
f <- function(x)  abs(x)+cos(x)
fitness <- function(x) -f(x)
GA <- ga(type = "real-valued", fitness = fitness, seed = 1, lower = -20, upper = 20, maxiter = 10)

set.seed(123)
(out1 <- gareal_laCrossover_R(GA, 1:2))
set.seed(123)
(out2 <- gareal_laCrossover_Rcpp(GA, 1:2))
identical(out1, out2)

microbenchmark(gareal_laCrossover_R(GA, 1:2),
               gareal_laCrossover_Rcpp(GA, 1:2),
               unit = "relative")
*/

// [[Rcpp::export]]
List gareal_blxCrossover_Rcpp(RObject object, IntegerVector parents, double a)
{
  // double a = 0.5;
  NumericMatrix pop = object.slot("population");
  int           n = pop.ncol();
  NumericVector fitness(2, NA_REAL);
  NumericVector lower = object.slot("lower");
  NumericVector upper = object.slot("upper");
  
  NumericMatrix children(2, n);
  for(int i=0; i < n; i++) 
  {
    NumericVector x = NumericVector::create(pop(parents[0]-1,i), 
                                            pop(parents[1]-1,i));
    std::sort(x.begin(), x.end());
    double xl = std::max(x[0] - a*(x[1]-x[0]), lower[i]);
    double xu = std::min(x[1] + a*(x[1]-x[0]), upper[i]);
    children(_,i) = Rcpp::runif(2, xl, xu);
  }

  List out = List::create(Rcpp::Named("children") = children,
                          Rcpp::Named("fitness") = fitness);
  return out;  
}

/***
library(GA)
library(microbenchmark)
f <- function(x) -(abs(x)+cos(x))
GA <- ga(type = "real-valued", fitness = f, seed = 1, lower = -20, upper = 20, maxiter = 10)

set.seed(123)
(out1 <- gareal_blxCrossover_R(GA, 1:2))
set.seed(123)
(out2 <- gareal_blxCrossover_Rcpp(GA, 1:2))
identical(out1$children, out2$children)

i <- sample(1:GA@popSize,2)
set.seed(123)
(out1 <- gareal_blxCrossover_R(GA, i))
set.seed(123)
(out2 <- gareal_blxCrossover_Rcpp(GA, i))
identical(out1$children, out2$children)

microbenchmark(gareal_blxCrossover_R(GA, 1:2),
               gareal_blxCrossover_Rcpp(GA, 1:2),
               unit = "relative")
*/

// [[Rcpp::export]]
List gareal_laplaceCrossover_Rcpp(RObject object, 
                                  IntegerVector parents, 
                                  NumericVector a,
                                  NumericVector b)
{
  NumericMatrix pop = object.slot("population");
  int           n = pop.ncol();
  if(a.length() == 1)  a = rep(a(0), n);
  if(b.length() == 1)  b = rep(b(0), n);
  NumericVector lower = object.slot("lower");
  NumericVector upper = object.slot("upper");
  NumericVector fitness(2, NA_REAL);
  NumericMatrix children(2, n);
  NumericVector r = Rcpp::runif(n, 0, 1);
  NumericVector u = Rcpp::runif(n, 0, 1);
  NumericVector beta = a + ifelse(r > 0.5, b*log(u), -b*log(u));
  NumericVector bpar = beta*abs( pop(parents[0]-1,_) - pop(parents[1]-1,_) );
  children(0,_) = pmin(pmax(lower, pop(parents[0]-1,_) + bpar), upper);
  children(1,_) = pmin(pmax(lower, pop(parents[1]-1,_) + bpar), upper);

  List out = List::create(Rcpp::Named("children") = children,
                          Rcpp::Named("fitness") = fitness);
  return out;
}

/***
library(GA)
library(microbenchmark)
f <- function(x) -(abs(x)+cos(x))
GA <- ga(type = "real-valued", fitness = f, lower = -20, upper = 20, maxiter = 10, seed = 1)
f <- function(x) { 20 + x[1]^2 + x[2]^2 - 10*(cos(2*pi*x[1]) + cos(2*pi*x[2])) }
GA <- ga(type = "real-valued", fitness = f, lower = c(-5.12,-5.12), upper = c(5.12,5.12), maxiter = 10, seed = 1)

set.seed(123)
(out1 = GA:::gareal_laplaceCrossover_R(GA, c(11,44), a = 1, b = 0.1))
set.seed(123)
(out2 = gareal_laplaceCrossover_Rcpp(GA, c(11,44), a = 1, b = 0.1))
identical(out1$children, out2$children)

microbenchmark(GA:::gareal_laplaceCrossover_R(GA, 1:2, a = 1, b = 0.1),
               gareal_laplaceCrossover_Rcpp(GA, 1:2, a = 1, b = 0.1),
               unit = "relative")
*/

// [[Rcpp::export]]
NumericVector gareal_raMutation_Rcpp(RObject object, int parent)
{
  NumericMatrix pop = object.slot("population");
  int           n = pop.ncol();
  NumericVector lower = object.slot("lower");
  NumericVector upper = object.slot("upper");
  NumericVector mutate = pop(parent-1,_);
  IntegerVector seq = Rcpp::seq(0, n-1);
  IntegerVector j = Rcpp::sample(seq, 1, true);
  NumericVector u = Rcpp::runif(1, as<double>(lower[j]), 
                                   as<double>(upper[j]));
  mutate[j] = u;
  return mutate;  
}  

/***
library(GA)
library(microbenchmark)
f <- function(x) -(abs(x)+cos(x))
GA <- ga(type = "real-valued", fitness = f, lower = -20, upper = 20, maxiter = 10, seed = 1)
f <- function(x) { 20 + x[1]^2 + x[2]^2 - 10*(cos(2*pi*x[1]) + cos(2*pi*x[2])) }
GA <- ga(type = "real-valued", fitness = f, lower = c(-5.12,-5.12), upper = c(5.12,5.12), maxiter = 10, seed = 1)

set.seed(123)
(out1 = gareal_raMutation_R(GA, 10))
set.seed(123)
(out2 = gareal_raMutation_Rcpp(GA, 10))
identical(out1, out2)

microbenchmark(gareal_raMutation_R(GA, 10),
               gareal_raMutation_Rcpp(GA, 10),
               unit = "relative")
*/

// [[Rcpp::export]]
NumericVector gareal_nraMutation_Rcpp(RObject object, int parent)
{
  NumericMatrix pop = object.slot("population");
  int           n = pop.ncol();
  NumericVector lower = object.slot("lower");
  NumericVector upper = object.slot("upper");
  double        iter = object.slot("iter");
  double        maxiter = object.slot("maxiter");
  NumericVector mutate = pop(parent-1,_);
  double        g = 1 - iter/maxiter; // dempening factor
  IntegerVector seq = Rcpp::seq(0, n-1);
  IntegerVector j = Rcpp::sample(seq, 1, true);
  NumericVector u = Rcpp::runif(2);
  NumericVector m = mutate[j];
  if(u[0] < 0.5)
    { NumericVector sa = ( mutate[j] - lower[j] )*(1 - pow(u[1],g));
      m += -sa;
    }
  else
    { NumericVector sa = ( upper[j] - mutate[j])*(1 - pow(u[1],g));
      m +=  sa;
    }
  mutate[j] = m;
  return mutate;  
}  

/***
library(GA)
library(microbenchmark)
# f <- function(x) -(abs(x)+cos(x))
# GA <- ga(type = "real-valued", fitness = f, lower = -20, upper = 20, maxiter = 10, seed = 1)
f <- function(x) { 20 + x[1]^2 + x[2]^2 - 10*(cos(2*pi*x[1]) + cos(2*pi*x[2])) }
GA <- ga(type = "real-valued", fitness = f, lower = c(-5.12,-5.12), upper = c(5.12,5.12), maxiter = 10, seed = 1)
GA@maxiter <- 100

set.seed(1234)
(out1 = GA:::gareal_nraMutation_R(GA, 10))
set.seed(1234)
(out2 = GA:::gareal_nraMutation_Rcpp(GA, 10))
identical(out1, out2)

microbenchmark(GA:::gareal_nraMutation_R(GA, 10),
               GA:::gareal_nraMutation_Rcpp(GA, 10),
               unit = "relative")
*/

// [[Rcpp::export]]
NumericVector gareal_rsMutation_Rcpp(RObject object, int parent)
{
  NumericMatrix pop = object.slot("population");
  int           n = pop.ncol();
  NumericVector lower = object.slot("lower");
  NumericVector upper = object.slot("upper");
  NumericVector mutate = pop(parent-1,_);
  double        iter = object.slot("iter");
  double        maxiter = object.slot("maxiter");
  double        dempeningFactor = 1 - iter/maxiter;
  double        direction = 0.0;
  if(R::runif(0,1) < 0.5)  
    direction = -1.0; else direction = 1.0;
  NumericVector value = (upper - lower)*0.67;

  for(int j=0; j < n; j++) 
  {
    mutate[j] += direction*dempeningFactor*value[j];
    if((mutate[j] < lower[j]) | (mutate[j] > upper[j]))
      { 
        NumericVector m = Rcpp::runif(1, lower[j], upper[j]);
        mutate[j] = m[0];
      }
  }
  
  return mutate;
}  

/***
library(GA)
library(microbenchmark)
# f <- function(x) -(abs(x)+cos(x))
# GA <- ga(type = "real-valued", fitness = f, lower = -20, upper = 20, maxiter = 10, seed = 1)
f <- function(x) { 20 + x[1]^2 + x[2]^2 - 10*(cos(2*pi*x[1]) + cos(2*pi*x[2])) }
GA <- ga(type = "real-valued", fitness = f, lower = c(-5.12,-5.12), upper = c(5.12,5.12), maxiter = 10, seed = 1)
GA@maxiter <- 100
i = 10
set.seed(1)
(out1 = gareal_rsMutation_R(GA, i))
set.seed(1)
(out2 = gareal_rsMutation_Rcpp(GA, i))
identical(out1, out2)

microbenchmark(gareal_rsMutation_R(GA, 10),
               gareal_rsMutation_Rcpp(GA, 10),
               unit = "relative")
*/

// [[Rcpp::export]]
NumericVector gareal_powMutation_Rcpp(RObject object, 
                                      int parent, 
                                      NumericVector pow)
{
  NumericMatrix  pop = object.slot("population");
  int            n = pop.ncol();
  if(pow.length() == 1)  pow = rep(pow(0), n);
  NumericVector  lower = object.slot("lower");
  NumericVector  upper = object.slot("upper");
  NumericVector  mutate = pop(parent-1,_);
  NumericVector  t = (mutate - lower)/(upper - mutate);
  double         u = R::runif(0,1);
  double         s;
  for(int j=0; j < n; j++) 
  {
    s = std::pow(u, pow(j));
    if(R::runif(0,1) < t[j])
      mutate[j] += -s*(mutate[j] - lower[j]);
    else
      mutate[j] += s*(upper[j] - mutate[j]);
  }
  return mutate;
}  

/***
library(GA)
library(microbenchmark)
# f <- function(x) -(abs(x)+cos(x))
# GA <- ga(type = "real-valued", fitness = f, lower = -20, upper = 20, maxiter = 10, seed = 1)
f <- function(x) { 20 + x[1]^2 + x[2]^2 - 10*(cos(2*pi*x[1]) + cos(2*pi*x[2])) }
GA <- ga(type = "real-valued", fitness = f, lower = c(-5.12,-5.12), upper = c(5.12,5.12), maxiter = 10, seed = 1)
GA@maxiter <- 100

i = 30
set.seed(123)
(out1 = GA:::gareal_powMutation_R(GA, i, 10))
set.seed(123)
(out2 = gareal_powMutation_Rcpp(GA, i, 10))
identical(out1, out2)

microbenchmark(GA:::gareal_powMutation_R(GA, i, 10),
               gareal_powMutation_Rcpp(GA, i, 10),
               unit = "relative")
*/

// 
// "permutation" GA genetic operators
// 

// [[Rcpp::export]]
IntegerMatrix gaperm_Population_Rcpp(RObject object)
{
  int popSize = object.slot("popSize");
  int lower   = object.slot("lower");
  int upper   = object.slot("upper");
  IntegerVector s = Rcpp::seq(lower, upper);
  int           n = s.length();
  IntegerMatrix pop(popSize, n);

  for(int i=0; i < popSize; i++)
  {
    pop(i,_) = Rcpp::sample(s, n, false);
  }
    
  return pop;
}

/***
library(GA)
library(microbenchmark)
data(eurodist)
D <- as.matrix(eurodist)
tourLength <- function(tour, distMatrix)
{ 
  tour <- c(tour, tour[1])
  route <- embed(tour, 2)[,2:1]
  sum(distMatrix[route])
}
tspFitness <- function(tour, ...) { 1/tourLength(tour, ...) }
        
GA <- ga(type = "permutation", 
         fitness = tspFitness, distMatrix = D,
         lower = 1, upper = attr(eurodist, "Size"), 
         # popSize = 50, maxiter = 5000, run = 500, 
         popSize = 10, maxiter = 10)

set.seed(1)
(out1 = gaperm_Population_R(GA))
set.seed(1)
(out2 = gaperm_Population_Rcpp(GA))
identical(out1, out2)

microbenchmark(gaperm_Population_R(GA),
               gaperm_Population_Rcpp(GA),
               unit = "relative")
*/

// [[Rcpp::export]]
List gaperm_cxCrossover_Rcpp(RObject object, IntegerVector parents)
{
  IntegerMatrix pop = object.slot("population");
  int           n = pop.ncol();
  int           k = 1; // cx point
  NumericVector fitness(2, NA_REAL);
  IntegerMatrix parentsPop(2,n);
                parentsPop(0,_) = pop(parents[0]-1,_);
                parentsPop(1,_) = pop(parents[1]-1,_);
  IntegerMatrix children(2, n);
  IntegerVector ALL = Rcpp::seq(0,n-1);
  while(ALL.length() > 0)
  {
    int i    = ALL(0);
    // perform a cycle
    int base = parentsPop(0,i);
    int vi   = parentsPop(1,i);
    IntegerVector I = IntegerVector::create(i);
    
    while(vi != base)
    {
      i  = as<int>(which_asR(parentsPop(0,_) == parentsPop(1,i)));
      vi = parentsPop(1,i);
      I.push_back(i);
    }
    ALL = setdiff_asR(ALL,I);
    if((k % 2) == 1) 
    { 
      for(int j=0; j < I.length(); j++) 
      {
        children(0,I(j)) = parentsPop(0,I(j));
        children(1,I(j)) = parentsPop(1,I(j));
      }
    }
    else
    {
      for(int j=0; j < I.length(); j++) 
      {
        children(0,I(j)) = parentsPop(1,I(j));
        children(1,I(j)) = parentsPop(0,I(j));
      }
    }
    k += 1;
  } 

  List out = List::create(Rcpp::Named("children") = children,
                          Rcpp::Named("fitness")  = fitness);
  return out;  
} 
 
/***
library(GA)
library(microbenchmark)

data(eurodist)
D <- as.matrix(eurodist)
tourLength <- function(tour, distMatrix)
{ 
  tour <- c(tour, tour[1])
  route <- embed(tour, 2)[,2:1]
  sum(distMatrix[route])
}
tspFitness <- function(tour, ...) { 1/tourLength(tour, ...) }
        
GA <- ga(type = "permutation", 
         fitness = tspFitness, distMatrix = D,
         lower = 1, upper = attr(eurodist, "Size"), 
         # popSize = 50, maxiter = 5000, run = 500, 
         popSize = 10, maxiter = 1)

i = c(1,5); GA@population[i,]
set.seed(1)
out1 = gaperm_cxCrossover_R(GA, i)
set.seed(1)
out2 = gaperm_cxCrossover_Rcpp(GA, i)
identical(out1, out2)
out1$children - out2$children
 
microbenchmark(gaperm_cxCrossover_R(GA, i),
               gaperm_cxCrossover_Rcpp(GA, i),
               unit = "relative")
*/

// [[Rcpp::export]]
List gaperm_pmxCrossover_Rcpp(RObject object, IntegerVector parents)
{
  IntegerMatrix pop = object.slot("population");
  int           n = pop.ncol();
  IntegerVector seq = Rcpp::seq(0, n-1);
  IntegerVector cxPoints = Rcpp::sample(seq, 2, false);
                cxPoints = Rcpp::seq(min(cxPoints), max(cxPoints));
  // IntegerVector seqdiff = setdiff(seq, cxPoints);
  IntegerVector seqdiff = setdiff_asR(seq, cxPoints);
  IntegerMatrix parentsPop(2,n);
                parentsPop(0,_) = pop(parents[0]-1,_);
                parentsPop(1,_) = pop(parents[1]-1,_);
  NumericVector fitness(2, NA_REAL);
  IntegerMatrix children(2, n);
  std::fill(children.begin(), children.end(), NumericVector::get_na());
  for(int j=0; j < cxPoints.length(); j++)
  {
    children(0,cxPoints[j]) = parentsPop(0,cxPoints[j]);
    children(1,cxPoints[j]) = parentsPop(1,cxPoints[j]);
  }
  IntegerVector ch(n);
  IntegerVector pp(n);
  for(int j=0; j < seqdiff.length(); j++)
  {
    ch = children(0,_); ch = ch[cxPoints];
    if(is_false(any(parentsPop(1,seqdiff(j)) == ch)))
      children(0,seqdiff(j)) = parentsPop(1,seqdiff(j));
    ch = children(1,_); ch = ch[cxPoints];
    if(is_false(any(parentsPop(0,seqdiff(j)) == ch)))
      children(1,seqdiff(j)) = parentsPop(0,seqdiff(j));
  }
  ch = children(0,_);
  pp = parentsPop(1,_);
  ch[is_na(ch)] = setdiff_asR(pp, ch);
  children(0,_) = ch;
  ch = children(1,_);
  pp = parentsPop(0,_);
  ch[is_na(ch)] = setdiff_asR(pp, ch);
  children(1,_) = ch;
  
  List out = List::create(Rcpp::Named("children") = children,
                          Rcpp::Named("fitness")  = fitness);
  return out;  
} 
 
/***
library(GA)
library(microbenchmark)

data(eurodist)
D <- as.matrix(eurodist)
tourLength <- function(tour, distMatrix)
{ 
  tour <- c(tour, tour[1])
  route <- embed(tour, 2)[,2:1]
  sum(distMatrix[route])
}
tspFitness <- function(tour, ...) { 1/tourLength(tour, ...) }
        
GA <- ga(type = "permutation", 
         fitness = tspFitness, distMatrix = D,
         lower = 1, upper = attr(eurodist, "Size"), 
         # popSize = 50, maxiter = 5000, run = 500, 
         popSize = 10, maxiter = 1)

i = c(1,5); GA@population[i,]
set.seed(1)
out1 = gaperm_pmxCrossover_R(GA, i)
set.seed(1)
out2 = gaperm_pmxCrossover_Rcpp(GA, i)
identical(out1, out2)
out1$children - out2$children

microbenchmark(gaperm_pmxCrossover_R(GA, i),
               gaperm_pmxCrossover_Rcpp(GA, i),
               unit = "relative")
*/

// [[Rcpp::export]]
List gaperm_oxCrossover_Rcpp(RObject object, IntegerVector parents)
{
  IntegerMatrix pop = object.slot("population");
  int           n = pop.ncol();
  IntegerVector seq = Rcpp::seq(1, n-2);
  IntegerVector cxPoints = Rcpp::sample(seq, 2, false);
                cxPoints = Rcpp::seq(min(cxPoints), max(cxPoints));
  IntegerMatrix parentsPop(2,n);
                parentsPop(0,_) = pop(parents[0]-1,_);
                parentsPop(1,_) = pop(parents[1]-1,_);
  NumericVector fitness(2, NA_REAL);
  IntegerMatrix children(2, n);
  std::fill(children.begin(), children.end(), NumericVector::get_na());
  for(int j=0; j < cxPoints.length(); j++)
  {
    children(0,cxPoints[j]) = parentsPop(0,cxPoints[j]);
    children(1,cxPoints[j]) = parentsPop(1,cxPoints[j]);
  }
  for(int j=0; j < 2; j++)
  { 
    IntegerVector pos = c_int(Rcpp::seq(max(cxPoints)+1, n-1), 
                              Rcpp::seq(0, max(cxPoints)));
    IntegerVector ch = children(j,_); 
                  ch = ch[cxPoints];
    IntegerVector po = parentsPop(abs(j-1),_); 
                  po = po[pos];
    IntegerVector val = setdiff_asR(po, ch);
    ch = children(j,_); 
    IntegerVector ival = intersect_asR(pos, which_asR(is_na(ch)));
    ch[ival] = val;
    children(j,_) = ch;
  }

  List out = List::create(Rcpp::Named("children") = children,
                          Rcpp::Named("fitness")  = fitness);
  return out;  
} 
 
/***
library(GA)
library(microbenchmark)

data(eurodist)
D <- as.matrix(eurodist)
tourLength <- function(tour, distMatrix)
{ 
  tour <- c(tour, tour[1])
  route <- embed(tour, 2)[,2:1]
  sum(distMatrix[route])
}
tspFitness <- function(tour, ...) { 1/tourLength(tour, ...) }
        
GA <- ga(type = "permutation", 
         fitness = tspFitness, distMatrix = D,
         lower = 1, upper = attr(eurodist, "Size"), 
         # popSize = 50, maxiter = 5000, run = 500, 
         popSize = 10, maxiter = 1)

i = c(1,5); GA@population[i,]
set.seed(1)
out1 = gaperm_oxCrossover_R(GA, i)
set.seed(1)
out2 = gaperm_oxCrossover_Rcpp(GA, i)
out1$children - out2$children

microbenchmark(gaperm_oxCrossover_R(GA, i),
               gaperm_oxCrossover_Rcpp(GA, i),
               unit = "relative")
*/

// [[Rcpp::export]]
List gaperm_pbxCrossover_Rcpp(RObject object, IntegerVector parents)
{
  IntegerMatrix pop = object.slot("population");
  int           n = pop.ncol();
  IntegerVector seq = Rcpp::seq(0, n-1);
  IntegerVector cxPoints = Rcpp::sample(seq, n, true);
                cxPoints = Rcpp::unique(cxPoints);
  IntegerMatrix parentsPop(2,n);
                parentsPop(0,_) = pop(parents[0]-1,_);
                parentsPop(1,_) = pop(parents[1]-1,_);
  NumericVector fitness(2, NA_REAL);
  IntegerMatrix children(2, n);
  std::fill(children.begin(), children.end(), NumericVector::get_na());
  for(int j=0; j < cxPoints.length(); j++)
  {
    children(0,cxPoints[j]) = parentsPop(1,cxPoints[j]);
    children(1,cxPoints[j]) = parentsPop(0,cxPoints[j]);
  }
  for(int j=0; j < 2; j++)
  { 
    IntegerVector ch  = children(j,_);
    IntegerVector pos = which_asR(is_na(ch));
    IntegerVector val = as<IntegerVector>(setdiff_asR(parentsPop(j,_), ch[cxPoints]));
    ch[pos] = val;
    children(j,_) = ch;
  }

  List out = List::create(Rcpp::Named("children") = children,
                          Rcpp::Named("fitness")  = fitness);
  return out;  
} 
 
/***
library(GA)
library(microbenchmark)

data(eurodist)
D <- as.matrix(eurodist)
tourLength <- function(tour, distMatrix)
{ 
  tour <- c(tour, tour[1])
  route <- embed(tour, 2)[,2:1]
  sum(distMatrix[route])
}
tspFitness <- function(tour, ...) { 1/tourLength(tour, ...) }
        
GA <- ga(type = "permutation", 
         fitness = tspFitness, distMatrix = D,
         lower = 1, upper = attr(eurodist, "Size"), 
         # popSize = 50, maxiter = 5000, run = 500, 
         popSize = 10, maxiter = 2)

i = c(1,5); GA@population[i,]
set.seed(1)
out1 = GA:::gaperm_pbxCrossover_R(GA, i)
set.seed(1)
out2 = GA:::gaperm_pbxCrossover_Rcpp(GA, i)
out1$children - out2$children

microbenchmark(GA:::gaperm_pbxCrossover_R(GA, i),
               GA:::gaperm_pbxCrossover_Rcpp(GA, i),
               unit = "relative")
*/

// [[Rcpp::export]]
IntegerVector gaperm_simMutation_Rcpp(RObject object, int parent)
{
  IntegerMatrix pop = object.slot("population");
  int           n = pop.ncol();
  IntegerVector mutate = pop((parent-1),_);
  IntegerVector seq = Rcpp::seq(0, n-1);
  IntegerVector m = Rcpp::sample(seq, 2, false);
                m = Rcpp::seq(min(m), max(m));
  IntegerVector i(n);
  if((min(m) == 0) & (max(m) == n-1)) 
  { 
    i = rev(m);
  } else 
  {
    if(min(m) == 0) 
    { 
      i = c_int(rev(m), Rcpp::seq(max(m)+1, n-1)); 
    } else 
    { 
      if(max(m) == n-1) 
      { 
        i = c_int(Rcpp::seq(0, min(m)-1), rev(m)); 
      } else 
      { 
        i = c_int(Rcpp::seq(0, min(m)-1), rev(m));
        i = c_int(i, Rcpp::seq(max(m)+1, n-1)); 
      }
    }
  }
  mutate = mutate[i];
  return mutate;  
} 

/***
library(GA)
library(microbenchmark)

data(eurodist)
D <- as.matrix(eurodist)
tourLength <- function(tour, distMatrix)
{ 
  tour <- c(tour, tour[1])
  route <- embed(tour, 2)[,2:1]
  sum(distMatrix[route])
}
tspFitness <- function(tour, ...) { 1/tourLength(tour, ...) }
        
GA <- ga(type = "permutation", 
         fitness = tspFitness, distMatrix = D,
         lower = 1, upper = attr(eurodist, "Size"), 
         # popSize = 50, maxiter = 5000, run = 500, 
         popSize = 10, maxiter = 1)

i = 1; GA@population[i,]
set.seed(1)
out1 = gaperm_simMutation_R(GA, i)
set.seed(1)
out2 = gaperm_simMutation_Rcpp(GA, i)
out1 - out2

microbenchmark(gaperm_simMutation_R(GA, i),
               gaperm_simMutation_Rcpp(GA, i),
               unit = "relative")
*/

// [[Rcpp::export]]
IntegerVector gaperm_ismMutation_Rcpp(RObject object, int parent)
{
  IntegerMatrix pop = object.slot("population");
  int           n = pop.ncol();
  IntegerVector mutate = pop((parent-1),_);
  IntegerVector seq = Rcpp::seq(0, n-1);
  IntegerVector m   = Rcpp::sample(seq, 1, false);
                seq = Rcpp::seq(0, n-2);
  IntegerVector pos = Rcpp::sample(seq, 1, false);
                pos = Rcpp::seq(0, pos[0]);
  IntegerVector i = c_int(setdiff_asR(pos, m), m);
                pos = Rcpp::seq(max(pos)+1, n-1);
                i = c_int(i, setdiff_asR(pos, m));
  mutate = mutate[i];
  return mutate;  
} 

/***
library(GA)
library(microbenchmark)

data(eurodist)
D <- as.matrix(eurodist)
tourLength <- function(tour, distMatrix)
{ 
  tour <- c(tour, tour[1])
  route <- embed(tour, 2)[,2:1]
  sum(distMatrix[route])
}
tspFitness <- function(tour, ...) { 1/tourLength(tour, ...) }
        
GA <- ga(type = "permutation", 
         fitness = tspFitness, distMatrix = D,
         lower = 1, upper = attr(eurodist, "Size"), 
         # popSize = 50, maxiter = 5000, run = 500, 
         popSize = 10, maxiter = 1)

i = 1; GA@population[i,]
set.seed(1)
out1 = gaperm_ismMutation_R(GA, i)
set.seed(1)
out2 = gaperm_ismMutation_Rcpp(GA, i)
out1 - out2

microbenchmark(gaperm_simMutation_R(GA, i),
               gaperm_simMutation_Rcpp(GA, i),
               unit = "relative")
*/

// [[Rcpp::export]]
IntegerVector gaperm_swMutation_Rcpp(RObject object, int parent)
{
  IntegerMatrix pop = object.slot("population");
  int           n = pop.ncol();
  IntegerVector parentPop = pop((parent-1),_);
  IntegerVector mutate    = pop((parent-1),_);
  IntegerVector seq = Rcpp::seq(0, n-1);
  IntegerVector m   = Rcpp::sample(seq, 2, false);
  mutate[m[0]] = parentPop[m[1]];
  mutate[m[1]] = parentPop[m[0]];
  return mutate;  
} 

/***
library(Rcpp)
sourceCpp("/Users/luca/R/GA/misc/genope.cpp")
library(GA)
library(microbenchmark)

data(eurodist)
D <- as.matrix(eurodist)
tourLength <- function(tour, distMatrix)
{ 
  tour <- c(tour, tour[1])
  route <- embed(tour, 2)[,2:1]
  sum(distMatrix[route])
}
tspFitness <- function(tour, ...) { 1/tourLength(tour, ...) }
        
GA <- ga(type = "permutation", 
         fitness = tspFitness, distMatrix = D,
         lower = 1, upper = attr(eurodist, "Size"), 
         # popSize = 50, maxiter = 5000, run = 500, 
         popSize = 10, maxiter = 1)

i = 1; GA@population[i,]
set.seed(1)
out1 = gaperm_swMutation_R(GA, i)
set.seed(1)
out2 = gaperm_swMutation_Rcpp(GA, i)
out1 - out2

microbenchmark(gaperm_swMutation_R(GA, i),
               gaperm_swMutation_Rcpp(GA, i),
               unit = "relative")
*/

// [[Rcpp::export]]
IntegerVector gaperm_dmMutation_Rcpp(RObject object, int parent)
{
  IntegerMatrix pop = object.slot("population");
  int           n = pop.ncol();
  IntegerVector parentPop = pop((parent-1),_);
  IntegerVector seq = Rcpp::seq(1, n);
  IntegerVector m = Rcpp::sample(seq, 2, false);
                m = Rcpp::seq(min(m), max(m));
  int           l = max(m)-min(m)+1;
                seq = Rcpp::seq(1, std::max(1,n-l));
  int           pos = as<int>(Rcpp::sample(seq, 1, false));
                seq = setdiff_asR(Rcpp::seq(1,n), m);
  IntegerVector mutate(n);
  if(seq.length() == 0) {
    mutate = parentPop;
  } else {
    IntegerVector i = c_int(seq[Rcpp::seq(0,pos-1)], m);
    if(pos < seq.length())
      i = c_int(i, seq[Rcpp::seq(pos,seq.length()-1)]);
    mutate = parentPop[i-1];
  }
  return mutate;  
} 

/***
library(GA)
library(microbenchmark)

data(eurodist)
D <- as.matrix(eurodist)
tourLength <- function(tour, distMatrix)
{ 
  tour <- c(tour, tour[1])
  route <- embed(tour, 2)[,2:1]
  sum(distMatrix[route])
}
tspFitness <- function(tour, ...) { 1/tourLength(tour, ...) }
        
GA <- ga(type = "permutation", 
         fitness = tspFitness, distMatrix = D,
         lower = 1, upper = attr(eurodist, "Size"), 
         # popSize = 50, maxiter = 5000, run = 500, 
         popSize = 10, maxiter = 1)

i = 3; GA@population[i,]
set.seed(5)
out1 = gaperm_dmMutation_R(GA, i)
set.seed(5)
out2 = gaperm_dmMutation_Rcpp(GA, i)
out1 - out2

microbenchmark(gaperm_dmMutation_R(GA, i),
               gaperm_dmMutation_Rcpp(GA, i),
               unit = "relative")
*/

// [[Rcpp::export]]
IntegerVector gaperm_scrMutation_Rcpp(RObject object, int parent)
{
  IntegerMatrix pop = object.slot("population");
  int           n = pop.ncol();
  IntegerVector parentPop = pop((parent-1),_);
  IntegerVector seq = Rcpp::seq(1, n);
  IntegerVector m = Rcpp::sample(seq, 2, false);
                m = Rcpp::seq(min(m), max(m));
                m = Rcpp::sample(m, m.length(), false);
  IntegerVector i(n);
                i = setdiff_asR(Rcpp::seq(1, min(m)), m);
                i = c_int(i, m);
                i = c_int(i, setdiff_asR(Rcpp::seq(max(m), n), m));
  IntegerVector mutate = parentPop[i-1];
  return mutate;  
} 

/***
library(GA)
library(microbenchmark)

data(eurodist)
D <- as.matrix(eurodist)
tourLength <- function(tour, distMatrix)
{ 
  tour <- c(tour, tour[1])
  route <- embed(tour, 2)[,2:1]
  sum(distMatrix[route])
}
tspFitness <- function(tour, ...) { 1/tourLength(tour, ...) }
        
GA <- ga(type = "permutation", 
         fitness = tspFitness, distMatrix = D,
         lower = 1, upper = attr(eurodist, "Size"), 
         # popSize = 50, maxiter = 5000, run = 500, 
         popSize = 10, maxiter = 1)

i = 5; GA@population[i,]
set.seed(3)
out1 = gaperm_scrMutation_R(GA, i)
set.seed(3)
out2 = gaperm_scrMutation_Rcpp(GA, i)
out1 - out2

microbenchmark(gaperm_scrMutation_R(GA, i),
               gaperm_scrMutation_Rcpp(GA, i),
               unit = "relative")
*/

//
// Miscellaneous
//

// [[Rcpp::export]]
double ga_pmutation_Rcpp(RObject object, double p0 = NA_REAL, 
                                         double p = NA_REAL, 
                                         double T = NA_REAL)
{
  double        maxiter = object.slot("maxiter");
  double        iter = object.slot("iter");
  if(std::isnan(p0))  p0 = 0.5;
  if(std::isnan(p))   p = 0.01;
  if(std::isnan(T))   T = round_double(maxiter/2,0);
  // linear decay ----------------------
  // double pmutation = p;
  // if(iter <= T) pmutation = p0 - (p0-p)/T * (iter-1)
  // exponential decay -----------------
  double        pmutation = (p0 - p)*exp(-2*(iter-1)/T) + p;
  return pmutation;  
}  

/***
library(GA)
library(microbenchmark)
f <- function(x) { 20 + x[1]^2 + x[2]^2 - 10*(cos(2*pi*x[1]) + cos(2*pi*x[2])) }
GA <- ga(type = "real-valued", fitness = f, lower = c(-5.12,-5.12), upper = c(5.12,5.12), maxiter = 10, seed = 1)
GA@maxiter <- 1000
identical(ga_pmutation_R(GA), ga_pmutation_Rcpp(GA))

microbenchmark(ga_pmutation_R(GA), 
               ga_pmutation_Rcpp(GA),
               unit = "relative")
*/

// [[Rcpp::export]]
NumericVector optimProbsel_Rcpp(NumericVector x, double q = NA_REAL)
{
  if(std::isnan(q))  q = 0.25;
  double        eps = sqrt(std::numeric_limits<double>::epsilon());
  // selection pressure parameter
                q = std::min(std::max(eps, q), 1.0-eps);
  NumericVector r = as<NumericVector>(rank_asR(x, true));
  NumericVector prob = exp(log(q) + (r-1.0)*log(1.0-q));
  prob[!is_finite(prob)] = 0.0;
  prob = prob/sum(prob);
  return prob;  
}  

/***
library(GA)
x <- runif(100, 1, 10)
x[sample(1:100,10)] <- NA
prob1 <- GA:::optimProbsel_R(x, q = 0.1)
prob2 <- GA:::optimProbsel_Rcpp(x, q = 0.1)
qcc::describe(abs(prob1 - prob2))
plot(x, prob1)
points(x, prob2, col = 4, pch = 0) 
 
microbenchmark(GA:::optimProbsel_R(x),
               GA:::optimProbsel_Rcpp(x),
               unit = "relative")
*/

//
//  DIFFERENTIAL EVOLUTION OPERATORS
//

// [[Rcpp::export]]
List gareal_de_Rcpp(RObject object, Function fitness, double F = NA_REAL, double p = NA_REAL)
{
  NumericMatrix pop     = object.slot("population");
  NumericMatrix newpop  = clone(pop);
  NumericVector f       = object.slot("fitness");
  NumericVector newf    = clone(f);
  int           popSize = pop.nrow();
  IntegerVector popseq  = Rcpp::seq(0, popSize-1);
  int           n       = pop.ncol();
  IntegerVector nseq    = Rcpp::seq(0, n-1); // seq_len(n);
  NumericVector lb      = object.slot("lower");
  NumericVector ub      = object.slot("upper");
  // if(std::isnan(F)) 
  //   F = 0.8;
  // else                
  //   F = std::max(0.0, std::min(F, 2.0));
  if(!std::isnan(F)) 
    F = std::max(0.0, std::min(F, 2.0));
  if(std::isnan(p)) 
    p = 0.5;
  else
    p = std::max(0.0, std::min(p, 1.0));
  
  NumericVector x(n);
  IntegerVector r(3);
  NumericVector v(n);
  int J; double fx; double Fi;
  for(int i=0; i < popSize; i++)
  {
    r = Rcpp::sample(popseq, 3, false);
    if(std::isnan(F))
      Fi = R::runif(0.5,1.0);
    else
      Fi = F;
    v = pop(r[0],_) + Fi*(pop(r[1],_) - pop(r[2],_));
    J = Rcpp::sample(nseq, 1, false)[0];
    x = pop(i,_);
    for(int j=0; j < n; j++)
    {
      double u = R::runif(0,1);
      if( (u < p) || (j == J) ) 
        x[j] = v[j];
      // reset to random amount if outside the bounds
      if(x[j] < lb[j])
        x[j] = (lb[j] + R::runif(0.0,1.0)*(ub[j] - lb[j]));
      if(x[j] > ub[j])
        x[j] = (ub[j] - R::runif(0.0,1.0)*(ub[j] - lb[j]));
    }
    fx = as<double>(fitness(x));
    // update pop and fitness if improved solution
    if(fx > f[i])
    {
      newf[i] = fx;
      newpop(i,_) = x;
    }
  }

  List out = List::create(Rcpp::Named("population") = newpop,
                          Rcpp::Named("fitness") = newf);
  return out;  
}

/***

library(GA)
library(microbenchmark)
Rastrigin <- function(x1, x2) 20 + x1^2 + x2^2 - 10*(cos(2*pi*x1) + cos(2*pi*x2))
fitness <- function(x) -Rastrigin(x[1], x[2])
obj <- de(fitness = fitness, lower = c(-5.12, -5.12), upper = c(5.12, 5.12), 
          popSize = 10, maxiter = 1, seed = 1)

set.seed(123)
out1 = GA:::gareal_de_R(obj, fitness, F = 0.8, p = 0.5)
set.seed(123)
out2 = GA:::gareal_de_Rcpp(obj, fitness, F = 0.8, p = 0.5)
identical(out1, out2)

microbenchmark(GA:::gareal_de_R(GA, fitness, F = 0.8, p = 0.5),
               GA:::gareal_de_Rcpp(GA, fitness, F = 0.8, p = 0.5),
               unit = "relative")
*/
