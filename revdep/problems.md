# Setup

## Platform

|setting  |value                        |
|:--------|:----------------------------|
|version  |R version 3.3.0 (2016-05-03) |
|system   |x86_64, darwin13.4.0         |
|ui       |RStudio (0.99.1166)          |
|language |(EN)                         |
|collate  |en_US.UTF-8                  |
|tz       |Europe/Rome                  |
|date     |2016-05-10                   |

## Packages

|package |*  |version |date       |source              |
|:-------|:--|:-------|:----------|:-------------------|
|GA      |   |3.0     |2016-05-10 |local (luca-scr/GA) |

# Check results
3 packages with problems

## datafsm (0.1.0)
Maintainer: Nay John J. <john.j.nay@gmail.com>  
Bug reports: https://github.com/JohnNay/datafsm/issues

2 errors | 1 warning  | 0 notes

```
checking examples ... ERROR
Running examples in ‘datafsm-Ex.R’ failed
The error most likely occurred in:

> base::assign(".ptime", proc.time(), pos = "CheckExEnv")
> ### Name: evolve_model
> ### Title: Use a Genetic Algorithm to Estimate a Finite-state Machine Model
> ### Aliases: evolve_model
> 
> ### ** Examples
> 
> # Create data:
> cdata <- data.frame(period = rep(1:10, 1000),
+                    outcome = rep(1:2, 5000),
+                    my.decision1 = sample(1:0, 10000, TRUE),
+                    other.decision1 = sample(1:0, 10000, TRUE))
> (res <- evolve_model(cdata, cv=FALSE))
We set a seed for you to make this reproducible. It is 12. If you want the same results, next time you run this with the same settings, also set the seed argument of this function to 12.
Error in fitness(Pop[i, ], ...) : unused argument (maxfitness = 1)
Calls: evolve_model -> <Anonymous> -> fitness
Execution halted

checking tests ... ERROR
Running the tests in ‘tests/testthat.R’ failed.
Last 13 lines of output:
  Looking up history[4], period 1: (0, 536870925, 0)
  Looking up history[1], period 1: (0, 536870925, 0)
  Looking up history[2], period 0: (1, 536870925, 0)
  Looking up history[3], period 1: (0, 536870925, 0)
  Looking up history[4], period 1: (0, 536870925, 0)
  Looking up history[1], period 1: (1, 536870925, 0)
  testthat results ================================================================
  OK: 3 SKIPPED: 0 FAILED: 2
  1. Error: evolve_model() returns correct type of object (@test_mainfunc.R#7) 
  2. Error: evolve_model() returns warnings and errors (@test_mainfunc.R#14) 
  
  Error: testthat unit tests failed
  Execution halted

checking re-building of vignette outputs ... WARNING
Error in re-building vignettes:
  ...
Quitting from lines 105-106 (datafsmVignette.Rmd) 
Error: processing vignette 'datafsmVignette.Rmd' failed with diagnostics:
unused argument (maxfitness = 1)
Execution halted

```

## GAabbreviate (1.2)
Maintainer: Baljinder K. Sahdra <baljinder.sahdra@acu.edu.au>

1 error  | 0 warnings | 0 notes

```
checking examples ... ERROR
Running examples in ‘GAabbreviate-Ex.R’ failed
The error most likely occurred in:

> base::assign(".ptime", proc.time(), pos = "CheckExEnv")
> ### Name: GAabbreviate
> ### Title: Abbreviating items (from questionnaire or other) measures using
> ###   Genetic Algorithms (GAs)
> ### Aliases: GAabbreviate
> ### Keywords: optimize multivariate survey
... 9 lines ...
+                nrow = nsubject, ncol = nitems)
> scales = cbind(rowSums(items[,1:10]), rowSums(items[,11:15]))
> 
> GAA = GAabbreviate(items, scales, itemCost = 0.01, maxItems = 5, 
+                    popSize = 50, maxiter = 300, run = 100,
+                    verbose = TRUE)
> plot(GAA)
Error in image.default(y = (0:iter), x = 1:l, z = t(object$results$solution),  : 
  increasing 'x' and 'y' values expected
Calls: plot ... plot -> plot.GAabbreviate -> image -> image.default
Execution halted
```

## mcga (3.0)
Maintainer: Mehmet Hakan Satman <mhsatman@istanbul.edu.tr>

1 error  | 0 warnings | 0 notes

```
checking examples ... ERROR
Running examples in ‘mcga-Ex.R’ failed
The error most likely occurred in:

> base::assign(".ptime", proc.time(), pos = "CheckExEnv")
> ### Name: mcga2
> ### Title: Performs a machine-coded genetic algorithm search for a given
> ###   optimization problem
> ### Aliases: mcga2
> 
> ### ** Examples
> 
> f <- function(x){ 
+   return(-sum( (x-5)^2 ) )
+ }
> myga <- mcga2(fitness = f, popSize = 100, maxiter = 300, 
+               min = rep(-50,5), max = rep(50,5))
Error in fitness(Pop[i, ], ...) : unused argument (maxfitness = Inf)
Calls: mcga2 -> ga -> fitness
Execution halted
```

