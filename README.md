[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/GA)](https://cran.r-project.org/package=GA)
[![](http://cranlogs.r-pkg.org/badges/GA)](https://cran.r-project.org/package=GA)
[![](http://cranlogs.r-pkg.org/badges/grand-total/GA)](https://cran.r-project.org/package=GA)


# GA

An R package for optimization using genetic algorithms. The package provides a flexible general-purpose set of tools for implementing genetic algorithms search in both the continuous and discrete case, whether constrained or not. Users can easily define their own objective function depending on the problem at hand. Several genetic operators are available and can be combined to explore the best settings for the current task. Furthermore, users can define new genetic operators and easily evaluate their performances. Local search using general-purpose optimisation algorithms can be applied stochastically to exploit interesting regions. GAs can be run sequentially or in parallel, using an explicit master-slave parallelisation or a coarse-grain islands approach.


## Installation

Get the released version from CRAN:

```R
install.packages("GA")
```

Or the development version from github:

```R
# install.packages("devtools")
devtools::install_github("luca/GA")
```

## How to Use This Package

See the papers in the references section below. 
A quick intro vignette is also available, which can be accessed using

```R
vignette("GA")
```

## References:

Scrucca, L. (2013) GA: A Package for Genetic Algorithms in R. **Journal of Statistical Software**, 53(4), 1-37. URL http://www.jstatsoft.org/v53/i04/.

Scrucca, L. (2016) On some extensions to GA package: hybrid optimisation, parallelisation and islands evolution. Submitted to **R Journal**. Pre-print available at http://arxiv.org/abs/1605.01931.


