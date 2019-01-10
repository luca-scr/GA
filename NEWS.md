# GA 3.2 (2019-01)

- Added website using pkgdown.
- Added a section to vignette on integer optimisation.
- Added function `de()` implementing Differential Evolution based on the description in Simon (2013) Evolutionary Optimization Algorithms. 
- Bug fix in C++ code calling `pow()` function.
- Bug fix in `gails()` when suggestions are provided.
- Long outputs in `summary()` function calls are shortened.
- Fix a bug in `gaperm_pbxCrossover()`.

# GA 3.1.1 (2018-05)

- Bug fix in C++ code calling `pow()` function.

# GA 3.1 (2018-05)

- Genetic operators available in C++ using **Rcpp** package.
- Add parameter `"useRcpp"` in `gaControl()` to control if the C++ implementation of genetic operators should be used. By default is set to TRUE.
- Renamed input parameters `min` and `max` to `lower` and `upper`. The old nomenclature is still accepted but deprecated.
- Added `stopParallel()` function to stop a cluster if `parallel = TRUE`, including `registerDoSEQ()` to avoid problems if foreach loop is used after.
- A parallel cluster is automatically stopped unless a registered parallel back end is provided as argument to `parallel` argument in the `ga()` or `gaisl()` function call.
- Bug fix in `gaperm_Population()` when `min` > 1. Thanks to Romero Barata.
- Update "A quick tour of GA" vignette.
- Add `README.Rmd` that is shown in GitHub and CRAN.
- Create a sticker for the **GA** package.
  
# GA 3.0.2 (2016-06)

- In interactive sessions use different monitor functions depending on whether or not is an RStudio session. This is due to different capabilities of the several consoles from which R can be run. 

# GA 3.0.1 (2016-05)

- **doRNG** package moved from Depends to Suggests due to the unavailability of the binary version of the package on Windows. As a consequence, **GA** package can be installed even without **doRNG** being installed, but in this case results are not reproducible by setting the seed if GAs are executed in parallel.    

# GA 3.0 (2016-05)

- Added option to provide a cluster for parallelisation.
- Added GA-hybrid using `optim()` for local search.
- Added GA Islands model in `gaisl()`.
- Rewrite of `gaMonitor()` to clear the previous output before printing the info about the current iteration. Old version is available in `gaMonitor2()`. Same behaviour for `gaislMonitor()` and `gaislMonitor2()`.

# GA 2.3 (2015-07)

- Added documentation vignettes.

# GA 2.2 (2014-10)

- Slight modification to `plot.ga()` method.
- Included link in main ga documentation to genetic operators.
- Added `ga_pmutation()` function to allow GAs having variable mutation probability.
- pdf files (previously included as vignettes) moved to `inst/doc` with corresponding index.html.

# GA 2.1 (2014-05)

- `.printShortMatrix()` is a function to print part of rows/columns of a matrix
- `print.summary.ga()` accept arguments to be passed to `.printShortMatrix()`. This allow to shorten the printed output in case of large dimensions for the matrices containing the suggestions and the final solutions.
- `plot.ga()` now shadowed the area between the max and the median of fitness values at each iteration. Changed to pch from 17 to 1 for means.
- Function `gaSummary()` is embedded in the code, so it cannot be defined by the user.
- Added function `ga_pmutation()` for computing variable mutation probability. 
- Modified main function `ga()` to allow for pmutation to be a function. This enables the use of variable mutation rate.
- Export `plot.ga` class in `NAMESPACE` and used as S4 method.
- Computing summaries at each step is done by function `gaSummary()`, but it may be defined by the user.
- Export `summary.ga` class in `NAMESPACE` and used as S4 method.
- Default for suggestions argument in `ga()` is set to NULL.
- Fixed a bug in `gabin_uCrossover()`.
- Added description of slots in the help page for `ga-class`.
- Bug fix on setting maxfitness.
- Bug fix on passing seed argument when "snow" parallel is used.
- Add explicitly stop clusters if parallel is used.
- Parallel computing improving set up of clusters.

# GA 2.0 (2013-08)

- Option for parallel computing.
- Added argument keepBest to include the best solution at each iteration.
- Bug fix when using suggestions argument.

# GA 1.1 (2013-04)

- Update citation and references to JSS paper.

# GA 1.0 (2012-06)

- First release on CRAN.
