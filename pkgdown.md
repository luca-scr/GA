# How to create a package website using `pkgdown`

<https://lbusettspatialr.blogspot.it/2017/08/building-website-with-pkgdown-short.html>

## Prerequisites

```{r}
require(devtools)
use_readme_rmd()
use_news_md()
use_github_links()
use_cran_badge() 
```

## Getting Started

```{r}
devtools::install_github("hadley/pkgdown")
library(pkgdown)
# build_site()
pkgdown:::build_site_rstudio()
```

## Customizing appearence and structure of the website

Crate and customize a text file named `_pkgdown.yaml` in the root folder of your project.

## Deploying the website to Git Hub and the web

- Commit and push your changes to the remote;
- Login on GitHub, navigate to your repo and select "Settings".
- Scroll down to find the "GitHub pages" subsection and, under "Source", select "master branch/docs folder" (from this, you can see that it is fundamental that your website material is pushed to the master branch) and look for the url provided above "master branch/docs folder".
- Click on the link of your website and... congratulations: you just published your new pkgdown website !
