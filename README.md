# High-Dimensional Time Series
We evaluate methods and conduct simulation studies for high-dimensional time series models. (more sophisticated description) 

## Continous integration with Travis CI to ensure Reproducibility
To ensure reproducibility, we have integrated Travis CI. The build's history can be found here [![Build Status](https://travis-ci.org/HumanCapitalAnalysis/microeconometrics-course-project-manuhuth.svg?branch=master)](https://travis-ci.org/github/manuhuth/hidiTS)

##Installation

Either install it directly from GitHub.
or clone the repository (do changes) and run

```javascript
devtools::install(dependencies = TRUE)
library(hidiTS)
roxygen2::roxygenise()

pack <- "hidiTS"
path <- find.package(pack)
system(paste(shQuote(file.path(R.home("bin"), "R")),"CMD", "Rd2pdf", shQuote(path)))
```

