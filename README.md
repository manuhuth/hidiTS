# Dynamic Factor Modeling
We review the standard exact dynamic factor model (DFM) and its sample inference via maximum likelihood (ML) and principal component analysis (PCA). We survey four information criteria that aim to estimate the number of factors and discuss how the factor space can be identified. Subsequently, we conduct simulation studies to asses the accuracy of the information criteria, the factor estimates, compare the performance of the ML and the PCA estimator and apply PCA to the FRED-QD database. In this context, we focus on how the performance of the estimator is affected by the structure of the sample covariance matrix's eigenvalues. Given large cut-offs in the eigenvalues, we find that the information criteria yield accurate estimates and that the factors can be estimated consistently in large samples. Regarding real-world application, our results indicate that DFMs can serve as a substantial dimension reduction method carrying most of the sample's variance.

## Continous integration with Travis CI to ensure Reproducibility
To ensure reproducibility, we have integrated Travis CI. The build's history can be found here [![Build Status](https://travis-ci.org/HumanCapitalAnalysis/microeconometrics-course-project-manuhuth.svg?branch=master)](https://travis-ci.org/github/manuhuth/hidiTS)

## Installation

Either install it directly from GitHub (later for the professors and students) 
```javascript
#install.packages("devtools") #only once
devtools::install_github("manuhuth/hidiTS")
```

or clone the repository (do changes) and run
```javascript
#install.packages("devtools") #needed only once
#install.packages("roxygen2") #needed only once
devtools::install(dependencies = TRUE) #installs package
library(hidiTS) #loads package
roxygen2::roxygenise() #creates documentation files (.rd) in man folder

pack <- "hidiTS"
path <- find.package(pack)
system(paste(shQuote(file.path(R.home("bin"), "R")),"CMD", "Rd2pdf", shQuote(path))) #creates Vignette
```
