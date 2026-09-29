# Project 2 : adegenet GWAS practical (simulated data)

Working through the Day 4 GWAS practical from the adegenet Glasgow 2015 course (Thibaut Jombart), using the simulated dataset `simGWAS`.

## Source
- Practical: https://adegenet.r-forge.r-project.org/files/Glasgow2015/practical-GWAS_day4.pdf
- Data: https://adegenet.r-forge.r-project.org/files/Glasgow2015/simGWAS.Rdata


## Setup
```r
url <- "https://adegenet.r-forge.r-project.org/files/Glasgow2015/simGWAS.Rdata"
download.file(url, destfile = "simGWAS.Rdata", mode = "wb")
load("simGWAS.Rdata")
```

## Citation
Jombart T (2008) adegenet: a R package for the multivariate analysis of genetic markers. Bioinformatics 24: 1403-1405.
