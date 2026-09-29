# Project 1 :  ROSETTA R/qtl QTL mapping (epiRIL data)

QTL mapping walkthrough in R using the qtl package, with ROSETTA flowering time phenotypes measured in epiRILs.

## Source

- Tutorial video: https://youtu.be/vRgDnjXBnnc ```

Input files: a genotype file (markers, chromosomes, map positions) and a phenotype CSV (phenotypes first, RIL IDs in the last column, missing data as `NA`).

## Citation

Broman KW, Wu H, Sen S, Churchill GA (2003) R/qtl: QTL mapping in experimental crosses. Bioinformatics 19: 889-890.
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
