##Set working directory!

##install packages

install.packages("lme4")
install.packages("qtl")
install.packages("devtools")
if (!requireNamespace("qtlTools", quietly = TRUE)) {
  if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
  remotes::install_github("jtlovell/qtlTools")
}
library(qtlTools)

##load libraries
library(lme4)
library(qtl)
library(qtlDesign)
library(qtlTools)

##Open a PDF for the plots
pdf("qtl_results.pdf")

#import the genotype data with the phenotype data

qtlfile <- read.cross (
  format = "csvs",
  dir = "~/Desktop/QTL-basics/qtl files",
  genfile = "Cortijo-epigenotypes-QTL.csv",
  phefile = "qtlfile_pheno.csv",
  crosstype = "riself" )

# removing outliers for first thne 2nd phenotype

qtlfile$pheno[c(which(
  qtlfile$pheno[,1] > 
      (mean(qtlfile$pheno[,1], na.rm = TRUE) + 3*(sd(qtlfile$pheno[,1], na.rm = TRUE))))
  ,
  
  which(
    qtlfile$pheno[,1] < 
      (mean(qtlfile$pheno[,1], na.rm = TRUE) - 3*(sd(qtlfile$pheno[,1], na.rm = TRUE))))),1] <- NA

qtlfile$pheno[c(which(
  qtlfile$pheno[,2] > 
    (mean(qtlfile$pheno[,2], na.rm = TRUE) + 3*(sd(qtlfile$pheno[,2], na.rm = TRUE))))
  ,
  
  which(
    qtlfile$pheno[,2] < 
      (mean(qtlfile$pheno[,2], na.rm = TRUE) - 3*(sd(qtlfile$pheno[,2], na.rm = TRUE))))),2] <- NA

##Uses the hidden Markov model model to calculate the probabilities 
##of the true underlying genotypes given the observed multipoint marker data, 
##with possible allowance for genotyping errors.
qtlfile <- calc.genoprob(qtlfile, step=2, error.prob = 0.001)

#phenotype names.

phenames(qtlfile)

phes <- phenames(qtlfile)[1:2]

##Perform the QTL analysis and significance threshold

s1<-scanone(qtlfile, pheno.col = phes, method = "hk")
modelperms <-scanone(qtlfile, pheno.col = phes, method = "hk", n.perm = 1000)

#plot the results.One color per phenotype.
graphics.off()
cols <- c("blue", "red") 
plot(s1, type="n", ylim = c(0,max(as.matrix(s1[,-c(1:2)]))), 
     ylab = "LOD Score", main = "all phenotypes plotted together")
for(i in 1:length(phes)) plot(s1, add=T, lodcolumn = i, col = cols[i])

#OR
cols <- c("blue", "red")
plot(s1, lodcolumn = 1:length(phes), col = cols, ylab = "LOD Score",
  main = "all phenotypes plotted together")

##Add the permutation thresholds. For every phenotype, add a line below.
abline(h = summary(modelperms)[1,1], col = "blue", lty = "dotted", lwd = 2)
abline(h = summary(modelperms)[1,2], col = "red", lty = "dotted", lwd = 2)

##Look at the significant QTLs
print(pullSigQTL(qtlfile, pheno.col=phes,
                 s1.output = s1,  perm.output = modelperms, 
                 returnQTLModel=FALSE, alpha = 0.05,
                 controlAcrossCol = TRUE))

##Create an R object with the significant QTLs
##This is to build the results file with the QTL intervals.
mods <- pullSigQTL(qtlfile, pheno.col=phes,
                 s1.output = s1,  perm.output = modelperms, 
                 returnQTLModel=TRUE, alpha = 0.05,
                 controlAcrossCol = TRUE)

mods<-mods[sapply(mods,length)!=1] # a NULL QTL model has length = 1, toss these

##Note: The qtlStats() function is using "lodint" to calculate the confidence intervals.
##You need to add the argument "lodint = FALSE" to use bayesint instead. 
stats<-lapply(names(mods), function(x){
  mod<-mods[[x]]
  form<-paste("y ~",paste(mod$altname, collapse = " + "))
  mod<-refineqtl(qtlfile, qtl = mod, formula = form, 
                 verbose=F, method="hk", pheno.col = x)
  qtlStats(qtlfile, 
           pheno.col = x, 
           form = form, 
           mod = mod)
})
stats <- do.call(rbind, stats)
print(head(stats))

##Write to file
write.table(stats, "stats.csv", sep = ",", row.names = FALSE)

segmentsOnMap(cross=qtlfile, phe=stats$phenotype,
              col = c("red", "blue"), ##add more colors if more phenotypes
              chr=stats$chr, 
              l = stats$lowposition, 
              lwd=5,
              chrBuffer = c(0.2, 0.2),
              h = stats$highposition,
              legendPosition = "bottomleft",
              legendCex= 1,
              leg.lwd = 5,
              leg.inset=.05)









