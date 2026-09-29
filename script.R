
library(devtools)

install.packages("glmnet", dep=TRUE)
install.packages("adegenet")

library(glmnet)
library(adegenet)

if (!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")
remotes::install_github("thibautjombart/adegenet")

list.files(system.file(package = "adegenet"))
url <- "https://adegenet.r-forge.r-project.org/files/Glasgow2015/simGWAS.Rdata"
download.file(url, destfile = "simGWAS.Rdata", mode = "wb")
load("simGWAS.Rdata")
ls()

class(simGWAS)
names(simGWAS)

print(object.size(simGWAS$snps), unit="Mb")
class(simGWAS$snps) #genetic variables, matrix of Single Nucleotide Polymorphism (SNPs)
dim(simGWAS$snps)

class(simGWAS$phen)
table(simGWAS$phen) #the phenotype. R vs S (antibiotic resistance)

#simplify use
snps <- simGWAS$snps
phen <- factor(simGWAS$phen)


#PCA-genetice diversity. all pcas
pca1 <- dudi.pca(snps, scale = FALSE, scannf = FALSE, nf = 5)

graphics.off()
devAskNewPage(FALSE)

pdf("all_plots.pdf")

pca1
pca1$eig #amount of variance per pca
barplot(pca1$eig, main = "eigenvalues")

head(pca1$li)

s.label(pca1$li, sub="PCA - PC 1 and 2")

#distance matrix 
D <- dist(pca1$li[,1:3])^2
clust <- hclust(D, method="complete")

temp <- as.data.frame(as.matrix(D))
temp <- t(as.matrix(D))
temp <- temp[,ncol(temp):1]
par(mar=c(1,5,5,1))
image(x = 1:95, y = 1:95, temp,
      col = rev(heat.colors(nlevels(as.factor(D)))),
      xaxt = "n", yaxt = "n",
      xlab = "", ylab = "")
axis(side = 2, at = 1:95, lab = rownames(snps), las = 2, cex.axis = .46)
axis(side = 3, at = 1:95, lab = rownames(snps), las = 2, cex.axis = .46)
title("Genetic distances between isolates", outer = TRUE, line = -1)


#
plot(clust, main="Clustering (complete linkage) first 5 pcs")

#table of major clusters
pop <- factor(cutree(clust, k=5))
table(pop)

#adding major clusters to the pca
s.class(pca1$li, fac=pop, col=transp(funky(5)), cpoint=2,
        sub="PCA 1-5")
?s.class

#plotting pcas
s.class(pca1$li, fac=pop, xax = 3, yax = 4, col=transp(funky(5)), cpoint=2,
        sub="PCA - axes 3 and 4")

s.class(pca1$li, fac=pop, xax = 1, yax = 2, col=transp(funky(5)), cpoint=2,
        sub="PCA - axes 1 and 2")

#comparing with phenotype
s.class(pca1$li, fac=phen, xax = 1, yax = 2, col=transp(c("royalblue","red")), cpoint=2,
        sub="PCA - axes 1 and 2")
s.class(pca1$li, fac=phen, xax = 3, yax = 4, col=transp(c("royalblue","red")), cpoint=2,
        sub="PCA - axes 3 and 4")

#checking if phenotype are correlated with genetic variations
chisq.test(table(phen, pop), simulate=TRUE) #too low pvalue indicate uneven spread across populations


#TEST FOR ASSOCIATION
#1)Fishers exact test - multiple correction needed sinc it carries the test for every snp
#2) Bonferroni correction

pval <- apply(snps, 2, function(e)
  fisher.test(table(factor(e, levels=c(0,1)), phen))$p.value)
pval.corrected.bonf <- p.adjust(pval, method="bonferroni")

#manhattan plot
log.pval <- -log10(pval.corrected.bonf)
set.seed(1)
log.pval <- jitter(log.pval, amount=0.5)
plot(log.pval,
     col = transp(azur(5)),
     pch = 19,
     cex = 1.5,
     ylim=c(-0.5, 25),
     main="Fisher's exact test \n(Bonferroni correction)",
     xlab="SNP loci", ylab="Bonferroni-corrected -log10(p-value)")
thresh <- -log10(0.05)
abline(h=thresh, col = "red", lwd=2)

#POPULATION STRATIFICATION CORRECTION
#regressing along pcas to correct snp matrix
snps.corrected <- apply(snps, 2, function(e)
  residuals(lm(e~pca1$li[,1]+pca1$li[,2]+pca1$li[,3]+
                 pca1$li[,4]+pca1$li[,5])))
dim(snps.corrected)

pca2 <- dudi.pca(snps.corrected, scale=FALSE, scannf=FALSE, nf=5)
barplot(pca2$eig, main="PCA2 eigenvalues")

s.class(pca2$li, fac=pop, xax = 1, yax = 2, col=transp(funky(5)), cpoint=2,
        sub="PCA2 - axes 1 and 2")
s.class(pca2$li, fac=pop, xax = 3, yax = 4, col=transp(funky(5)), cpoint=2,
        sub="PCA2 - axes 3 and 4")

#GWAS AFTER POPULATION STRATIFICATION CORRECTION WITH PCA
#from 1,0 to continous snpa variable eg 0.8.. so fisher is no longer suitable
#LASSO multivariate method
set.seed(1)
LASSO <- cv.glmnet(snps, phen,
                   family="binomial",
                   lambda.min.ratio=0.01, alpha=1)

beta <- as.vector(t(coef(LASSO, s="lambda.min"))) #extrcating coefficents

res <- which(beta[-1] !=0)
coefs.LASSO <- beta[-1][res]
names(coefs.LASSO) <- colnames(snps.corrected)[res]

fit <- LASSO$glmnet.fit

plot(fit, xvar = "dev", label = TRUE)
grid()
title("Deviance explained by LASSO coefficients (corrected)", line=3)

#resulting snps are assumed to be real since population is corrected
res <- which(beta[-1] !=0)
length(res)

result <- snpzip(snps.corrected, phen,
                 xval.plot=TRUE, plot=TRUE, loading.plot=TRUE,
                 method="ward")
par(ask=FALSE)
#DAPC and LASSO validation needs to be carried out before and after populaion startification

dev.off()
