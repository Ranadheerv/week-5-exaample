############################################################
##### Segregation distortion test of molecular markers #####
############################################################

## Example 1
## Datasets from Miscanthus genetic mapping study
## Dong et al (2018) Global Change Biology Bioenergy. 10:165-185: https://onlinelibrary.wiley.com/doi/full/10.1111/gcbb.12472 
## Two parents: Miscanthus sinensis x M. sacchariflorus

setwd("~/Dropbox/PSS 8990 QG/1. Lecture Presentation/Learning R/Week 5")

# classify markers into three allele types
Miscanthus = read.csv("Miscanthus dataset.csv", header = TRUE, row.names = 1)
dim(Miscanthus) # 281 indiviuals in rows; 1926 SNP markers in columns

nrow(Miscanthus) # check number of rows using nrow
ncol(Miscanthus) # check number of columns using ncol

# Try chisq.test in R, use the first marker (TP144):
table(Miscanthus[,1]) # 0: 128; 1:148
chisq.test(c(128, 148), p=c(0.5,0.5))$p.value


#Chi-square test for segregation distortion within each allele freq type:
chisq_pval_Miscanthus=c()
for (i in 1:ncol(Miscanthus)){
  homo_ref_sum = sum(Miscanthus[,i]==0, na.rm=T)
  hetero_sum = sum(Miscanthus[,i]==1, na.rm=T)
  observed_geno_freqs = c(c(homo_ref_sum,hetero_sum))
  chisq_pval_Miscanthus[i] = chisq.test(observed_geno_freqs, p=c(0.5,0.5))$p.value
}

plot(chisq_pval_Miscanthus), xlab="Marker", ylab="chisq p Value")
abline(h=0.05, col="red")
hist(chisq_pval_Miscanthus, breaks=20)
table(chisq_pval_Miscanthus < 0.0001) # All 1926 retained for nondistorted 
table(chisq_pval_Miscanthus < 0.001) # 1919 retained for nondistorted, 7 snps distorted
table(chisq_pval_Miscanthus < 0.05) # 1565 retained for nondistorted, 361 distorted 



## Example 2
## Datasets from Rqtl book
## Broman & Sen: A Guide to QTL Mapping with R/qtl
## F2 population with 300 individuals, genotyped at 100 markers
## Three genotype classes for each marker: AA, AB, BB
F2 = read.csv("F2 dataset.csv", header = TRUE, row.names = 1)

#Chi-square test for segregation distortion within each allele freq type:
chisq_pval_F2=c()
for (i in 1:nrow(F2)){
  AA_sum = F2[i,3]
  AB_sum = F2[i,4]
  BB_sum = F2[i,5]
  observed_geno_freqs = c(AA_sum, AB_sum, BB_sum)
  chisq_pval_F2[i] = chisq.test(observed_geno_freqs, p=c(0.25, 0.5, 0.25))$p.value
}

plot(chisq_pval_F2)
hist(chisq_pval_F2, breaks=20)
table(chisq_pval_F2 < 0.0001) # 95 retained for nondistorted, 5 distorted 
table(chisq_pval_F2 < 0.001) # 92 vs 8
table(chisq_pval_F2 < 0.05) # 79 vs 21 

