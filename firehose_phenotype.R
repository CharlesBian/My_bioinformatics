library(openxlsx)
library(readr)
library(survival)
library(survminer)
library(dplyr)

setwd(choose.dir())
getwd()

all_normal_countdata <- read.csv(choose.files())
all_tumor_countdata <- read.csv(choose.files())
pair_normal_countdata <- read.csv(choose.files())
pair_tumor_countdata <- read.csv(choose.files())

b <- t(all_tumor_countdata[-1])
colnames(b)<-all_tumor_countdata$geneID

pheno <- read.csv(choose.files(),sep = '\t', header = TRUE)
c<-t(pheno[-1])
colnames(c) <- pheno$Hybridization.REF


name_tumor <- colnames(all_tumor_countdata[-1])
a <- substr(name_tumor,1,12)
a <-  c("GeneID",a)


name <- colnames(pheno[-1])

finally <- full_join()

































