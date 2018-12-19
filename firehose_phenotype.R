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

b <- as.data.frame(t(all_tumor_countdata[-1]))
colnames(b)<-all_tumor_countdata$geneID

pheno <- read.csv(choose.files(),sep = '\t', header = TRUE)
c<-as.data.frame(t(pheno[-1]))
colnames(c) <- pheno$Hybridization.REF
c$sample <-toupper(colnames(pheno[-1]))

name_tumor <- colnames(all_tumor_countdata[-1])
a <- substr(name_tumor,1,12)
#a <-  c("GeneID",a)
b$sample <- a
b <- b[!duplicated(b[,20532]),]

finally <- inner_join(c,b,by =c("sample"="sample"))
rownames(finally)<-finally$sample
finally3 <- filter(finally,days_to_last_followup ==0)

finally <- as.data.frame(t(finally))
colnames(finally3)
finally3 <- as.data.frame(t(finally3))

finally2 <- subset(finally,selected = setdiff(finally$sample,finally3$sample))
write.csv(finally,"all_tumor_coadread.csv",
           row.names = TRUE)































