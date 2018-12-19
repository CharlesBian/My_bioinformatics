library(openxlsx)
library(readr)
library(survival)
library(survminer)
library(dplyr)

setwd(choose.dir())
getwd()

all_tumor_countdata <- read.csv(choose.files())
pair_normal_countdata <- read.csv(choose.files())
pair_tumor_countdata <- read.csv(choose.files())

countdata <- all_normal_countdata
countdata <- all_tumor_countdata
countdata <- pair_normal_countdata
countdata <- pair_tumor_countdata

b <- as.data.frame(t(countdata[-1]))
colnames(b)<-countdata$geneID

pheno <- read.csv(choose.files(),sep = '\t', header = TRUE)
c<-as.data.frame(t(pheno[-1]))
colnames(c) <- pheno$Hybridization.REF
c$sample <-toupper(colnames(pheno[-1]))

name_tumor <- colnames(countdata[-1])
a <- substr(name_tumor,1,12)
#a <-  c("GeneID",a)
b$sample <- a
b <- b[!duplicated(b[,20532]),]

finally <- inner_join(c,b,by =c("sample"="sample"))
rownames(finally)<-finally$sample
finally3 <- filter(finally,days_to_last_followup ==0)
rownames(finally3)<-finally3$sample

finally <- as.data.frame(t(finally))
finally3 <- as.data.frame(t(finally3))

finally2 <- subset(finally,select = setdiff(colnames(finally),colnames(finally3)))

write.csv(finally2[-20,],"pair_tumor_coadread.csv",
           row.names = TRUE)


#survival analysis




















































