setwd("G:/charles-code/code_project/Depmap")
getwd()
library(tidyr)

library(openxlsx)
library(tibble)
library(pheatmap)
library(ggplot2)
file <- read.csv("gene_effect_corrected.csv",header = TRUE,check.names = FALSE)
file_t <-read.csv("gene_effect_corrected_2.csv",header = TRUE,row.names = "X")
RNAi <- read.csv("D2_combined_gene_dep_scores.csv",header = TRUE,row.names = "X")

rownames(file)<-file[,1]
file <- file[-1]

'file2 <- as.tibble(file)
file3<-remove_rownames(file2)
has_rownames(file3)
column_to_rownames(as.data.frame(file3),var = "X")'

file_t <- t(file)
a <- rownames(file_t[1:6,])
b <- as.data.frame(strsplit(a,"[()]"))[1,]
d <- strsplit(a,"[()]")
d <- unlist(d)



c <- as.data.frame(strsplit(a,"[()]"))[2,]
rownames(file_t)<-b

write.csv(file_t,"gene_effect_corrected_2.csv")

a<-file_t[1:6,1:6]
c<- rownames(a)

b <- t(as.data.frame(strsplit(c,"[()]")))[,1]

d <- t(as.data.frame(strsplit(c,"[()]")))[,3]


hist(file_t[1,])

ggplot(file_t, aes())+geom_density(fill="pink")










