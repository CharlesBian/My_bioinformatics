#1--open excel
install.packages("openxlsx")
library("openxlsx")

#2--GO analysis and kegg analysis
install.packages("BioInstaller",repos="https://bioconductor.org/packages/3.7/bioc")
source("http://bioconductor.org/biocLite.R")
biocLite()
biocLite("limma")
biocLite("DOSE")
biocLite("clusterProfiler")
biocLite("pathview")
biocLite("RCytoscape")
biocLite("GOstats")
biocLite("org.Hs.eg.db")
biocLite("org.Mm.eg.db")
biocLite("recount")
biocLite(c("GenomicFeatures", "AnnotationDbi"))
biocLite("Category")
biocLite("humanid")
biocLite("RSQLite")
biocLite("GEOquery")
biocLite("annaffy")

library("GSEABase")
library("GOstats")
library("DOSE")
library("clusterProfiler")
library("RCytoscape")
library("Category")
library("limma")
library("org.Mm.eg.db")
library("org.Hs.eg.db")

#DE analysis
biocLite("DESeq2")
biocLite("edgeR")


#ms
install.packages("XML")
biocLite("xcms")  # ms
biocLite("faahKO")
biocLite("MassSpecWavelet")
biocLite("msdata")
biocLite("mzR")
biocLite("ncdf")
biocLite("Rcpp")
biocLite("waveslim")
biocLite("MSnbase")
biocLite("multtest")
biocLite("CAMERA")
biocLite("RNetCDF")



#3--cairo plot
install.packages("Cairo")
library("Cairo")

#4--survival packages
install.packages("survival")
library("survival")

#5--correlation
install.packages("corrplot")
library("corrplot")

install.packages("corrgram")
library("corrgram")

library("corrplot")
library("lattice")
library("ggplot2")
library("gplots")

#6--venn 
install.packages("devtools")
library("devtools")
devtools::install_github("ramnathv/rCharts")
devtools::install_github("neuhausi/canvasXpress") #用于venn图
library(canvasXpress)

install.packages("samr")
library("samr")

install.packages("VennDiagram")  #venn图
library("grid")
library("futile.logger")
library(VennDiagram)

install.packages("colorfulVennPlot")
library(colorfulVennPlot)

#7-- heatmap plot
install.packages("pheatmap")  #heatmap
library(pheatmap)


#8--machine learning
#t-sne
install.packages("Rtsne")
install.packages("ggcorrplot")
install.packages("ggfortify")
install.packages("scatterplot3d")
install.packages("Seurat")


#neural networks
install.packages("nnet")

#Recursive Partitioning
install.packages("rpart")
install.packages("party")

#Random Forests
install.packages("randomForest")

#Support Vector Machines and Kernel Methods
install.packages("e1071")

#Regularized
install.packages("glmnet")

#Model selection
install.packages("ROCR")

#caret,mlr
install.packages("caret")
install.packages("mlr")

#rattle+RGtk2
install.packages("RGtk2")
install.packages("rattle")
library(RGtk2)
library(rattle)
rattle()


#update R
install.packages("installr")
install.packages("stringi")
install.packages("stringr")
library("stringi")
library("stringr")
library(installr)
updateR()



.packages(all.available = TRUE)







