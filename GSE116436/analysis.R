#   Differential expression analysis with limma
library("GEOquery")
library("annaffy")
library("limma")
library("affy")
library("Biobase")

# load series and platform data from GEO

gset <- getGEO("GSE116436", GSEMatrix =TRUE, AnnotGPL=TRUE)
if (length(gset) > 1) idx <- grep("GPL571", attr(gset, "names")) else idx <- 1
gset <- gset[[idx]]
