library(DESeq2)
library(readr)

coadread_rawdata <- read.table(choose.files(),sep = '\t', header = TRUE)
ESCA_rawdata <- read.table(choose.files(),sep = '\t', header = TRUE)
STAD_rawdata <- read.table(choose.files(),sep = '\t', header = TRUE)
STES_rawdata <- read.table(choose.files(),sep = '\t', header = TRUE)

raw_counts <-ESCA_rawdata

#find the normal samples
geneID <- raw_counts[,1]
x<- colnames(raw_counts)
x <- x[-1]
x <- substr(x, 14,15)
x <- as.integer(x)
x <- x %in% (10:19)
y <- colnames(raw_counts)[-1]
y <- y[x]
normal_counts <- cbind(geneID, subset(raw_counts, select = y))

#cancer sample
x <- colnames(normal_counts)[-1]
y <- substr(x,9,12)
z <- subset(raw_counts,
            select = grep(paste(y, collapse = "|"), 
                          colnames(raw_counts),
                          value = TRUE))
tumor_counts <- subset(z, select = setdiff(colnames(z),colnames(normal_counts)))
tumor_counts <- cbind(geneID, tumor_counts)


