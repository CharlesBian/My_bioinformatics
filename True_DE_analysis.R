library("openxlsx")
library("limma")
library("DESeq2")
library("ggplot2")
library("Cairo")

rna_seq_data <- read.xlsx(file.choose())

#DE
#rna_seq_data<-subset(rna_seq_data,select = c(1,3:8) )

countdata <- as.matrix(floor(rna_seq_data[1:nrow(rna_seq_data),3:ncol(rna_seq_data)]))
condition <- as.factor(c(rep("sh1",2),rep("sh2",2),rep("Ctrl",2)))
col_data <-data.frame(row.names = colnames(countdata),
                      group_list=condition)

dds <- DESeqDataSetFromMatrix(countData = countdata,
                              colData = col_data,
                              design = ~group_list)
rownames(dds)<-rna_seq_data$gene_id
dds2 <- DESeq(dds)
rld <- rlogTransformation(dds2)
exprSet = assay(rld)

boxplot(countdata)
boxplot(exprSet)

resultsNames(dds2)
res <- results(dds2)
res_ordered <- as.data.frame(res[order(res$padj),])


