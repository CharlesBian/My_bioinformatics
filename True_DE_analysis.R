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
res$log2FoldChange
res$padj


ggplot(res_ordered,
       aes(x=res_ordered$log2FoldChange,
           y=-log10(res_ordered$padj)))+
  scale_color_manual(values=c("blue","grey","red"))+
  geom_point()+
  xlim(c(-4,4))+
  ylab("-log10 p-value")+
  xlab("log2 fold change")+
  labs(title="Volcano of CBSsh22 vs CBSsh32")+
  geom_vline(xintercept = c(-1.5,1.5),lty=4,col="black",lwd=0.6)+
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.6)+
  theme(legend.position = "right",
        panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(face = "bold",color = "black",size = 12),
        plot.title = element_text(hjust = 0.5))









