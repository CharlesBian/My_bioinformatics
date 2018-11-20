library("openxlsx")
library("limma")
library("DESeq2")
library("edgeR")
library("ggplot2")
library("Cairo")
library("dplyr")

rna_seq_data <- read.xlsx(file.choose())
Genename <- read.xlsx(file.choose())
rna_seq_data <- merge(rna_seq_data,Genename,by.x="gene_id",by.y="ENTREZ_GENE_ID",all.x = TRUE)
sum(is.na(rna_seq_data$symbol))

#DE
#rna_seq_data<-subset(rna_seq_data,select = c(1,3:8) )
colnamber<-ncol(rna_seq_data)
countdata<- as.matrix(floor(rna_seq_data[1:nrow(rna_seq_data),4:colnamber]))
countdata <- subset(countdata,select=c(1,2,5,6))
condition <- as.factor(c(rep("sh1",2),rep("Ctrl",2)))
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

res_ordered$threshold <- as.factor(ifelse(res_ordered$pvalue<0.05&abs(res_ordered$log2FoldChange)>=1,
                                            ifelse(res_ordered$log2FoldChange>1,"Up","Down"),"Not"))

downname <- head(subset(res_ordered,threshold=="Down"),13)  #subset the important gene
downname <- merge(downname,rna_seq_data,by.x=rownames(downname),by.y="gene_id",all.x = TRUE)
name <- subset(rna_seq_data,gene_id==rownames(downname))

sum(is.na(rna_seq_data$symbol))
upname <- head(subset(res_ordered,threshold=="Up"),13)  #subset the important gene
upname <- merge(upname,rna_seq_data,by.x=Row.names,by.y="gene_id",all.x = TRUE)

dename <- subset(res_ordered,pvalue<0.05&threshold!="Not")
write.xlsx(rownames(dename),file.choose())

ggplot(res_ordered,
       aes(x=res_ordered$log2FoldChange,
           y=-log10(res_ordered$pvalue),
           colour = threshold))+
  scale_color_manual(values=c("skyblue","grey","red"))+
  geom_point(size = 2.5)+
  xlim(c(-4,4))+
  ylab("-log10 p-value")+
  xlab("log2 fold change")+
  labs(title="Volcano of CBSsh")+
  geom_vline(xintercept = c(-1,1),lty=4,col="black",lwd=0.6)+
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.6)+
  theme(legend.position = "right",
        panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(face = "bold",color = "black",size = 12),
        plot.title = element_text(hjust = 0.5))+
  annotate(geom ="text",
           x=downname$log2FoldChange,
           y=-log10(downname$pvalue),
           label=rownames(downname),
           size = 2.3)+
  annotate(geom ="text",
           x=upname$log2FoldChange,
           y=-log10(upname$pvalue),
           label=rownames(upname),
           size = 2.3)

##venn plot

library(VennDiagram)

venn1 <- read.xlsx(file.choose())
venn2 <- read.xlsx(file.choose())
venn3 <- read.xlsx(file.choose())

a <- as.character(intersect(venn2$geneid,venn3$Gene.ID))
b <- as.character(union(venn3$Gene.ID,venn2$geneid))


venn.diagram(list(CBSsh2vsctrl =venn1$gene_id,CBSsh23vsctrl=venn2$geneid),
             filename = "1_1.tif",
             col = "transparent",
             fill = c("cornflowerblue","darkorchid1"),
             label.col = "black",
             cat.col = c("blue","red"),
             cat.dist = c(0.03,0.03),
             cat.pos = c(-8,8),
             cat.cex = 1.5,
             cex = 2,
             main = "CBSsh DE genes cross",
             main.cex = 2)

venn.diagram(list(CBSsh2vsctrl_new =venn1$gene_id,CBSsh2vsctrl_old = venn3$Gene.ID),
             filename = "BGI_1.tif",
             col = "transparent",
             fill = c("cornflowerblue","darkorchid1"),
             label.col = "black",
             cat.col = c("blue","red"),
             cat.dist = c(0.03,0.03),
             cat.pos = c(-8,8),
             cat.cex = 1.5,
             cex = 2,
             main = "CBSsh DE genes cross",
             main.cex = 2)

#go
library("DOSE")
library("clusterProfiler")
ego3 <- enrichGO(gene = c$GENE.ID,
                 OrgDb = "org.Mm.eg.db",
                 ont = "CC",
                 pvalueCutoff = 0.05,
                 readable = TRUE)
summary(as.data.frame(ego))

write.csv(summary(ego),"liu_GO.csv",row.names = FALSE)
head(ego3,n = 12L)
barplot(ego3, drop = TRUE, showCategory = 30)
dotplot(ego3,showCategory = 40)

cnetplot(ego3)
emapplot(ego3)

#KEGG analysis
ekk <- enrichKEGG(gene =c$GENE.ID,
                  organism = "mmu",
                  pvalueCutoff = 0.05,
                  pAdjustMethod = "BH",
                  use_internal_data = TRUE)

ekk2 <- enrichKEGG(gene =ang_data$Symbol,
                   organism = "ang",
                   pvalueCutoff = 0.05,
                   pAdjustMethod = "BH",
                   use_internal_data = FALSE)
ekk <- enrichKEGG(gene = gene_list2$gene_id,
                  organism = "hsa",
                  pvalueCutoff = 0.05,
                  pAdjustMethod = "BH",
                  use_internal_data = TRUE)

summary(as.data.frame(ekk))

write.csv(summary(ekk),"KEGG1-enrich-down.csv",row.names = FALSE)

head(ekk,n = 12L)

barplot(ekk, drop = TRUE, showCategory = 40,color = "qvalue")

dotplot(ekk,
        showCategory = 40,
        color = "qvalue",
        title = "kegg pathway")

cnetplot(ekk)
emapplot(ekk)

#pathway view
require("pathview")
pathview(gene.data = c[,2],
         gene.idtype = "entrez",
         pathway.id = ekk,
         species = "mmu",
         # limit = list(gene = max(log2(LH_data$`Act.HSC.vs.qHSC.(Fold)`)),cpd = 1),
         bins = list(gene = 10,cpd=10),
         out.suffix = "LH_pathview_1",
         kegg.native = TRUE)


pathview(gene.data = c$GENE.ID,
         gene.idtype = "entrez",
         pathway.id = "mmu00230",
         species = "mmu",
         out.suffix = "LH_pathview",
         limit = list(gene = max(log2(LH_data$`Act.HSC.vs.qHSC.(Fold)`)),cpd = 1),
         bins = list(gene=10),
         kegg.native = FALSE,
         expand.node = FALSE,
         split.group = TRUE,
         map.null = TRUE,
         sign.pos = "bottomleft")













