library(DESeq2)
library(openxlsx)
library(readr)
library(ggplot2)
library(dplyr)

coadread_rawdata <- read.table(choose.files(),sep = '\t', header = TRUE)
ESCA_rawdata <- read.table(choose.files(),sep = '\t', header = TRUE)
STAD_rawdata <- read.table(choose.files(),sep = '\t', header = TRUE)
STES_rawdata <- read.table(choose.files(),sep = '\t', header = TRUE)

raw_counts <-coadread_rawdata
raw_counts <-ESCA_rawdata
raw_counts <-STAD_rawdata
raw_counts <-STES_rawdata

#find the normal samples and the tumor samples
geneID <- raw_counts[,1]
x<- colnames(raw_counts)
x <- x[-1]
x <- substr(x, 14,15)
x <- as.integer(x)
x <- x %in% (10:19)
y <- colnames(raw_counts)[-1]
y <- y[x]
all_normal_counts <- cbind(geneID, subset(raw_counts, select = y))
all_tumor_counts <- subset(raw_counts,select = setdiff(colnames(raw_counts),colnames(all_normal_counts)))
all_tumor_counts <- cbind(geneID,all_tumor_counts)

#find the pair sample
a <- colnames(all_normal_counts)[-1]
a <- substr(a,9,12)
b <- subset(all_tumor_counts,
            select = sort(grep(paste(a, collapse = "|"), 
                          colnames(all_tumor_counts),
                          value = TRUE)))
  #pair_tumor_counts <- subset(z, select = setdiff(colnames(z),colnames(all_normal_counts)))
pair_tumor_counts <- cbind(geneID, b)


c <- colnames(pair_tumor_counts)[-1]
c <- substr(c,9,12)
d <- subset(all_normal_counts,
            select = grep(paste(c, collapse = "|"),
                          colnames(all_normal_counts),
                          value = TRUE))
  #pair_normal_counts <- subset(z, select = setdiff(colnames(z),colnames(pair_tumor_counts)))
pair_normal_counts <- cbind(geneID, d)
all <- cbind(pair_normal_counts,pair_tumor_counts[-1])


#write the files
setwd(choose.dir())
getwd()
write.csv(pair_normal_counts,"pair_coadread_normal_rawcounts.csv")
write.csv(pair_tumor_counts,"pair_coadread_tumor_rawcounts.csv")
write.csv(all_normal_counts,"all_coadread_normal_rawcounts.csv")
write.csv(all_tumor_counts,"all_coadread_tumor_rawcounts.csv")

#pheatmap
library(pheatmap)
pheatmap(all[1:20531,2:ncol(all)],
         scale = "column",
         fontsize = 9,
         color = colorRampPalette(c("royalblue2","white","firebrick1"))(10),
         #cluster_cols = FALSE,
         show_rownames = FALSE)


#prcomp

setwd(choose.dir())
getwd()
a<-read.csv(file.choose())
b<-read.csv(file.choose())
c<-read.csv(file.choose())
d<-read.csv(file.choose())
e<-read.csv(file.choose())
f<-read.csv(file.choose())
g<-read.csv(file.choose())
h<-read.csv(file.choose())  
all<-cbind(a,b[-1],c[-1],d[-1],e[-1],f[-1])
all<-cbind(b, d[-1], f[-1])

pair_normal_counts <- read.csv(file.choose())
pair_tumor_counts <- read.csv(file.choose())
all <- cbind(pair_normal_counts,pair_tumor_counts[-1])

p4 <-prcomp(all[-1],scale. = T,rank. = 8,retx = T)
p4
p4$x
p4$rotation
p4$sdev
p4$center

plot(p4$rotation[,1],                 #plot
     p4$rotation[,2],
     col="red",
     pch=20,
     xlab = "PC1",
     ylab = "PC2")
biplot(p4)

prc_plot <- as.data.frame(p4$rotation[,1:2])                     #ggplot

lab_x <- paste("PC1","(",round((summary(p4))$importance[2,1]*100,1),"%)",sep="")
lab_y <- paste("PC2","(",round((summary(p4))$importance[2,2]*100,1),"%)",sep="")
lab_z <- paste("PC3","(",round((summary(p4))$importance[2,3]*100,1),"%)",sep="")

group <- rownames(p4$rotation)
group <- factor(c(rep("coadread normal",50),
                  rep("coadread tumor",50),
                  rep("ESCA normal",11),
                  rep("ESCA tumor",11),
                  rep("STAD normal",32),
                  rep("STAD tumor",32)))
group <- factor(c(rep("coadread tumor",50),
                  rep("ESCA tumor",11),
                  rep("STAD tumor",32)))

ggplot(prc_plot,aes(p4$rotation[,1],p4$rotation[,2],colour=group))+
  geom_point(size=3)+
  # scale_colour_manual(values=c("blue","white","red","yellow","black","green"))+
  scale_color_manual(values = rainbow(n=6))+
  stat_ellipse(level = 0.92)+
  ylab(lab_y)+
  xlab(lab_x)+
  labs(title="normal vs tumor")

#plot 3D PCA
library(rgl)
library("RColorBrewer")


plot3d(p4$rotation[,1],
       p4$rotation[,2],
       p4$rotation[,3],
       xlab = lab_x,
       ylab = lab_y,
       zlab = lab_z,
       col =c(rep("red",50),
              rep("green",11),
              rep("blue",32)))

plot3d(p4$rotation[,1],
       p4$rotation[,2],
       p4$rotation[,3],
       xlab = lab_x,
       ylab = lab_y,
       zlab = lab_z,
       col =c(rep("red",50),
              rep("yellow",50),
              rep("green",11),
              rep("deepskyblue",11),
              rep("blue",32),
              rep("hotpink",32)))

movie3d(spin3d(axis = c(0,0,1),rpm=3),duration = 10, fps = 50,dir = "F:/project/bian_data/GDC/picture")

#use the cmd for magick
system("F:/project/bian_data/GDC/deal/already/pair_data/3d/magick convert -delay 1 *.png ESCA.gif")

ellips <- ellipse3d(cov(cbind(p4$rotation[,1],p4$rotation[,2],p4$rotation[,3])),
                    centre = c(mean(p4$rotation[,1]),mean(p4$rotation[,2]),mean(p4$rotation[,3])),level=0.99)

plot3d(ellips,col="blue",alpha=0.2,add=TRUE,box=FALSE)

#DESeq2
library(vioplot)
setwd(choose.dir())
getwd()

pair_normal_rawcounts<-read.csv(file.choose())
pair_tumor_rawcounts<-read.csv(file.choose())
geneID <-as.vector(pair_normal_rawcounts[,1])
i=1
genenames<-c()
geneid <-c()
for(i in 1:20531)
{
  genenames[i] <- strsplit(geneID[i],'[|]')[[1]][1]
  geneid[i] <- strsplit(geneID[i],'[|]')[[1]][2]
  i=i+1
}

rna_seq_data <- cbind(pair_normal_rawcounts,pair_tumor_rawcounts[-1])
rna_seq_data <-cbind(rna_seq_data[-1],geneid)

colnamber<-ncol(rna_seq_data)-1
countdata<- as.matrix(floor(rna_seq_data[1:nrow(rna_seq_data),1:colnamber]))
countdata <- subset(countdata,select=c(1,2,5,6))
condition <- as.factor(c(rep("normal",32),rep("tumor",32)))
col_data <-data.frame(row.names = colnames(countdata),
                      group_list=condition)

dds <- DESeqDataSetFromMatrix(countData = countdata,
                              colData = col_data,
                              design = ~group_list)

rownames(dds)<-geneid
dds2 <- DESeq(dds)

#<50 samples
rld <- rlogTransformation(dds2)
exprSet = assay(rld)

#>50 samples
vsd<-vst(dds2) 
exprSet = assay(vsd)

boxplot(countdata)
boxplot(exprSet)

resultsNames(dds2)
res <- results(dds2)
res_ordered <- as.data.frame(res[order(res$padj),])
res$log2FoldChange
res$padj

res_ordered$threshold <- as.factor(ifelse(res_ordered$pvalue<0.05&abs(res_ordered$log2FoldChange)>=2,
                                          ifelse(res_ordered$log2FoldChange>2,"Up","Down"),"Not"))

downname <- head(subset(res_ordered,threshold=="Down"),13)  #subset the important gene
downname <- merge(downname,rna_seq_data,by.x=rownames(downname),by.y="geneid",all.x = TRUE)
name <- subset(rna_seq_data,gene_id==rownames(downname))

sum(is.na(rna_seq_data$symbol))
upname <- head(subset(res_ordered,threshold=="Up"),13)  #subset the important gene
upname <- merge(upname,rna_seq_data,by.x=Row.names,by.y="geneid",all.x = TRUE)

dename <- subset(res_ordered,pvalue<0.05&threshold!="Not")
write.xlsx(dename,"STAD_DEgene.xlsx")

ggplot(res_ordered,
       aes(x=res_ordered$log2FoldChange,
           y=-log10(res_ordered$padj),
           colour = threshold))+
  scale_color_manual(values=c("skyblue","grey","red"))+
  geom_point(size = 2.5)+
  xlim(c(-4,4))+
  ylab("-log10 p-value")+
  xlab("log2 fold change")+
  labs(title="Volcano of TCGA_coadread")+
  geom_vline(xintercept = c(-2,2),lty=4,col="black",lwd=0.6)+
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.6)+
  theme(legend.position = "right",
        panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(face = "bold",color = "black",size = 12),
        plot.title = element_text(hjust = 0.5))


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


  
#GO and KEGG
library("DOSE")
library("clusterProfiler")

dename <- read.xlsx(choose.files())
up_gene<-subset(dename,threshold=="Up")
down_gene<-subset(dename,threshold=="Down")

gene <-dename

ego <- enrichGO(gene = rownames(gene),
                 OrgDb = "org.Hs.eg.db",
                 ont = "ALL",
                 pvalueCutoff = 0.05,
                 readable = TRUE)
summary(as.data.frame(ego))

head(ego,n = 12L)
barplot(ego, drop = TRUE, showCategory = 30)
dotplot(ego,
        showCategory = 40,
        color = "qvalue",
        title="Gene ontology of TCGA_coadread")
  
cnetplot(ego)
emapplot(ego)
plotGOgraph(ego)

#KEGG analysis
ekk <- enrichKEGG(gene =rownames(gene),
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
        title = "kegg pathway of TCGA_coadread_upregulate")
  
cnetplot(ekk)
emapplot(ekk)

#GSEA
genelist<-subset(dename,select=log2FoldChange)
genelist <-cbind(geneID=rownames(genelist),genelist)
#write.xlsx(genelist,"coadread_DEgene.xlsx")
genelist <- read.xlsx(choose.files())
c5 <-read.gmt(choose.files())

gl <- genelist[,2]
names(gl) <- as.character((genelist[,1]))
glist <- sort(gl,decreasing = TRUE)

#gsea analysis
gs<-GSEA(glist,
         TERM2GENE = c5,
         verbose = FALSE,
         pvalueCutoff = 0.1)
head(gs)
dotplot(gs)

gseaplot(gs,
         gs[1]$ID,
         title = paste("GO difference:",gs[1]$Description))

#gsea GO analysis
gsea.go <- gseGO(glist,
                 OrgDb= org.Hs.eg.db,
                 pvalueCutoff = 0.5)
summary(as.data.frame(gsea.go))
head(gsea.go)

dotplot(gsea.go,
        showCategory = 40,
        title = "GSEA.GO of TCGA_coadread_upregulate")
gseaplot(gsea.go,
         gsea.go[1]$ID,
         title = paste("Gene ontology:",gsea.go[1]$Description))


#gsea KEGG pathways
gsea.kegg <- gseKEGG(glist,
                     organism = "hsa",
                     pvalueCutoff = 0.5)
summary(as.data.frame(gsea.kegg))
head(gsea.kegg,40)

dotplot(gsea.kegg,
        showCategory = 50,
        title = "GSEA.KEGG of TCGA_coadread")

gseaplot(gsea.kegg,
         gsea.kegg[1]$ID,
         title = paste("KEGG pathways:",gsea.kegg[1]$Description))


#pathway view
require("pathview")
pathview(gene.data = glist,
         pathway.id = 'hsa04310',
         species = "hsa",
         limit = list(gene = max(abs(glist)),cpd=1))

  
pathview(gene.data = glist,
         pathway.id = 'hsa04390',
         species = "hsa",
         limit = list(gene = max(abs(glist)),cpd=1),
         kegg.native = FALSE,
         expand.node = FALSE,
         split.group = TRUE,
         map.null = TRUE,
         sign.pos = "bottomleft")


