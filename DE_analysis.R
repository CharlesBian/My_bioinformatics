library("openxlsx")
library("limma")
library("DESeq2")
library("ggplot2")
library("Cairo")

liu_rnaseq <- read.xlsx("C:/Users/77214/Desktop/seq_information/BGI/cbskd-expression.xlsx",1)

liu_rnaseq2 <- read.xlsx("C:/Users/77214/Desktop/seq_information/BGI/cbskd-expression.xlsx",2)

liu_rnaseq3 <- read.xlsx("C:/Users/77214/Desktop/seq_information/BGI/cbskd-expression.xlsx",3)

liu_rnaseq4 <- read.xlsx("C:/Users/77214/Desktop/seq_information/BGI/cbskd-expression.xlsx",4)

liu_rnaseq5 <- read.xlsx("C:/Users/77214/Desktop/seq_information/BGI/cbskd-expression.xlsx",5)

BGI_liu <- read.xlsx("C:/Users/77214/Desktop/seq_information/BGI/DE_liu.xlsx",1)

LH_data <- read.xlsx("C:/Users/77214/Desktop/seq_information/LH/hep26429-sup-0010-supptable1.xlsx",1)

ang_data <- read.xlsx("C:/Users/77214/Desktop/seq_information/micro/ang.xlsx",1)

#list of survival
gene_list <- read.xlsx(choose.files())
gene_list2 <-subset(gene_list, gene_list$`CTNNB1-up_p.value_survival`<0.05)

#volcano
liu_rnaseq2$threshold21 <- as.factor(liu_rnaseq2$`(CBSsh2/CTRL)p-value`<0.05
                                     &abs(liu_rnaseq2$`log2(CBSsh21/CTRL1)`)>=1.5)

liu_rnaseq2$threshold21 <- as.factor(ifelse(liu_rnaseq2$`(CBSsh2/CTRL)p-value`<0.05
                                            &abs(liu_rnaseq2$`log2(CBSsh21/CTRL1)`)>=1.5,
                                            ifelse(liu_rnaseq2$`log2(CBSsh21/CTRL1)`>1.5,"Up","Down"),"Not"))

liu_rnaseq2$threshold22 <- as.factor(ifelse(liu_rnaseq2$`(CBSsh2/CTRL)p-value`<0.05
                                            &abs(liu_rnaseq2$`log2(CBSsh22/CTRL2)`)>=1.5,
                                            ifelse(liu_rnaseq2$`log2(CBSsh22/CTRL2)`>1.5,"Up","Down"),"Not"))

liu_rnaseq2$threshold2_3 <- as.factor(ifelse(liu_rnaseq2$`(CBSsh2/CBSsh3)p-value`<0.05
                                             &abs(liu_rnaseq2$`log2(CBSsh21/CBSsh31)`)>=1.5,
                                             ifelse(liu_rnaseq2$`log2(CBSsh21/CBSsh31)`>1.5,"Up","Down"),"Not"))

liu_rnaseq2$threshold2_3_2 <- as.factor(ifelse(liu_rnaseq2$`(CBSsh2/CBSsh3)p-value`<0.05
                                             &abs(liu_rnaseq2$`log2(CBSsh22/CBSsh32)`)>=1.5,
                                             ifelse(liu_rnaseq2$`log2(CBSsh22/CBSsh32)`>1.5,"Up","Down"),"Not"))

liu_rnaseq2$threshold31 <- as.factor(ifelse(liu_rnaseq2$`(CBSsh3/CTRL)p-value`<0.05
                                            &abs(liu_rnaseq2$`log2(CBSsh31/CTRL1)`)>=1.5,
                                            ifelse(liu_rnaseq2$`log2(CBSsh31/CTRL1)`>1.5,"Up","Down"),"Not"))

liu_rnaseq2$threshold32 <- as.factor(ifelse(liu_rnaseq2$`(CBSsh3/CTRL)p-value`<0.05
                                            &abs(liu_rnaseq2$`log2(CBSsh32/CTRL2)`)>=1.5,
                                            ifelse(liu_rnaseq2$`log2(CBSsh32/CTRL2)`>1.5,"Up","Down"),"Not"))


LH_data$threshold1 <- as.factor(ifelse(LH_data$`P.value.Act.HSC+HM.vs.Act.HSC`<0.05
                                       &abs(log2(LH_data$`Act.HSC.+HM.vs.Act.HSC.(Fold)`))>=1.5,
                                       ifelse(log2(LH_data$`Act.HSC.+HM.vs.Act.HSC.(Fold)`)>1.5,"Up","Down"),"Not"))

LH_data$threshold2 <- as.factor(ifelse(LH_data$adj.P.Val.Act.HSC.vs.Q<0.05
                                       &abs(log2(LH_data$`Act.HSC.vs.qHSC.(Fold)`))>=1,
                                       ifelse(log2(LH_data$`Act.HSC.vs.qHSC.(Fold)`)>1,"Up","Down"),"Not"))



Cairo(file="volcano_LH_1.png",
      type = "png",
      units = "in",
      bg = "white",
      width = 5.5,
      height = 5,
      pointsize = 12,
      dpi = 300)


ggplot(liu_rnaseq2,
       aes(x=liu_rnaseq2$`log2(CBSsh22/CBSsh32)`,
           y=-log10(liu_rnaseq2$`(CBSsh2/CBSsh3)p-value`),
           colour=threshold2_3_2))+
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


ggplot(LH_data,
       aes(x=log2(LH_data$`Act.HSC.+HM.vs.Act.HSC.(Fold)`),
           y=-log10(LH_data$`P.value.Act.HSC+HM.vs.Act.HSC`),
           colour=threshold1))+
  scale_color_manual(values=c("blue","grey","red"))+
  geom_point()+
  xlim(c(-4,4))+
  ylab("-log10 p-value")+
  xlab("log2 fold change")+
  labs(title="Volcano plot")+
  geom_vline(xintercept = c(-1.5,1.5),lty=4,col="black",lwd=0.6)+
  geom_hline(yintercept = -log10(0.05),lty=4,col="black",lwd=0.6)+
  theme(legend.position = "right",
        panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(face = "bold",color = "black",size = 12),
        plot.title = element_text(hjust = 0.5))

dev.off()


#venn+go
library(VennDiagram)

x<-ang_data$GeneID

a<-subset(LH_data,LH_data$threshold2=="Up")
b<-subset(LH_data,LH_data$threshold2=="Down")
c<- rbind(a,b)

a<-subset(liu_rnaseq2,liu_rnaseq2$threshold21=="Down")
b<-subset(liu_rnaseq2,liu_rnaseq2$threshold21=="Up")
c <- rbind(a,b)

d<-subset(liu_rnaseq2,liu_rnaseq2$threshold32=="Down")
e<-subset(liu_rnaseq2,liu_rnaseq2$threshold2_3_2=="Up")
f <- rbind(d,e)

g<- rbind(c,f)

h <- as.character(intersect(c$Gene.ID,f$Gene.ID))

i <- as.character(intersect(b$Gene.ID,e$Gene.ID))

venn.diagram(list(CBSsh3vs = liu_rnaseq4$`Gene.ID-1`,CBSsh2_2down=liu_rnaseq5$`Gene.ID-2`),
             filename = "1_1.tif",
             col = "transparent",
             fill = c("cornflowerblue","darkorchid1"),
             label.col = "black",
             cat.col = c("blue","red"),
             cat.dist = c(0.03,0.03),
             cat.pos = c(-8,8),
             cat.cex = 1.5,
             cex = 2,
             main = "sh3 downregulate cross",
             main.cex = 2)

venn.diagram(list(home=liu_rnaseq3$Gene.ID,BGI_sh2vsctrl=BGI_liu$Gene.ID),
             filename = "BGI_1.tif",
             col = "transparent",
             fill = c("cornflowerblue","darkorchid1"),
             label.col = "black",
             cat.col = c("blue","red"),
             cat.dist = c(0.03,0.03),
             cat.pos = c(0,0),
             cat.cex = 2,
             cex = 2)

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

#GSEA analysis
gsearesult <- gseGO(sort(BGI_liu$Gene.ID,decreasing = TRUE), 
                    OrgDb = "org.Hs.eg.db",
                    ont = "CC", 
                    nPerm = 100, 
                    minGSSize = 5,
                    pvalueCutoff = 0.1)
gseaplot(gsearesult,gene)











