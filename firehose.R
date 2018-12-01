library(DESeq2)
library(readr)

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
a <- sort(substr(a,9,12))
b <- subset(all_tumor_counts,
            select = sort(grep(paste(a, collapse = "|"), 
                          colnames(all_tumor_counts),
                          value = TRUE)))
  #pair_tumor_counts <- subset(z, select = setdiff(colnames(z),colnames(all_normal_counts)))
pair_tumor_counts <- cbind(geneID, b)


c <- colnames(pair_tumor_counts)[-1]
c <- sort(substr(c,9,12))
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
write.csv(pair_normal_counts,"pair_ESCA_normal_counts.csv")
write.csv(pair_tumor_counts,"pair__ESCA_tumor_counts.csv")
write.csv(all_normal_counts,"all__ESCA_normal_counts.csv")
write.csv(all_tumor_counts,"all__ESCA_tumor_counts.csv")

#pheatmap
library(pheatmap)
pheatmap(all[1:20531,2:ncol(all)],
         scale = "column",
         fontsize = 9,
         color = colorRampPalette(c("royalblue2","white","firebrick1"))(10),
         #cluster_cols = FALSE,
         show_rownames = FALSE)


#prcomp
library(ggplot2)

pair_normal_counts <- read.csv(file.choose())
pair_tumor_counts <- read.csv(file.choose())
all <- cbind(pair_normal_counts,pair_tumor_counts[-1])

p4 <- prcomp(all[-1],scale. = T,rank. = 8,retx = T)
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
group <- factor(c(rep("normal",50),
                  rep("tumor",50)))

ggplot(prc_plot,aes(p4$rotation[,1],p4$rotation[,2],colour=group))+
  geom_point(size=3)+
  # scale_colour_manual(values=c("blue","white","red","yellow","black","green"))+
  scale_color_manual(values = rainbow(n=13))+
  stat_ellipse(level = 0.92)+
  ylab(lab_y)+
  xlab(lab_x)+
  labs(title="ESCA normal vs tumor")

library(rgl)
library("RColorBrewer")


plot3d(p4$rotation[,1],
       p4$rotation[,2],
       p4$rotation[,3],
       xlab = lab_x,
       ylab = lab_y,
       zlab = lab_z,
       col =c(rep("red",50),rep("orange",50)))
movie3d(spin3d(axis = c(0,0,1),rpm=3),duration = 10, fps = 50,dir = "F:/project/bian_data/GDC/picture")


#use the cmd for magick
system("F:/project/bian_data/GDC/deal/already/pair_data/3d/magick convert -delay 1 *.png ESCA.gif")

ellips <- ellipse3d(cov(cbind(p4$rotation[,1],p4$rotation[,2],p4$rotation[,3])),
                    centre = c(mean(p4$rotation[,1]),mean(p4$rotation[,2]),mean(p4$rotation[,3])),level=0.99)

plot3d(ellips,col="blue",alpha=0.2,add=TRUE,box=FALSE)





