#Rtsne
library("openxlsx")
library("ggplot2")
library("Rtsne")
library("pcaMethods")
library("ggcorrplot")
library("ggfortify")
library("scatterplot3d")
train <- read.xlsx(file.choose())
train2 <- read.xlsx(file.choose(),2)

a <- as.data.frame(train[,3:8]) #LJJ RAN-seq data
a <- as.data.frame(train[,5:274]) #all array data
a <- as.data.frame(train[,3:433]) #coadread data from TCGA

arrage <- train2[,10]

#correlation
pca_correlation <- cor(a)
pca_correlation <- as.data.frame(pca_correlation)
ggcorrplot(pca_correlation,
           hc.order = TRUE,
           lab = TRUE)+
  scale_fill_gradientn(limits=c(0.93,1),
                       colors = c("blue","white","red"))

#princomp
p2 <- princomp(a)
p2
p2$scores
p2$loadings
p2$center

screeplot(p2,type = "line")                             #plot
abline(1,0)

biplot(p2)
plot(p2$loadings[,1],
     p2$loadings[,2],
     col="red",
     pch=20,
     xlab = "PC1",
     ylab = "PC2")
plot(p2$scores[,1],
     p2$scores[,2],
     col="red",
     pch=20,
     xlab = "PC1",
     ylab = "PC2")

pca_plot <- as.data.frame(p2$loadings[,1:2])                    #ggplot
pca_plot <- as.data.frame(p2$scores[,1:2])

lab_x_1 <- paste("PC1","(",round((summary(p2))$importance[2,1]*100,1),"%)",sep="")
lab_y_1 <- paste("PC1","(",round((summary(p2))$importance[2,2]*100,1),"%)",sep="")
group_1 <- rownames(p2$loadings)
group_1 <- factor(c(rep("normal",16),
                  rep("crypt",5),
                  rep("colon_cell_line",44),
                  rep("breast",9),
                  rep("B stage",41),
                  rep("C stage",25),
                  rep("D stage",50),
                  rep("Met stage",39),
                  rep("live",4),
                  rep("xeno",37)))

ggplot(pca_plot,aes(pca_plot[,1],pca_plot[,2],colour=group_1))+
  scale_color_manual(values = rainbow(n=10))+
  stat_ellipse(level = 0.9)+
  geom_point(size=4)+
  ylab(lab_y_1)+
  xlab(lab_x_1)+
  labs(title="PCA plot")

#prcomp
p4 <- prcomp(a,scale. = T,rank. = 8,retx = T)
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
lab_y <- paste("PC1","(",round((summary(p4))$importance[2,2]*100,1),"%)",sep="")

group <- rownames(p4$rotation)
group <- factor(c(rep("SW480CBSsh2",2),
                  rep("SW480CBSsh3",2),
                  rep("SW480CBSCTRL",2)))

group <- factor(c(rep("normal",16),
                  rep("crypt",5),
                  rep("colon_cell_line",44),
                  rep("breast",9),
                  rep("B stage",41),
                  rep("C stage",25),
                  rep("D stage",50),
                  rep("Met stage",39),
                  rep("live",4),
                  rep("xeno",37)))

group <- factor(arrage)

ggplot(prc_plot,aes(p4$rotation[,1],p4$rotation[,2],colour=group))+
  geom_point(size=3)+
# scale_colour_manual(values=c("blue","white","red","yellow","black","green"))+
  scale_color_manual(values = rainbow(n=13))+
  stat_ellipse(level = 0.9)+
  ylab(lab_y)+
  xlab(lab_x)+
  labs(title="DNA-binding transcription factor PCA plot")

#autoplot
autoplot(p4,frame=TRUE)

#scatterplotsd
scatterplot3d(p4$rotation[,1],p4$rotation[,2],p4$rotation[,3])

#tsne
tsne <- Rtsne(a[,1:100], dims = 2, perplexity = 30, verbose = TRUE, max_iter = 1000) 












