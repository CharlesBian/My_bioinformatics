library(WGCNA)
library(data.table)
library(stringr)
library(openxlsx)

allowWGCNAThreads()
ALLOW_WGCNA_THREADS = 4
memory.limit(size = 20000) #???

#设定阈值范围
powers = c(c(1:10),seq(from = 12,to = 20, by = 2))

nGenes <- length(a$ACH.000004) 

#获得各个阈值下的R方和平均连接度
sft = pickSoftThreshold(t(a),powerVector = powers,verbose = 5)

#作图
sizeGrWindow(9,5)
par(mfrow = c(1,2))
cexl = 0.9

#Scale-free topology fit index as a function of the soft-th
plot(sft$fitIndices[,1],
     -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab = "Soft Threshold (power)",
     ylab = "Scale Free Topology Model Fit,signed R^2",
     type = "n",
     main = paste("Scale independence"))
text(sft$fitIndices[,1],
     -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,
     cex = cexl,
     col = "red")
plot(sft$fitIndices[,1],
     sft$fitIndices[,5],
     xlab = "Soft Threshold(power)",
     ylab = "Mean Connectivity",
     type = "n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1],
     sft$fitIndices[,5],
     labels = powers,
     cex = cexl,
     col = "red")

#获得临近矩阵
softPower <- sft$powerEstimate
adj <- adjacency(t(a),power = softPower)

TOM <- TOMsimilarity(adj)

dissTOM <- 1- TOM

hierTOM <- hclust(as.dist(dissTOM),method = "average")

#检验选定的β值下记忆网络是否逼近scale free
k <- softConnectivity(datE = t(a),power = softPower)
sizeGrWindow(10,5)
par(mfrow = c(1,2))
hist(k)
scaleFreePlot(k,main = "Check Scale free topology\n")

#使用相异度来聚类为gene tree
geneTree <- hclust(as.dist(dissTOM),method = "average")

#Plot the resulting clustering tree
windows()
sizeGrWindow(12,9)
plot(geneTree,
     xlab="",
     sub="",
     main = "Gene Clusting on TOM-based dissimilarity",
     labels = FALSE,
     hang = 0.04)

#使用动态剪切树挖掘模块
minModuleSize = 30

#动态切割树
dynamicMods = cutreeDynamic(dendro = geneTree,
                            distM = dissTOM,
                            deepSplit = 2,
                            pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize)
table(dynamicMods)

#拓扑热图
nSelect = 400

#For reproducibility, we set the random seed
set.seed(10)
select = sample(nGenes, size = nSelect)
selectTOM = dissTOM[select,select]

#There's no simple way of restricting a clustering tree to a subset of genes,so we must re-cluster
selectTree = hclust(as.dist(selectTOM),method = "average")

dynamicColors = labels2colors(dynamicMods)

selectColors = dynamicColors[select]

#Open a graphical window
sizeGrWindow(9,9)

plotDiss <- selectTOM^softPower
diag(plotDiss) = NA

TOMplot(plotDiss,
        selectTree,
        selectColors,
        main = "Network heatmap plot, selected genes")




























