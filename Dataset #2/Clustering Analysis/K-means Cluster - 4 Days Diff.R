#install.packages("readxl")
#install.packages("factoextra")
#install.packages("viridis")
#install.packages("rafalib")
#source("https://bioconductor.org/biocLite.R")
#biocLite("limma")
#biocLite("gplots")
library(rafalib)
library(readxl)
library("RColorBrewer")
library(limma)
library("gplots")
library(factoextra)
library(dplyr)
library(tidyverse)
library(viridis)
library(DESeq2)
library(genefilter)

setwd("/Users/ap14958/OneDrive - University of Bristol/Genomics Facility Bioinformatics/Project #1 LuciaCellSignalling 21_2_19/Dataset #2/")
full_table = data.matrix(read.csv("full_table_Marucci_curated_v2.csv", header = TRUE, row.names = 1, stringsAsFactors = FALSE))

metadata = read.table("MetaData_Exp2.txt", sep = "\t", header = TRUE)
metadata


timepoint = "4 Days"
metadata_t=metadata[metadata$TimePoint==timepoint,]
ColourStrategy = metadata_t$Sample
ColourStrategy2 = metadata_t$Condition


order1 = as.character(metadata_t$Original_Sample_Name)
iso_table = full_table[, match(order1, colnames(full_table))]
iso_table[,1:ncol(iso_table)] <- sapply(round(iso_table[,1:ncol(iso_table)], 0), as.integer)
head(iso_table)

####################################
#### Heatmap Generation
####################################

vst = varianceStabilizingTransformation(iso_table, blind = TRUE)
select <- vst[order(rowMeans(vst),decreasing=T)[1:100],]

cols <- palette(magma(256))[as.fumeric(as.character(ColourStrategy2))]

heatmap.2(select, col=viridis(256), Colv=TRUE,
          scale="row", key=T, keysize=1, symkey=T,
          density.info="none", trace="none",
          cexCol=1.3, labRow=F, labCol = ColourStrategy, ColSideColors=cols,
          main="100 Top Expressed Genes Heatmap")
dev.copy(png, width = 1000,height = 1000, paste0(timepoint, "-Heatmap.png"))
dev.off()

#library(RColorBrewer) 
#hmcol <- colorRampPalette(brewer.pal(9, "GnBu"))(100)
#head(cbind(colnames(iso_table),cols))

####################################
#### Hierarchical Cluster Generation
####################################

clusters= hclust(dist(t(iso_table),method = "euclidean"), method = "complete")
plot(clusters, cex = 0.6, hang = -1)
#rect.hclust(clusters, k = 6, border = 2:5)


myplclust(clusters, labels=metadata_t$Treatment, lab.col=cols, cex=0.6, hang = 0.06, main = paste0(timepoint, " - Cluster Dendrogram"))
#rect.hclust(clusters, k = 6, border = 2:5)

dev.copy(png, width = 1000,height = 1000, paste0(timepoint, "-cluster_samples.png"))
dev.off()

#Optimise level to cluster
ClusterLevel = 5000
abline(h=ClusterLevel)
hclusters <- cutree(clusters, h=ClusterLevel)
table(true=metadata_t$Treatment, cluster=hclusters)

hclusters <- cutree(clusters, k=3)
table(true=metadata_t$Treatment, cluster=hclusters)

distance <- get_dist(t(iso_table))
fviz_dist(distance, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))


#select <- iso_table[order(rowMeans(iso_table),decreasing=T)[1:10],]
#clusters= hclust(dist(select,method = "euclidean"), method = "complete")
#plot(clusters, cex = 0.6, hang = -1)
#rect.hclust(clusters, k = 6, border = 2:5)
#dev.copy(png, width = 1000,height = 1000, paste0(timepoint, "-cluster_genes.png"))
#dev.off()

#distance <- get_dist(iso_table)
#fviz_dist(distance, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))


####################################
#### K-Means Cluster Generation
####################################

k = 2

iso_table2 = t(iso_table)
iso_table3=iso_table2[,which(colSums(iso_table2)>1)]
fviz_nbclust(iso_table3, kmeans, method = "gap_stat")
dev.copy(png, width = 1000,height = 1000, paste0(timepoint, "-K-means_optimisation.png"))
dev.off()

kmeanscluster = kmeans(iso_table3, centers = k, nstart = 25)
str(kmeanscluster)
fviz_cluster(kmeanscluster, data = iso_table3)
dev.copy(png, width = 1000,height = 1000, paste0(timepoint, "-K-means_cluster.png"))
dev.off()

#k3 <- kmeans(iso_table3, centers = 3, nstart = 25)
#k4 <- kmeans(iso_table3, centers = 4, nstart = 25)
#k5 <- kmeans(iso_table3, centers = 5, nstart = 25)

# plots to compare
#p1 <- fviz_cluster(k2, geom = "point", data = iso_table3) + ggtitle("k = 2")
#dev.off()
#p2 <- fviz_cluster(k3, geom = "point",  data = iso_table3) + ggtitle("k = 3")
#dev.off()
#p3 <- fviz_cluster(k4, geom = "point",  data = iso_table3) + ggtitle("k = 4")
#dev.off()
#p4 <- fviz_cluster(k5, geom = "point",  data = iso_table3) + ggtitle("k = 5")
#dev.off()

#library(gridExtra)
#grid.arrange(p1, p2, p3, p4, nrow = 2)



