# Set working directory
setwd("/Volumes/Nam_HDD/Ref data (sc_bulk RNAseq)/CIA")

# CPU core
options(mc.cores = parallel::detectCores())
options(future.globals.maxSize = 1000 * 1024^2)

# Library
library(RColorBrewer)
library(pheatmap)
library(tidyverse)
library(data.table)
library(dplyr)
library(plyr)
library(ggplot2)
library(multtest)
library(Seurat)
library(R.utils)
library(rhdf5)
library(Ecfun)
library(sctransform)
library(scater)
library(scran)
library(SingleCellExperiment)
library(gplots)
library(ggbeeswarm)
library(ggthemes)
library(slingshot)
library(gam)
library(viridis)
library(autoimage)
library(clusterExperiment)
library(readxl)

set.seed(1234567)

# Load data
merge.integrated_4 <- readRDS("./Save_new/seuset_merge_Ctr_CIA_JAK1i_30_3.RDS")
## DE of each cluster
DimPlot(merge.integrated_4, label = TRUE, label.size = 2) + NoLegend()

merge.integrated_4$clusters.group <- paste(Idents(merge.integrated_4), merge.integrated_4$group, sep = "_")
merge.integrated_4$clusters <- Idents(merge.integrated_4)
clusters <- as.data.frame(table(merge.integrated_4$clusters))

# OSM (mouse)
list1 <- read_excel("/Volumes/Nam_HDD/Ref data (sc_bulk RNAseq)/phs001529.v1.p1/Results/JAK1-cytokine/Ref/Osm/jbc.M116.748483-1.xlsx",sheet=1, col_names = F)
colnames(list1) <- list1[2,]
list1 <- list1[-c(1:2),]
list2 <- read_excel("/Volumes/Nam_HDD/Ref data (sc_bulk RNAseq)/phs001529.v1.p1/Results/JAK1-cytokine/Ref/Osm/jbc.M116.748483-1.xlsx",sheet=2, col_names = F)
colnames(list2) <- list2[2,]
list2 <- list2[-c(1:2),]
list2 <- subset(list2, list2$logFC>0)
list3 <- read_excel("/Volumes/Nam_HDD/Ref data (sc_bulk RNAseq)/phs001529.v1.p1/Results/JAK1-cytokine/Ref/Osm/Table1.xlsx",sheet=1, col_names = T)
list3 <- list3[-1,]
list <- data.frame(c(list1$Symbol,list2$Symbol,list3$Rank))
colnames(list) <- "Gene.symbol"
list[list=="NA"]<-NA
listn_osm <- na.omit(list)
listn_osm <- data.frame(listn_osm[!duplicated(listn_osm$Gene.symbol),]) # remove duplicated genes
colnames(listn_osm) <- "Gene.symbol"

#IL6 (mouse)
list <- read.table("/Volumes/Nam_HDD/Ref data (sc_bulk RNAseq)/phs001529.v1.p1/Results/JAK1-cytokine/Ref/IL6/result.txt",sep="\t", header = T)
list <- subset(list, list$adj.P.Val<0.05)
list <- subset(list, list$logFC>0)
list[list==""]<-NA
listn_il6 <- na.omit(list)
listn_il6 <- listn_il6[!duplicated(listn_il6$Gene.symbol),] # remove duplicated genes

# For Fibro_2
Idents(merge.integrated_4) <- "clusters.group"
marker <- FindMarkers(merge.integrated_4, ident.1 = "Fibro_2_B_CIA", ident.2 = "Fibro_2_C_JAK1i", 
                      only.pos = F, min.pct = 0.1, logfc.threshold = 0.25, verbose = FALSE)
marker <- subset(marker, marker$p_val_adj < 0.05)
marker$gene <- rownames(marker)
marker <- dplyr ::arrange(marker, desc(marker$avg_log2FC))
Idents(merge.integrated_4) <- "clusters"
merge.integrated_4_sub <- subset(merge.integrated_4, idents = "Fibro_2")
Idents(merge.integrated_4_sub) <- "group"
avg <- data.frame(log1p(AverageExpression(merge.integrated_4_sub, verbose = FALSE)$RNA))
avg$gene <- rownames(avg)

#
cluster.averages <- AverageExpression(merge.integrated_4_sub, return.seurat = T)
n <- which(marker$gene=="Tnfsf11", arr.ind=TRUE)
gene_for_heatmap <- marker[c(1:25,n,(nrow(marker)-24):nrow(marker)),]
pdf("./Output_new4_30/Try13_merge_Ctr_CIA_JAK1i/avg.Fibro2_CIAvsJAK1i_heatmap.pdf", useDingbats = F, height = 10, width = 4)
DoHeatmap(subset(cluster.averages, idents = c("B_CIA","C_JAK1i")), features = gene_for_heatmap$gene, size = 3, draw.lines = F, raster =F)
dev.off()

#
list.merge_osm <- merge(x=gene_for_heatmap, y=listn_osm, by.x="gene", by.y="Gene.symbol")
list.merge_osm <- dplyr :: arrange(list.merge_osm, desc(list.merge_osm$avg_log2FC))
write.table(list.merge_osm,"./Output_new4_30/Try13_merge_Ctr_CIA_JAK1i/avg.Fibro2_CIAvsJAK1i_heatmap_OSM_target.txt", sep="\t", quote=FALSE, row.names=F)

list.merge_il6 <- merge(x=gene_for_heatmap, y=listn_il6, by.x="gene", by.y="Gene.symbol")
list.merge_il6 <- dplyr :: arrange(list.merge_il6, desc(list.merge_il6$avg_log2FC))
write.table(list.merge_il6,"./Output_new4_30/Try13_merge_Ctr_CIA_JAK1i/avg.Fibro2_CIAvsJAK1i_heatmap_Il6_target.txt", sep="\t", quote=FALSE, row.names=F)

# For Fibro_3
Idents(merge.integrated_4) <- "clusters.group"
marker <- FindMarkers(merge.integrated_4, ident.1 = "Fibro_3_B_CIA", ident.2 = "Fibro_3_C_JAK1i", 
                      only.pos = F, min.pct = 0.1, logfc.threshold = 0.25, verbose = FALSE)
marker <- subset(marker, marker$p_val_adj < 0.05)
marker$gene <- rownames(marker)
marker <- dplyr ::arrange(marker, desc(marker$avg_log2FC))
Idents(merge.integrated_4) <- "clusters"
merge.integrated_4_sub <- subset(merge.integrated_4, idents = "Fibro_3")
Idents(merge.integrated_4_sub) <- "group"
avg <- log1p(AverageExpression(merge.integrated_4_sub, verbose = FALSE)$RNA)
avg$gene <- rownames(avg)

#
cluster.averages <- AverageExpression(merge.integrated_4_sub, return.seurat = T)
n <- which(marker$gene=="Tnfsf11", arr.ind=TRUE)
gene_for_heatmap <- marker[c(1:25,n,(nrow(marker)-24):nrow(marker)),]
pdf("./Output_new4_30/Try13_merge_Ctr_CIA_JAK1i/avg.Fibro3_CIAvsJAK1i_heatmap.pdf", useDingbats = F, height = 10, width = 4)
DoHeatmap(subset(cluster.averages, idents = c("B_CIA","C_JAK1i")), features = gene_for_heatmap$gene, size = 3, draw.lines = F, raster =F)
dev.off()

#
list.merge_osm <- merge(x=gene_for_heatmap, y=listn_osm, by.x="gene", by.y="Gene.symbol")
list.merge_osm <- dplyr :: arrange(list.merge_osm, desc(list.merge_osm$avg_log2FC))
write.table(list.merge_osm,"./Output_new4_30/Try13_merge_Ctr_CIA_JAK1i/avg.Fibro3_CIAvsJAK1i_heatmap_OSM_target.txt", sep="\t", quote=FALSE, row.names=F)

list.merge_il6 <- merge(x=gene_for_heatmap, y=listn_il6, by.x="gene", by.y="Gene.symbol")
list.merge_il6 <- dplyr :: arrange(list.merge_il6, desc(list.merge_il6$avg_log2FC))
write.table(list.merge_il6,"./Output_new4_30/Try13_merge_Ctr_CIA_JAK1i/avg.Fibro3_CIAvsJAK1i_heatmap_Il6_target.txt", sep="\t", quote=FALSE, row.names=F)

# For myel_c2
Idents(merge.integrated_4) <- "clusters.group"
marker <- FindMarkers(merge.integrated_4, ident.1 = "Myel_c2_B_CIA", ident.2 = "Myel_c2_C_JAK1i", 
                      only.pos = F, min.pct = 0.1, logfc.threshold = 0.25, verbose = FALSE)
marker <- subset(marker, marker$p_val_adj < 0.05)
marker$gene <- rownames(marker)
marker <- dplyr ::arrange(marker, desc(marker$avg_log2FC))
Idents(merge.integrated_4) <- "clusters"
merge.integrated_4_sub <- subset(merge.integrated_4, idents = "Myel_c2")
Idents(merge.integrated_4_sub) <- "group"
#
cluster.averages <- AverageExpression(merge.integrated_4_sub, return.seurat = T)
n <- which(marker$gene=="Tnfsf11", arr.ind=TRUE)
gene_for_heatmap <- marker[c(1:25,n,(nrow(marker)-24):nrow(marker)),]
pdf("./Output_new4_30/Try13_merge_Ctr_CIA_JAK1i/avg.Myel_c2_CIAvsJAK1i_heatmap.pdf", useDingbats = F, height = 10, width = 4)
DoHeatmap(subset(cluster.averages, idents = c("B_CIA","C_JAK1i")), features = gene_for_heatmap$gene, size = 3, draw.lines = F, raster =F)
dev.off()
#
list.merge_osm <- merge(x=gene_for_heatmap, y=listn_osm, by.x="gene", by.y="Gene.symbol")
list.merge_osm <- dplyr :: arrange(list.merge_osm, desc(list.merge_osm$avg_log2FC))
write.table(list.merge_osm,"./Output_new4_30/Try13_merge_Ctr_CIA_JAK1i/avg.Myel_c2_CIAvsJAK1i_heatmap_OSM_target.txt", sep="\t", quote=FALSE, row.names=F)
#
list.merge_il6 <- merge(x=gene_for_heatmap, y=listn_il6, by.x="gene", by.y="Gene.symbol")
list.merge_il6 <- dplyr :: arrange(list.merge_il6, desc(list.merge_il6$avg_log2FC))
write.table(list.merge_il6,"./Output_new4_30/Try13_merge_Ctr_CIA_JAK1i/avg.Myel_c2_CIAvsJAK1i_heatmap_Il6_target.txt", sep="\t", quote=FALSE, row.names=F)

# For myel_c4
Idents(merge.integrated_4) <- "clusters.group"
marker <- FindMarkers(merge.integrated_4, ident.1 = "Myel_c4_B_CIA", ident.2 = "Myel_c4_C_JAK1i", 
                      only.pos = F, min.pct = 0.1, logfc.threshold = 0.25, verbose = FALSE)
marker <- subset(marker, marker$p_val_adj < 0.05)
marker$gene <- rownames(marker)
marker <- dplyr ::arrange(marker, desc(marker$avg_log2FC))
Idents(merge.integrated_4) <- "clusters"
merge.integrated_4_sub <- subset(merge.integrated_4, idents = "Myel_c4")
Idents(merge.integrated_4_sub) <- "group"
#
cluster.averages <- AverageExpression(merge.integrated_4_sub, return.seurat = T)
n <- which(marker$gene=="Tnfsf11", arr.ind=TRUE)
gene_for_heatmap <- marker[c(1:25,n,(nrow(marker)-24):nrow(marker)),]
pdf("./Output_new4_30/Try13_merge_Ctr_CIA_JAK1i/avg.Myel_c4_CIAvsJAK1i_heatmap.pdf", useDingbats = F, height = 10, width = 4)
DoHeatmap(subset(cluster.averages, idents = c("B_CIA","C_JAK1i")), features = gene_for_heatmap$gene, size = 3, draw.lines = F, raster =F)
dev.off()
#
list.merge_osm <- merge(x=gene_for_heatmap, y=listn_osm, by.x="gene", by.y="Gene.symbol")
list.merge_osm <- dplyr :: arrange(list.merge_osm, desc(list.merge_osm$avg_log2FC))
write.table(list.merge_osm,"./Output_new4_30/Try13_merge_Ctr_CIA_JAK1i/avg.Myel_c4_CIAvsJAK1i_heatmap_OSM_target.txt", sep="\t", quote=FALSE, row.names=F)
#
list.merge_il6 <- merge(x=gene_for_heatmap, y=listn_il6, by.x="gene", by.y="Gene.symbol")
list.merge_il6 <- dplyr :: arrange(list.merge_il6, desc(list.merge_il6$avg_log2FC))
write.table(list.merge_il6,"./Output_new4_30/Try13_merge_Ctr_CIA_JAK1i/avg.Myel_c4_CIAvsJAK1i_heatmap_Il6_target.txt", sep="\t", quote=FALSE, row.names=F)

# Number of Osm, Il6 target genes in sig genes
# Fibro_1
marker_CIAvsJAK1i <- read.table("./Output_new4_30/Try13_merge_Ctr_CIA_JAK1i/marker_Fibro_1_CIAvsJAK1i.txt", sep="\t", header = T)
list.merge_osm <- merge(x=marker_CIAvsJAK1i, y=listn_osm, by.x="gene", by.y="Gene.symbol")
list.merge_osm <- dplyr :: arrange(list.merge_osm, desc(list.merge_osm$avg_log2FC))
list.merge_il6 <- merge(x=marker_CIAvsJAK1i, y=listn_il6, by.x="gene", by.y="Gene.symbol")
list.merge_il6 <- dplyr :: arrange(list.merge_il6, desc(list.merge_il6$avg_log2FC))

marker_CIA <- subset(marker_CIAvsJAK1i, marker_CIAvsJAK1i$avg_log2FC>0)
list.merge_osm <- merge(x=marker_CIA, y=listn_osm, by.x="gene", by.y="Gene.symbol")
list.merge_osm <- dplyr :: arrange(list.merge_osm, desc(list.merge_osm$avg_log2FC))
list.merge_il6 <- merge(x=marker_CIA, y=listn_il6, by.x="gene", by.y="Gene.symbol")
list.merge_il6 <- dplyr :: arrange(list.merge_il6, desc(list.merge_il6$avg_log2FC))

marker_JAK1i<- subset(marker_CIAvsJAK1i, marker_CIAvsJAK1i$avg_log2FC<0)
list.merge_osm <- merge(x=marker_JAK1i, y=listn_osm, by.x="gene", by.y="Gene.symbol")
list.merge_osm <- dplyr :: arrange(list.merge_osm, desc(list.merge_osm$avg_log2FC))
list.merge_il6 <- merge(x=marker_JAK1i, y=listn_il6, by.x="gene", by.y="Gene.symbol")
list.merge_il6 <- dplyr :: arrange(list.merge_il6, desc(list.merge_il6$avg_log2FC))

# Fibro_2
marker_CIAvsJAK1i <- read.table("./Output_new4_30/Try13_merge_Ctr_CIA_JAK1i/marker_Fibro_2_CIAvsJAK1i.txt", sep="\t", header = T)
list.merge_osm <- merge(x=marker_CIAvsJAK1i, y=listn_osm, by.x="gene", by.y="Gene.symbol")
list.merge_osm <- dplyr :: arrange(list.merge_osm, desc(list.merge_osm$avg_log2FC))
write.table(list.merge_osm,"./Output_new4_30/Try13_merge_Ctr_CIA_JAK1i/avg.Fibro2_CIAvsJAK1i_total_OSM_target.txt", sep="\t", quote=FALSE, row.names=F)

list.merge_il6 <- merge(x=marker_CIAvsJAK1i, y=listn_il6, by.x="gene", by.y="Gene.symbol")
list.merge_il6 <- dplyr :: arrange(list.merge_il6, desc(list.merge_il6$avg_log2FC))

marker_CIA <- subset(marker_CIAvsJAK1i, marker_CIAvsJAK1i$avg_log2FC>0)
list.merge_osm <- merge(x=marker_CIA, y=listn_osm, by.x="gene", by.y="Gene.symbol")
list.merge_osm <- dplyr :: arrange(list.merge_osm, desc(list.merge_osm$avg_log2FC))
list.merge_il6 <- merge(x=marker_CIA, y=listn_il6, by.x="gene", by.y="Gene.symbol")
list.merge_il6 <- dplyr :: arrange(list.merge_il6, desc(list.merge_il6$avg_log2FC))

marker_JAK1i<- subset(marker_CIAvsJAK1i, marker_CIAvsJAK1i$avg_log2FC<0)
list.merge_osm <- merge(x=marker_JAK1i, y=listn_osm, by.x="gene", by.y="Gene.symbol")
list.merge_osm <- dplyr :: arrange(list.merge_osm, desc(list.merge_osm$avg_log2FC))
list.merge_il6 <- merge(x=marker_JAK1i, y=listn_il6, by.x="gene", by.y="Gene.symbol")
list.merge_il6 <- dplyr :: arrange(list.merge_il6, desc(list.merge_il6$avg_log2FC))

# Fibro_3
marker_CIAvsJAK1i <- read.table("./Output_new4_30/Try13_merge_Ctr_CIA_JAK1i/marker_Fibro_3_CIAvsJAK1i.txt", sep="\t", header = T)
list.merge_osm <- merge(x=marker_CIAvsJAK1i, y=listn_osm, by.x="gene", by.y="Gene.symbol")
list.merge_osm <- dplyr :: arrange(list.merge_osm, desc(list.merge_osm$avg_log2FC))
write.table(list.merge_osm,"./Output_new4_30/Try13_merge_Ctr_CIA_JAK1i/avg.Fibro3_CIAvsJAK1i_total_OSM_target.txt", sep="\t", quote=FALSE, row.names=F)

list.merge_il6 <- merge(x=marker_CIAvsJAK1i, y=listn_il6, by.x="gene", by.y="Gene.symbol")
list.merge_il6 <- dplyr :: arrange(list.merge_il6, desc(list.merge_il6$avg_log2FC))

marker_CIA <- subset(marker_CIAvsJAK1i, marker_CIAvsJAK1i$avg_log2FC>0)
list.merge_osm <- merge(x=marker_CIA, y=listn_osm, by.x="gene", by.y="Gene.symbol")
list.merge_osm <- dplyr :: arrange(list.merge_osm, desc(list.merge_osm$avg_log2FC))
list.merge_il6 <- merge(x=marker_CIA, y=listn_il6, by.x="gene", by.y="Gene.symbol")
list.merge_il6 <- dplyr :: arrange(list.merge_il6, desc(list.merge_il6$avg_log2FC))

marker_JAK1i<- subset(marker_CIAvsJAK1i, marker_CIAvsJAK1i$avg_log2FC<0)
list.merge_osm <- merge(x=marker_JAK1i, y=listn_osm, by.x="gene", by.y="Gene.symbol")
list.merge_osm <- dplyr :: arrange(list.merge_osm, desc(list.merge_osm$avg_log2FC))
list.merge_il6 <- merge(x=marker_JAK1i, y=listn_il6, by.x="gene", by.y="Gene.symbol")
list.merge_il6 <- dplyr :: arrange(list.merge_il6, desc(list.merge_il6$avg_log2FC))

# Fibro_4
marker_CIAvsJAK1i <- read.table("./Output_new4_30/Try13_merge_Ctr_CIA_JAK1i/marker_Fibro_4_CIAvsJAK1i.txt", sep="\t", header = T)
list.merge_osm <- merge(x=marker_CIAvsJAK1i, y=listn_osm, by.x="gene", by.y="Gene.symbol")
list.merge_osm <- dplyr :: arrange(list.merge_osm, desc(list.merge_osm$avg_log2FC))

list.merge_il6 <- merge(x=marker_CIAvsJAK1i, y=listn_il6, by.x="gene", by.y="Gene.symbol")
list.merge_il6 <- dplyr :: arrange(list.merge_il6, desc(list.merge_il6$avg_log2FC))

marker_CIA <- subset(marker_CIAvsJAK1i, marker_CIAvsJAK1i$avg_log2FC>0)
list.merge_osm <- merge(x=marker_CIA, y=listn_osm, by.x="gene", by.y="Gene.symbol")
list.merge_osm <- dplyr :: arrange(list.merge_osm, desc(list.merge_osm$avg_log2FC))
list.merge_il6 <- merge(x=marker_CIA, y=listn_il6, by.x="gene", by.y="Gene.symbol")
list.merge_il6 <- dplyr :: arrange(list.merge_il6, desc(list.merge_il6$avg_log2FC))

marker_JAK1i<- subset(marker_CIAvsJAK1i, marker_CIAvsJAK1i$avg_log2FC<0)
list.merge_osm <- merge(x=marker_JAK1i, y=listn_osm, by.x="gene", by.y="Gene.symbol")
list.merge_osm <- dplyr :: arrange(list.merge_osm, desc(list.merge_osm$avg_log2FC))
list.merge_il6 <- merge(x=marker_JAK1i, y=listn_il6, by.x="gene", by.y="Gene.symbol")
list.merge_il6 <- dplyr :: arrange(list.merge_il6, desc(list.merge_il6$avg_log2FC))

# Myel_a1
marker_CIAvsJAK1i <- read.table("./Output_new4_30/Try13_merge_Ctr_CIA_JAK1i/marker_Myel_a1_CIAvsJAK1i.txt", sep="\t", header = T)
list.merge_osm <- merge(x=marker_CIAvsJAK1i, y=listn_osm, by.x="gene", by.y="Gene.symbol")
list.merge_osm <- dplyr :: arrange(list.merge_osm, desc(list.merge_osm$avg_log2FC))

list.merge_il6 <- merge(x=marker_CIAvsJAK1i, y=listn_il6, by.x="gene", by.y="Gene.symbol")
list.merge_il6 <- dplyr :: arrange(list.merge_il6, desc(list.merge_il6$avg_log2FC))

marker_CIA <- subset(marker_CIAvsJAK1i, marker_CIAvsJAK1i$avg_log2FC>0)
list.merge_osm <- merge(x=marker_CIA, y=listn_osm, by.x="gene", by.y="Gene.symbol")
list.merge_osm <- dplyr :: arrange(list.merge_osm, desc(list.merge_osm$avg_log2FC))
list.merge_il6 <- merge(x=marker_CIA, y=listn_il6, by.x="gene", by.y="Gene.symbol")
list.merge_il6 <- dplyr :: arrange(list.merge_il6, desc(list.merge_il6$avg_log2FC))

marker_JAK1i<- subset(marker_CIAvsJAK1i, marker_CIAvsJAK1i$avg_log2FC<0)
list.merge_osm <- merge(x=marker_JAK1i, y=listn_osm, by.x="gene", by.y="Gene.symbol")
list.merge_osm <- dplyr :: arrange(list.merge_osm, desc(list.merge_osm$avg_log2FC))
list.merge_il6 <- merge(x=marker_JAK1i, y=listn_il6, by.x="gene", by.y="Gene.symbol")
list.merge_il6 <- dplyr :: arrange(list.merge_il6, desc(list.merge_il6$avg_log2FC))

# Myel_a2
marker_CIAvsJAK1i <- read.table("./Output_new4_30/Try13_merge_Ctr_CIA_JAK1i/marker_Myel_a2_CIAvsJAK1i.txt", sep="\t", header = T)
list.merge_osm <- merge(x=marker_CIAvsJAK1i, y=listn_osm, by.x="gene", by.y="Gene.symbol")
list.merge_osm <- dplyr :: arrange(list.merge_osm, desc(list.merge_osm$avg_log2FC))

list.merge_il6 <- merge(x=marker_CIAvsJAK1i, y=listn_il6, by.x="gene", by.y="Gene.symbol")
list.merge_il6 <- dplyr :: arrange(list.merge_il6, desc(list.merge_il6$avg_log2FC))

marker_CIA <- subset(marker_CIAvsJAK1i, marker_CIAvsJAK1i$avg_log2FC>0)
list.merge_osm <- merge(x=marker_CIA, y=listn_osm, by.x="gene", by.y="Gene.symbol")
list.merge_osm <- dplyr :: arrange(list.merge_osm, desc(list.merge_osm$avg_log2FC))
list.merge_il6 <- merge(x=marker_CIA, y=listn_il6, by.x="gene", by.y="Gene.symbol")
list.merge_il6 <- dplyr :: arrange(list.merge_il6, desc(list.merge_il6$avg_log2FC))

marker_JAK1i<- subset(marker_CIAvsJAK1i, marker_CIAvsJAK1i$avg_log2FC<0)
list.merge_osm <- merge(x=marker_JAK1i, y=listn_osm, by.x="gene", by.y="Gene.symbol")
list.merge_osm <- dplyr :: arrange(list.merge_osm, desc(list.merge_osm$avg_log2FC))
list.merge_il6 <- merge(x=marker_JAK1i, y=listn_il6, by.x="gene", by.y="Gene.symbol")
list.merge_il6 <- dplyr :: arrange(list.merge_il6, desc(list.merge_il6$avg_log2FC))

# Myel_a3
marker_CIAvsJAK1i <- read.table("./Output_new4_30/Try13_merge_Ctr_CIA_JAK1i/marker_Myel_a3_CIAvsJAK1i.txt", sep="\t", header = T)
list.merge_osm <- merge(x=marker_CIAvsJAK1i, y=listn_osm, by.x="gene", by.y="Gene.symbol")
list.merge_osm <- dplyr :: arrange(list.merge_osm, desc(list.merge_osm$avg_log2FC))

list.merge_il6 <- merge(x=marker_CIAvsJAK1i, y=listn_il6, by.x="gene", by.y="Gene.symbol")
list.merge_il6 <- dplyr :: arrange(list.merge_il6, desc(list.merge_il6$avg_log2FC))

marker_CIA <- subset(marker_CIAvsJAK1i, marker_CIAvsJAK1i$avg_log2FC>0)
list.merge_osm <- merge(x=marker_CIA, y=listn_osm, by.x="gene", by.y="Gene.symbol")
list.merge_osm <- dplyr :: arrange(list.merge_osm, desc(list.merge_osm$avg_log2FC))
list.merge_il6 <- merge(x=marker_CIA, y=listn_il6, by.x="gene", by.y="Gene.symbol")
list.merge_il6 <- dplyr :: arrange(list.merge_il6, desc(list.merge_il6$avg_log2FC))

marker_JAK1i<- subset(marker_CIAvsJAK1i, marker_CIAvsJAK1i$avg_log2FC<0)
list.merge_osm <- merge(x=marker_JAK1i, y=listn_osm, by.x="gene", by.y="Gene.symbol")
list.merge_osm <- dplyr :: arrange(list.merge_osm, desc(list.merge_osm$avg_log2FC))
list.merge_il6 <- merge(x=marker_JAK1i, y=listn_il6, by.x="gene", by.y="Gene.symbol")
list.merge_il6 <- dplyr :: arrange(list.merge_il6, desc(list.merge_il6$avg_log2FC))

# Myel_b
marker_CIAvsJAK1i <- read.table("./Output_new4_30/Try13_merge_Ctr_CIA_JAK1i/marker_Myeloid_b_CIAvsJAK1i.txt", sep="\t", header = T)
list.merge_osm <- merge(x=marker_CIAvsJAK1i, y=listn_osm, by.x="gene", by.y="Gene.symbol")
list.merge_osm <- dplyr :: arrange(list.merge_osm, desc(list.merge_osm$avg_log2FC))

list.merge_il6 <- merge(x=marker_CIAvsJAK1i, y=listn_il6, by.x="gene", by.y="Gene.symbol")
list.merge_il6 <- dplyr :: arrange(list.merge_il6, desc(list.merge_il6$avg_log2FC))

marker_CIA <- subset(marker_CIAvsJAK1i, marker_CIAvsJAK1i$avg_log2FC>0)
list.merge_osm <- merge(x=marker_CIA, y=listn_osm, by.x="gene", by.y="Gene.symbol")
list.merge_osm <- dplyr :: arrange(list.merge_osm, desc(list.merge_osm$avg_log2FC))
list.merge_il6 <- merge(x=marker_CIA, y=listn_il6, by.x="gene", by.y="Gene.symbol")
list.merge_il6 <- dplyr :: arrange(list.merge_il6, desc(list.merge_il6$avg_log2FC))

marker_JAK1i<- subset(marker_CIAvsJAK1i, marker_CIAvsJAK1i$avg_log2FC<0)
list.merge_osm <- merge(x=marker_JAK1i, y=listn_osm, by.x="gene", by.y="Gene.symbol")
list.merge_osm <- dplyr :: arrange(list.merge_osm, desc(list.merge_osm$avg_log2FC))
list.merge_il6 <- merge(x=marker_JAK1i, y=listn_il6, by.x="gene", by.y="Gene.symbol")
list.merge_il6 <- dplyr :: arrange(list.merge_il6, desc(list.merge_il6$avg_log2FC))

# Myel_c1
marker_CIAvsJAK1i <- read.table("./Output_new4_30/Try13_merge_Ctr_CIA_JAK1i/marker_Myel_c1_CIAvsJAK1i.txt", sep="\t", header = T)
list.merge_osm <- merge(x=marker_CIAvsJAK1i, y=listn_osm, by.x="gene", by.y="Gene.symbol")
list.merge_osm <- dplyr :: arrange(list.merge_osm, desc(list.merge_osm$avg_log2FC))

list.merge_il6 <- merge(x=marker_CIAvsJAK1i, y=listn_il6, by.x="gene", by.y="Gene.symbol")
list.merge_il6 <- dplyr :: arrange(list.merge_il6, desc(list.merge_il6$avg_log2FC))

marker_CIA <- subset(marker_CIAvsJAK1i, marker_CIAvsJAK1i$avg_log2FC>0)
list.merge_osm <- merge(x=marker_CIA, y=listn_osm, by.x="gene", by.y="Gene.symbol")
list.merge_osm <- dplyr :: arrange(list.merge_osm, desc(list.merge_osm$avg_log2FC))
list.merge_il6 <- merge(x=marker_CIA, y=listn_il6, by.x="gene", by.y="Gene.symbol")
list.merge_il6 <- dplyr :: arrange(list.merge_il6, desc(list.merge_il6$avg_log2FC))

marker_JAK1i<- subset(marker_CIAvsJAK1i, marker_CIAvsJAK1i$avg_log2FC<0)
list.merge_osm <- merge(x=marker_JAK1i, y=listn_osm, by.x="gene", by.y="Gene.symbol")
list.merge_osm <- dplyr :: arrange(list.merge_osm, desc(list.merge_osm$avg_log2FC))
list.merge_il6 <- merge(x=marker_JAK1i, y=listn_il6, by.x="gene", by.y="Gene.symbol")
list.merge_il6 <- dplyr :: arrange(list.merge_il6, desc(list.merge_il6$avg_log2FC))

# Myel_c2
marker_CIAvsJAK1i <- read.table("./Output_new4_30/Try13_merge_Ctr_CIA_JAK1i/marker_Myel_c2_CIAvsJAK1i.txt", sep="\t", header = T)
list.merge_osm <- merge(x=marker_CIAvsJAK1i, y=listn_osm, by.x="gene", by.y="Gene.symbol")
list.merge_osm <- dplyr :: arrange(list.merge_osm, desc(list.merge_osm$avg_log2FC))
write.table(list.merge_osm,"./Output_new4_30/Try13_merge_Ctr_CIA_JAK1i/avg.Myel_c2_CIAvsJAK1i_total_OSM_target.txt", sep="\t", quote=FALSE, row.names=F)

list.merge_il6 <- merge(x=marker_CIAvsJAK1i, y=listn_il6, by.x="gene", by.y="Gene.symbol")
list.merge_il6 <- dplyr :: arrange(list.merge_il6, desc(list.merge_il6$avg_log2FC))

marker_CIA <- subset(marker_CIAvsJAK1i, marker_CIAvsJAK1i$avg_log2FC>0)
list.merge_osm <- merge(x=marker_CIA, y=listn_osm, by.x="gene", by.y="Gene.symbol")
list.merge_osm <- dplyr :: arrange(list.merge_osm, desc(list.merge_osm$avg_log2FC))
list.merge_il6 <- merge(x=marker_CIA, y=listn_il6, by.x="gene", by.y="Gene.symbol")
list.merge_il6 <- dplyr :: arrange(list.merge_il6, desc(list.merge_il6$avg_log2FC))

marker_JAK1i<- subset(marker_CIAvsJAK1i, marker_CIAvsJAK1i$avg_log2FC<0)
list.merge_osm <- merge(x=marker_JAK1i, y=listn_osm, by.x="gene", by.y="Gene.symbol")
list.merge_osm <- dplyr :: arrange(list.merge_osm, desc(list.merge_osm$avg_log2FC))
list.merge_il6 <- merge(x=marker_JAK1i, y=listn_il6, by.x="gene", by.y="Gene.symbol")
list.merge_il6 <- dplyr :: arrange(list.merge_il6, desc(list.merge_il6$avg_log2FC))

# Myel_c3
marker_CIAvsJAK1i <- read.table("./Output_new4_30/Try13_merge_Ctr_CIA_JAK1i/marker_Myel_c3_CIAvsJAK1i.txt", sep="\t", header = T)
list.merge_osm <- merge(x=marker_CIAvsJAK1i, y=listn_osm, by.x="gene", by.y="Gene.symbol")
list.merge_osm <- dplyr :: arrange(list.merge_osm, desc(list.merge_osm$avg_log2FC))

list.merge_il6 <- merge(x=marker_CIAvsJAK1i, y=listn_il6, by.x="gene", by.y="Gene.symbol")
list.merge_il6 <- dplyr :: arrange(list.merge_il6, desc(list.merge_il6$avg_log2FC))

marker_CIA <- subset(marker_CIAvsJAK1i, marker_CIAvsJAK1i$avg_log2FC>0)
list.merge_osm <- merge(x=marker_CIA, y=listn_osm, by.x="gene", by.y="Gene.symbol")
list.merge_osm <- dplyr :: arrange(list.merge_osm, desc(list.merge_osm$avg_log2FC))
list.merge_il6 <- merge(x=marker_CIA, y=listn_il6, by.x="gene", by.y="Gene.symbol")
list.merge_il6 <- dplyr :: arrange(list.merge_il6, desc(list.merge_il6$avg_log2FC))

marker_JAK1i<- subset(marker_CIAvsJAK1i, marker_CIAvsJAK1i$avg_log2FC<0)
list.merge_osm <- merge(x=marker_JAK1i, y=listn_osm, by.x="gene", by.y="Gene.symbol")
list.merge_osm <- dplyr :: arrange(list.merge_osm, desc(list.merge_osm$avg_log2FC))
list.merge_il6 <- merge(x=marker_JAK1i, y=listn_il6, by.x="gene", by.y="Gene.symbol")
list.merge_il6 <- dplyr :: arrange(list.merge_il6, desc(list.merge_il6$avg_log2FC))

# Myel_c4
marker_CIAvsJAK1i <- read.table("./Output_new4_30/Try13_merge_Ctr_CIA_JAK1i/marker_Myel_c4_CIAvsJAK1i.txt", sep="\t", header = T)
list.merge_osm <- merge(x=marker_CIAvsJAK1i, y=listn_osm, by.x="gene", by.y="Gene.symbol")
list.merge_osm <- dplyr :: arrange(list.merge_osm, desc(list.merge_osm$avg_log2FC))
write.table(list.merge_osm,"./Output_new4_30/Try13_merge_Ctr_CIA_JAK1i/avg.Myel_c4_CIAvsJAK1i_total_OSM_target.txt", sep="\t", quote=FALSE, row.names=F)

list.merge_il6 <- merge(x=marker_CIAvsJAK1i, y=listn_il6, by.x="gene", by.y="Gene.symbol")
list.merge_il6 <- dplyr :: arrange(list.merge_il6, desc(list.merge_il6$avg_log2FC))

marker_CIA <- subset(marker_CIAvsJAK1i, marker_CIAvsJAK1i$avg_log2FC>0)
list.merge_osm <- merge(x=marker_CIA, y=listn_osm, by.x="gene", by.y="Gene.symbol")
list.merge_osm <- dplyr :: arrange(list.merge_osm, desc(list.merge_osm$avg_log2FC))
list.merge_il6 <- merge(x=marker_CIA, y=listn_il6, by.x="gene", by.y="Gene.symbol")
list.merge_il6 <- dplyr :: arrange(list.merge_il6, desc(list.merge_il6$avg_log2FC))

marker_JAK1i<- subset(marker_CIAvsJAK1i, marker_CIAvsJAK1i$avg_log2FC<0)
list.merge_osm <- merge(x=marker_JAK1i, y=listn_osm, by.x="gene", by.y="Gene.symbol")
list.merge_osm <- dplyr :: arrange(list.merge_osm, desc(list.merge_osm$avg_log2FC))
list.merge_il6 <- merge(x=marker_JAK1i, y=listn_il6, by.x="gene", by.y="Gene.symbol")
list.merge_il6 <- dplyr :: arrange(list.merge_il6, desc(list.merge_il6$avg_log2FC))


# Affected genes
Affected_gene <- read.table("./Output_new4_30/Try13_merge_Ctr_CIA_JAK1i/Affected_gene.txt", sep="\t", header = T)
row.names(Affected_gene) <- Affected_gene$Gene_number

#Affected_gene_total
Affected_gene_total <- Affected_gene[7:9,-1]

Affected_gene_total_Percentage <- rbind(Affected_gene_total, c(Affected_gene_total[1,]/Affected_gene_total[1,]),
                                        c(Affected_gene_total[2,]/Affected_gene_total[1,]),
                                        c(Affected_gene_total[3,]/Affected_gene_total[1,]))

Affected_gene_total <- Affected_gene_total%>%gather("Cell","Gene_number" ,1:12)
Affected_gene_total$Group <- c(Group,Group,Group,Group,Group,Group,Group,Group,Group,Group,Group,Group)
Group <-  c("Affected genes", "IL6 target genes", "OSM target genes")

Affected_gene_total_Percentage <- Affected_gene_total_Percentage[4:6,]%>%gather("Cell","Percentage" ,1:12)
Affected_gene_total$Percentage <- Affected_gene_total_Percentage$Percentage

# Barplot
p <- ggplot(data=Affected_gene_total, aes(x=Cell, y=Gene_number, fill=Group)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+
  geom_text(aes(label=Gene_number), vjust=-2, color="black",
            position = position_dodge(0.9), size=2)+
  geom_text(aes(label=scales::percent(Percentage,1)), vjust=-0.5, color="black",
            position = position_dodge(0.9), size=2)+
  scale_fill_manual(values=c("#424242",'#adadad','#FFFFFF')) + theme_classic()


pdf("./Output_new4_30/Try13_merge_Ctr_CIA_JAK1i/Affected_gene.pdf", useDingbats = F, height = 2.5, width = 10)
p
dev.off()

#Affected_gene_CIA
Affected_gene_CIA <- Affected_gene[1:3,-1]

Affected_gene_CIA_Percentage <- rbind(Affected_gene_CIA, c(Affected_gene_CIA[1,]/Affected_gene_CIA[1,]),
                                        c(Affected_gene_CIA[2,]/Affected_gene_CIA[1,]),
                                        c(Affected_gene_CIA[3,]/Affected_gene_CIA[1,]))

Affected_gene_CIA <- Affected_gene_CIA%>%gather("Cell","Gene_number" ,1:12)
Affected_gene_CIA$Group <- c(Group,Group,Group,Group,Group,Group,Group,Group,Group,Group,Group,Group)
Group <-  c("Affected genes", "IL6 target genes", "OSM target genes")

Affected_gene_CIA_Percentage <- Affected_gene_CIA_Percentage[4:6,]%>%gather("Cell","Percentage" ,1:12)
Affected_gene_CIA$Percentage <- Affected_gene_CIA_Percentage$Percentage

# Barplot
p <- ggplot(data=Affected_gene_CIA, aes(x=Cell, y=Gene_number, fill=Group)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+
  geom_text(aes(label=Gene_number), vjust=-2, color="black",
            position = position_dodge(0.9), size=2)+
  geom_text(aes(label=scales::percent(Percentage,1)), vjust=-0.5, color="black",
            position = position_dodge(0.9), size=2)+
  ylim(0, 400)+
  scale_fill_manual(values=c("#424242",'#adadad','#FFFFFF')) + theme_classic()


pdf("./Output_new4_30/Try13_merge_Ctr_CIA_JAK1i/Affected_gene_CIA.pdf", useDingbats = F, height = 2, width = 10)
p
dev.off()

#Affected_gene_JAK1i
Affected_gene_JAK1i <- Affected_gene[4:6,-1]
Affected_gene_JAK1i <- -Affected_gene_JAK1i

Affected_gene_JAK1i_Percentage <- rbind(Affected_gene_JAK1i, c(Affected_gene_JAK1i[1,]/Affected_gene_JAK1i[1,]),
                                      c(Affected_gene_JAK1i[2,]/Affected_gene_JAK1i[1,]),
                                      c(Affected_gene_JAK1i[3,]/Affected_gene_JAK1i[1,]))

Affected_gene_JAK1i <- Affected_gene_JAK1i%>%gather("Cell","Gene_number" ,1:12)
Affected_gene_JAK1i$Group <- c(Group,Group,Group,Group,Group,Group,Group,Group,Group,Group,Group,Group)
Group <-  c("Affected genes", "IL6 target genes", "OSM target genes")

Affected_gene_JAK1i_Percentage <- Affected_gene_JAK1i_Percentage[4:6,]%>%gather("Cell","Percentage" ,1:12)
Affected_gene_JAK1i$Percentage <- Affected_gene_JAK1i_Percentage$Percentage

# Barplot
p <- ggplot(data=Affected_gene_JAK1i, aes(x=Cell, y=Gene_number, fill=Group)) +
  geom_bar(stat="identity", color="black", position=position_dodge())+
  geom_text(aes(label=-Gene_number), vjust= 1.5 , color="black",
            position = position_dodge(0.9), size=2)+
  geom_text(aes(label=scales::percent(Percentage,1)), vjust= 3, color="black",
            position = position_dodge(0.9), size=2)+
  ylim(-400, 0)+
  scale_fill_manual(values=c("#424242",'#adadad','#FFFFFF')) + theme_classic()


pdf("./Output_new4_30/Try13_merge_Ctr_CIA_JAK1i/Affected_gene_JAK1i.pdf", useDingbats = F, height = 2, width = 10.07)
p
dev.off()

# Cell number
merge.integrated_3 <- readRDS("./Save_new/seuset_merge_Ctr_CIA_JAK1i_30_2.RDS")

DimPlot(merge.integrated_3, label = TRUE, label.size = 2) + NoLegend()
merge.integrated_3$clusters.group <- paste(Idents(merge.integrated_3), merge.integrated_3$group, sep = "_")
merge.integrated_3$clusters <- Idents(merge.integrated_3)
clusters <- as.data.frame(table(merge.integrated_3$clusters.group))

table(merge.integrated_3$group)
table(merge.integrated_3$clusters.group)

# 
library(ggmosaic)
library(plyr)

fig1 <- data.frame(Cell_type=factor(c("T_cell","B_cell","Myeloid_a","Myeloid_b","Myeloid_c",
                                        "Fibroblast","Neutrophil","Endothelial","Mural_cell")),
                   A_Ctr=c(6/181,3/181,0/181,15/181,50/181,60/181,2/181,42/181,3/181),
                   B_CIA=c(35/1278,22/1278,457/1278,43/1278,389/1278,225/1278,74/1278,22/1278,11/1278),
                   C_JAK1i=c(42/2571,45/2571,271/2571,48/2571,924/2571,890/2571,208/2571,123/2571,20/2571))

fig1$Cell_type <- ordered(fig1$Cell_type, levels = c("T_cell","B_cell","Myeloid_a","Myeloid_b","Myeloid_c",
                                                     "Fibroblast","Neutrophil","Endothelial","Mural_cell"))

fig1 <- fig1%>%gather("Group","Percentage",2:4)
fig1$Percentage <- fig1$Percentage*100
fig1$Group <- ordered(fig1$Group, levels = c("A_Ctr","B_CIA","C_JAK1i"))
fig1 <- ddply(fig1, .(Group), transform, pos = 100 - (cumsum(Percentage) - (0.5 * Percentage)))

fig1.1 <- ggplot(fig1) +
  geom_bar(aes(y = Percentage, 
                  x = Group,
                  fill = Cell_type), data = fig1, stat="identity")+
  geom_text(data=fig1, aes(x = Group, y = pos, label = paste0(round(Percentage,digits=1),"%")),
            size=3)+
  scale_fill_manual(values=c('#F8766D','#CD9600','#7CAE00','#0CB702','#00C19A',
                             '#00B8E7','#8494FF','#C77CFF','#FF61CC')) +
  guides(fill=guide_legend(title = "Cell_type")) +
  theme_classic()

pdf("./Output_new4_30/Try13_merge_Ctr_CIA_JAK1i/Cell_percentage.pdf", useDingbats = F, height = 4, width = 4)
fig1.1
dev.off()

#
fig1 <- data.frame(Cell_type=factor(c("T_cell","B_cell","Myel_a1","Myel_a2","Myel_a3","Myel_b",
                                      "Myel_c1","Myel_c2","Myel_c3","Myel_c4","Myel_c5",
                                      "Fibro_1","Fibro_2","Fibro_3","Fibro_4","Chondrocyte","Neutrophil",
                                      "Endothelial","Mural_cell")),
                   A_Ctr=c(6/181,3/181,0/181,0/181,0/181,15/181,6/181,18/181,5/181,17/181,4/181,
                           4/181,29/181,17/181,3/181,7/181,2/181,42/181,3/181),
                   B_CIA=c(35/1278,22/1278,201/1278,132/1278,124/1278,49/1278,81/1278,143/1278,20/1278,127/1278,12/1278,
                           9/1278,111/1278,68/1278,13/1278,24/1278,74/1278,22/1278,11/1278),
                   C_JAK1i=c(42/2571,45/2571,104/2571,80/2571,87/2571,58/2571,115/2571,287/2571,151/2571,345/2571,16/2571,
                             47/2571,307/2571,345/2571,136/2571,55/2571,208/2571,123/2571,20/2571))

fig1$Cell_type <- ordered(fig1$Cell_type, levels = c("T_cell","B_cell","Myel_a1","Myel_a2","Myel_a3","Myel_b",
                                                     "Myel_c1","Myel_c2","Myel_c3","Myel_c4","Myel_c5",
                                                     "Fibro_1","Fibro_2","Fibro_3","Fibro_4","Chondrocyte","Neutrophil",
                                                     "Endothelial","Mural_cell"))

fig1 <- fig1%>%gather("Group","Percentage",2:4)
fig1$Percentage <- fig1$Percentage*100
fig1$Group <- ordered(fig1$Group, levels = c("A_Ctr","B_CIA","C_JAK1i"))

fig1 <- ddply(fig1, .(Group), transform, pos = 100 - (cumsum(Percentage) - (0.5 * Percentage)))

getPalette <- colorRampPalette(brewer.pal(9, "Set1"))

fig1.1 <- ggplot(fig1) +
  geom_bar(aes(y = Percentage, 
               x = Group,
               fill = Cell_type), data = fig1, stat="identity")+
  geom_text(data=fig1, aes(x = Group, y = pos, label = paste0(round(Percentage,digits=1),"%")),
            size=3)+
  scale_fill_manual(values=getPalette(19)) +
  guides(fill=guide_legend(title = "Cell_type")) +
  theme_classic()

pdf("./Output_new4_30/Try13_merge_Ctr_CIA_JAK1i/Cell_percentage_2.pdf", useDingbats = F, height = 10, width = 3.5)
fig1.1
dev.off()

fig1.1 <- ggplot(fig1) +
  geom_bar(aes(y = Percentage, 
               x = Group,
               fill = Cell_type), data = fig1, stat="identity")+
  geom_text(data=fig1, aes(x = Group, y = pos, label = paste0(round(Percentage,digits=1),"%")),
            size=3)+
  scale_fill_manual(values=getPalette(19)) +
  guides(fill=guide_legend(title = "Cell_type")) +
  theme_classic() + coord_flip()

pdf("./Output_new4_30/Try13_merge_Ctr_CIA_JAK1i/Cell_percentage_3.pdf", useDingbats = F, height = 3.5, width = 10)
fig1.1
dev.off()

