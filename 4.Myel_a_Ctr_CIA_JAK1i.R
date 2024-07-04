# Set working directory
setwd("/Volumes/Nam_HDD/Ref data (sc_bulk RNAseq)/CIA")

# CPU core
options(mc.cores = parallel::detectCores())
options(future.globals.maxSize = 2000 * 1024^2)

# Library
library(RColorBrewer)
library(pheatmap)
library(tidyverse)
library(data.table)
library(dplyr)
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
library (clusterExperiment)

set.seed(1234567)

# Load data
merge.integrated_3 <- readRDS("./Save_new/seuset_merge_Ctr_CIA_JAK1i_30_2.RDS")

# Select Myel_a
Myel_a <- subset(merge.integrated_3, idents = c("Myeloid_a"))
Myel_a

# DefaultAssay
DefaultAssay(object = Myel_a) <- "RNA"
merge.list <- SplitObject(object = Myel_a, split.by = "group")

# SCT
for (i in 1:length(x = merge.list)) {
  merge.list[[i]] <- SCTransform(object = merge.list[[i]], verbose = FALSE
                                 , vars.to.regress = "percent.mt"
  )
}
merge.features <- SelectIntegrationFeatures(object.list = merge.list, nfeatures = 3000)
merge.list <- PrepSCTIntegration(object.list = merge.list, anchor.features = merge.features, verbose = FALSE)
merge.anchors <- FindIntegrationAnchors(object.list = merge.list, normalization.method = "SCT", 
                                        anchor.features = merge.features, verbose = FALSE)
merge.integrated <- IntegrateData(anchorset = merge.anchors, normalization.method = "SCT", verbose = FALSE)

# PCA, UMAP
merge.integrated <- RunPCA(object = merge.integrated,verbose = FALSE)
merge.integrated <- RunUMAP(object = merge.integrated, dims = 1:30)
merge.integrated_1 <- FindNeighbors(object = merge.integrated, dims = 1:30)

# Save
saveRDS(merge.integrated_1, "./Save_new/Myel_a_Ctr_CIA_JAK1i_30_1.RDS")
merge.integrated_1 <- readRDS("./Save_new/Myel_a_Ctr_CIA_JAK1i_30_1.RDS")

# Clustering
Myel_a <- FindClusters(merge.integrated_1, verbose = FALSE, resolution = 0.3)
table(Myel_a@active.ident)

Myel_a <- RenameIdents(object=Myel_a,"2"="Myel_a3")
Myel_a <- RenameIdents(object=Myel_a,"1"="Myel_a2")
Myel_a <- RenameIdents(object=Myel_a,"0"="Myel_a1")

Myel_a$saved.idents <- Idents(object = Myel_a)

# DefaultAssay
DefaultAssay(object = Myel_a) <- "RNA"
Myel_a <- NormalizeData(Myel_a, verbose = FALSE)

pdf("./Output_new4_30/Try4_Myel_a_Ctr_CIA_JAK1i/dimplot_0.5.pdf", useDingbats = F, height = 4, width = 6)
DimPlot(Myel_a, label = TRUE, label.size = 4) + NoLegend()
DimPlot(Myel_a, label = TRUE, label.size = 0, group.by = "group")
dev.off()

#
markers.to.plot <- c("Cd14","Ccl4","Il1b","Tnf",
                     "Osm","Tlr4","Cd274","Icam1",
                     "Itgam","Adgre1","Retnlg","Isg15","Irf7",
                     "Ifitm6","Rsad2","Ifit1","Fcgr4","Itga2","Ly6g")

pdf("./Output_new4_30/Try4_Myel_a_Ctr_CIA_JAK1i/featureplot_1.pdf", useDingbats = F, height = 14, width = 15)
FeaturePlot(Myel_a, markers.to.plot,
            cols = c("lightgrey", "#005493"))
dev.off()

#
pdf("./Output_new4_30/Try4_Myel_a_Ctr_CIA_JAK1i/DotPlot_0.5.pdf", useDingbats = F, height = 3, width = 8)
DotPlot(Myel_a, features = markers.to.plot, dot.scale = 8) + RotatedAxis()
dev.off()

#
pdf("./Output_new4_30/Try4_Myel_a_Ctr_CIA_JAK1i/VlnPlot_1.pdf", useDingbats = F, height = 25, width = 8)
markers.to.plot <- c("Cd14","Il1b","Tnf","Osm","Ccl3","Cd274","Retnlg","Lcn2","Lyz2","Ifitm6","Ifit1")
plots <- VlnPlot(Myel_a, features = markers.to.plot
                 , split.by = "group" , 
                 pt.size = 0.1, combine = FALSE, split.plot = F)
CombinePlots(plots = plots, ncol = 1)
dev.off()

# Markers
markers <- FindAllMarkers(object = Myel_a, only.pos = F, min.pct = 0.1, logfc.threshold = 0.1)
markers <- subset(markers, markers$p_val_adj < 0.05)
markers <- dplyr ::arrange(markers, desc(markers$avg_logFC))
markers <- dplyr ::arrange(markers, markers$cluster)
write.table(markers,"./Output_new4_30/Try4_Myel_a_Ctr_CIA_JAK1i/markers_0.5.txt", sep="\t", quote=FALSE, row.names=FALSE)

#
Myel_a$clusters.group <- paste(Idents(Myel_a), Myel_a$group, sep = "_")

#RNAseq Myel_a
avg.Myel_a <- log1p(AverageExpression(Myel_a, verbose = FALSE)$RNA)
avg.Myel_a$gene <- rownames(avg.Myel_a)
avg.markers <- FindMarkers(Myel_a, ident.1 = "B_CIA", ident.2 = "C_JAK1i", verbose = FALSE)
#avg.markers <- subset(avg.markers, avg.markers$p_val_adj < 0.05)
avg.markers$gene <- rownames(avg.markers)
avg.markers <- dplyr :: arrange(avg.markers, desc(avg.markers$avg_log2FC))
write.table(avg.markers,"./Output_new4_30/Try4_Myel_a_Ctr_CIA_JAK1i/markers_Myel_a_CIAvsJAK1i.txt", sep="\t", quote=FALSE, row.names=T)

markers.to.plot <- c("Cd14","Il1b","Tnf","Osm","Ccl3","Cxcl3","Cd274","Retnlg","Isg15","Irf7","Ifitm6","Ifit1")

pdf("./Output_new4_30/Try4_Myel_a_Ctr_CIA_JAK1i/avg.Myel_a_CIAvsJAK1i.pdf", useDingbats = F, height = 5, width = 5)
p1 <- ggplot(avg.Myel_a, aes(B_CIA, C_JAK1i)) + geom_point() +
  geom_point(data=avg.Myel_a[markers.to.plot, ], aes(B_CIA, C_JAK1i), colour="red") +
  ggtitle("Myel_a")
LabelPoints(plot = p1, points = markers.to.plot, repel = TRUE, color="red")
dev.off()

# RNAseq Myel_a1
Idents(Myel_a) <- "saved.idents"
Myel_asub <- subset(Myel_a, idents = "Myel_a1")
Idents(Myel_asub) <- "group"
avg.Myel_a <- data.frame(log1p(AverageExpression(Myel_asub, verbose = FALSE)$RNA))
avg.Myel_a$gene <- rownames(avg.Myel_a)
avg.markers <- FindMarkers(Myel_asub, ident.1 = "B_CIA", ident.2 = "C_JAK1i",
                           only.pos = F, min.pct = 0.1, logfc.threshold = 0.25, verbose = FALSE)
avg.markers <- subset(avg.markers, avg.markers$p_val_adj < 0.05)
avg.markers$gene <- rownames(avg.markers)
avg.markers <- dplyr :: arrange(avg.markers, desc(avg.markers$avg_log2FC))

markers.to.plot <- c("Tnfrsf1b","Nfkb1","Acod1","Trem1","Clic1","Prg4")

pdf("./Output_new4_30/Try4_Myel_a_Ctr_CIA_JAK1i/avg.Myel_a1_CIAvsJAK1i.pdf", useDingbats = F, height = 3.3, width = 3)
p1 <- ggplot(avg.Myel_a, aes(B_CIA, C_JAK1i)) + geom_point(size=0.1) +
  geom_point(data=avg.Myel_a[avg.markers$gene, ], aes(B_CIA, C_JAK1i), size=0.1, colour="red") +
  #geom_point(data=avg.Myel_b[markers.to.plot, ], aes(B_CIA, C_JAK1i), colour="red") +
  ggtitle("Myel_a1") + theme(text = element_text(size=15))
LabelPoints(plot = p1, points = markers.to.plot, repel = TRUE, color="blue")
dev.off()

#
Idents(Myel_a) <- "clusters.group"
cluster.averages <- AverageExpression(Myel_a, return.seurat = T)
gene_for_heatmap <- c(head(avg.markers$gene,25), tail(avg.markers$gene,25))
pdf("./Output_new4_30/Try4_Myel_a_Ctr_CIA_JAK1i/avg.Myel_a1_CIAvsJAK1i_heatmap.pdf", useDingbats = F, height = 11, width = 4.5)
DoHeatmap(subset(cluster.averages, idents=c("Myel_a1_B_CIA","Myel_a1_C_JAK1i")), features = gene_for_heatmap, size = 3, draw.lines = F, raster =F,
          group.colors = c("#F8766D","#D39200")) + scale_fill_gradientn(colors = c("blue","white","red"))
dev.off()

# RNAseq Myel_a2
Idents(Myel_a) <- "saved.idents"
Myel_asub <- subset(Myel_a, idents = "Myel_a2")
Idents(Myel_asub) <- "group"
avg.Myel_a <- data.frame(log1p(AverageExpression(Myel_asub, verbose = FALSE)$RNA))
avg.Myel_a$gene <- rownames(avg.Myel_a)
avg.markers <- FindMarkers(Myel_asub, ident.1 = "B_CIA", ident.2 = "C_JAK1i",
                           only.pos = F, min.pct = 0.1, logfc.threshold = 0.25, verbose = FALSE)
avg.markers <- subset(avg.markers, avg.markers$p_val_adj < 0.05)
avg.markers$gene <- rownames(avg.markers)
avg.markers <- dplyr :: arrange(avg.markers, desc(avg.markers$avg_log2FC))

markers.to.plot <- c("Plaur","Tnfaip2","Cd14","Trem1","Fos","Cxcr2","Lrg1")

pdf("./Output_new4_30/Try4_Myel_a_Ctr_CIA_JAK1i/avg.Myel_a2_CIAvsJAK1i.pdf", useDingbats = F, height = 3.3, width = 3)
p1 <- ggplot(avg.Myel_a, aes(B_CIA, C_JAK1i)) + geom_point(size=0.1) +
  geom_point(data=avg.Myel_a[avg.markers$gene, ], aes(B_CIA, C_JAK1i), size=0.1, colour="red") +
  #geom_point(data=avg.Myel_b[markers.to.plot, ], aes(B_CIA, C_JAK1i), colour="red") +
  ggtitle("Myel_a2") + theme(text = element_text(size=15))
LabelPoints(plot = p1, points = markers.to.plot, repel = TRUE, color="blue")
dev.off()

#
Idents(Myel_a) <- "clusters.group"
cluster.averages <- AverageExpression(Myel_a, return.seurat = T)
gene_for_heatmap <- c(head(avg.markers$gene,25), tail(avg.markers$gene,21))
pdf("./Output_new4_30/Try4_Myel_a_Ctr_CIA_JAK1i/avg.Myel_a2_CIAvsJAK1i_heatmap.pdf", useDingbats = F, height = 10, width = 4.5)
DoHeatmap(subset(cluster.averages, idents=c("Myel_a2_B_CIA","Myel_a2_C_JAK1i")), features = gene_for_heatmap, size = 3, draw.lines = F, raster =F,
          group.colors = c("#F8766D","#D39200")) + scale_fill_gradientn(colors = c("blue","white","red"))
dev.off()

# RNAseq Myel_a3
Idents(Myel_a) <- "saved.idents"
Myel_asub <- subset(Myel_a, idents = "Myel_a3")
Idents(Myel_asub) <- "group"
avg.Myel_a <- data.frame(log1p(AverageExpression(Myel_asub, verbose = FALSE)$RNA))
avg.Myel_a$gene <- rownames(avg.Myel_a)
avg.markers <- FindMarkers(Myel_asub, ident.1 = "B_CIA", ident.2 = "C_JAK1i",
                           only.pos = F, min.pct = 0.1, logfc.threshold = 0.25, verbose = FALSE)
avg.markers <- subset(avg.markers, avg.markers$p_val_adj < 0.05)
avg.markers$gene <- rownames(avg.markers)
avg.markers <- dplyr :: arrange(avg.markers, desc(avg.markers$avg_log2FC))

markers.to.plot <- c("Cstb","Acod1","Trem1","Junb","Prg4","Lrg1")

pdf("./Output_new4_30/Try4_Myel_a_Ctr_CIA_JAK1i/avg.Myel_a3_CIAvsJAK1i.pdf", useDingbats = F, height = 3.3, width = 3)
p1 <- ggplot(avg.Myel_a, aes(B_CIA, C_JAK1i)) + geom_point(size=0.1) +
  geom_point(data=avg.Myel_a[avg.markers$gene, ], aes(B_CIA, C_JAK1i), size=0.1, colour="red") +
  #geom_point(data=avg.Myel_b[markers.to.plot, ], aes(B_CIA, C_JAK1i), colour="red") +
  ggtitle("Myel_a3") + theme(text = element_text(size=15))
LabelPoints(plot = p1, points = markers.to.plot, repel = TRUE, color="blue")
dev.off()

#
Idents(Myel_a) <- "clusters.group"
cluster.averages <- AverageExpression(Myel_a, return.seurat = T)
gene_for_heatmap <- c(head(avg.markers$gene,25), tail(avg.markers$gene,25))
pdf("./Output_new4_30/Try4_Myel_a_Ctr_CIA_JAK1i/avg.Myel_a3_CIAvsJAK1i_heatmap.pdf", useDingbats = F, height = 11, width = 4.5)
DoHeatmap(subset(cluster.averages, idents=c("Myel_a3_B_CIA","Myel_a3_C_JAK1i")), features = gene_for_heatmap, size = 3, draw.lines = F, raster =F,
          group.colors = c("#F8766D","#D39200")) + scale_fill_gradientn(colors = c("blue","white","red"))
dev.off()

