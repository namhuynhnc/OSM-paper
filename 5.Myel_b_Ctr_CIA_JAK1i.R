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
library(clusterExperiment)
library(ggrepel)

set.seed(1234567)

options(ggrepel.max.overlaps = Inf)

# Load data
merge.integrated_3 <- readRDS("./Save_new/seuset_merge_Ctr_CIA_JAK1i_30_2.RDS")

# Select myeloid b
Myel_b <- subset(merge.integrated_3, idents = c("Myeloid_b", "Myeloid_c"))
Myel_b

# DefaultAssay
DefaultAssay(object = Myel_b) <- "RNA"
merge.list <- SplitObject(object = Myel_b, split.by = "group")

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
saveRDS(merge.integrated_1, "./Save_new/Myel_b_Ctr_CIA_JAK1i_30_1.RDS")
merge.integrated_1 <- readRDS("./Save_new/Myel_b_Ctr_CIA_JAK1i_30_1.RDS")

# Clustering
Myel_b <- FindClusters(merge.integrated_1, verbose = FALSE, resolution = 0.5)
table(Myel_b@active.ident)

Myel_b <- RenameIdents(object=Myel_b,"5"="Myel_c5")
Myel_b <- RenameIdents(object=Myel_b,"0"="Myel_c4")
Myel_b <- RenameIdents(object=Myel_b,"3"="Myel_c3")
Myel_b <- RenameIdents(object=Myel_b,"1"="Myel_c2")
Myel_b <- RenameIdents(object=Myel_b,"2"="Myel_c1")
Myel_b <- RenameIdents(object=Myel_b,"4"="Myel_b")

Myel_b$saved.idents <- Idents(object = Myel_b)

# DefaultAssay
DefaultAssay(object = Myel_b) <- "RNA"
Myel_b <- NormalizeData(Myel_b, verbose = FALSE)

pdf("./Output_new4_30/Try5_Myel_b_Ctr_CIA_JAK1i/dimplot_0.5.pdf", useDingbats = F, height = 4, width = 6)
DimPlot(Myel_b, label = TRUE, label.size = 4) + NoLegend()
DimPlot(Myel_b, label = TRUE, label.size = 4, group.by = "group")
dev.off()

DimPlot(Myel_b, label = TRUE, label.size = 4, split.by = "group")


markers.to.plot <- c("Cd14","Ccr2","Ly6c2","Pdpn",
                    "Mmp19","Malt1","Msr1","Lgmn","Tnfrsf11a","Cx3cr1","Mrc1","Itgam","Adgre1","Csf1r","Fcgr3","Mertk","Itgax","Cd209a")

pdf("./Output_new4_30/Try5_Myel_b_Ctr_CIA_JAK1i/featureplot_1.pdf", useDingbats = F, height = 20, width = 15)
FeaturePlot(Myel_b, markers.to.plot, pt.size = 0.5,
            cols = c("lightgrey", "#005493")
              , min.cutoff = "q10", max.cutoff = "q90"
)
dev.off()

#
markers.to.plot <- c("Mertk","Folr2","Trem2","Nupr1","Plaur","Hbegf","Cx3cr1","Lyve1","Vegfa",
                     "Tnf","Il1b","Il6","Osm","Acp5","Ctsk","Nfatc1","Tnfrsf11a","Csf1r",
                     "Ptges2","Egfr","Tgfbi","Itgb5","Adora3","Isg15","Ifitm1",
                     "Cxcl2","Cxcl10","Il10")

pdf("./Output_new4_30/Try5_Myel_b_Ctr_CIA_JAK1i/Dotplot.pdf", useDingbats = F, height = 6, width = 14)
DotPlot(Myel_b, features = markers.to.plot, dot.scale = 8, cols = c("blue","blue","blue"),
        split.by = "group") + RotatedAxis()
dev.off()

# Markers
markers <- FindAllMarkers(object = Myel_b, only.pos = F, min.pct = 0.1, logfc.threshold = 0.1)
markers <- subset(markers, markers$p_val_adj < 0.05)
markers <- dplyr ::arrange(markers, desc(markers$avg_log2FC))
markers <- dplyr ::arrange(markers, markers$cluster)
write.table(markers,"./Output_new4_30/Try5_Myel_b_Ctr_CIA_JAK1i/markers_0.5.txt", sep="\t", quote=FALSE, row.names=FALSE)

# RNAseq Myel_b
Idents(Myel_b) <- "saved.idents"
Myel_bsub <- subset(Myel_b, idents = "Myel_b")
Idents(Myel_bsub) <- "group"
avg.Myel_b <- data.frame(log1p(AverageExpression(Myel_bsub, verbose = FALSE)$RNA))
avg.Myel_b$gene <- rownames(avg.Myel_b)
avg.markers <- FindMarkers(Myel_bsub, ident.1 = "B_CIA", ident.2 = "C_JAK1i",
                           only.pos = F, min.pct = 0.1, logfc.threshold = 0.25, verbose = FALSE)
avg.markers <- subset(avg.markers, avg.markers$p_val_adj < 0.05)
avg.markers$gene <- rownames(avg.markers)
avg.markers <- dplyr :: arrange(avg.markers, desc(avg.markers$avg_log2FC))

markers.to.plot <- c("Cd14","Cxcr4","Lilr4b","Fosb","Gm42418")

pdf("./Output_new4_30/Try5_Myel_b_Ctr_CIA_JAK1i/avg.Myel_b_CIAvsJAK1i.pdf", useDingbats = F, height = 5, width = 5)
p1 <- ggplot(avg.Myel_b, aes(B_CIA, C_JAK1i)) + geom_point(size=0.1) +
  geom_point(data=avg.Myel_b[avg.markers$gene, ], aes(B_CIA, C_JAK1i), size=0.1, colour="red") +
 # geom_point(data=avg.Myel_b[markers.to.plot, ], aes(B_CIA, C_JAK1i), colour="red") +
  ggtitle("Myel_b")
LabelPoints(plot = p1, points = markers.to.plot, repel = TRUE, color="blue")
dev.off()

#
cluster.averages <- AverageExpression(Myel_bsub, return.seurat = T)
gene_for_heatmap <- c(head(avg.markers$gene,25), tail(avg.markers$gene,11))
pdf("./Output_new4_30/Try5_Myel_b_Ctr_CIA_JAK1i/avg.Myel_b_CIAvsJAK1i_heatmap.pdf", useDingbats = F, height = 8, width = 4)
DoHeatmap(subset(cluster.averages, idents = c("B_CIA","C_JAK1i")), features = gene_for_heatmap, size = 3, draw.lines = F, raster =F,
          group.colors = c("#F8766D","#D39200")) + scale_fill_gradientn(colors = c("blue","white","red"))
dev.off()

# RNAseq Myel_c1
Idents(Myel_b) <- "saved.idents"
Myel_bsub <- subset(Myel_b, idents = "Myel_c1")
Idents(Myel_bsub) <- "group"
avg.Myel_b <- data.frame(log1p(AverageExpression(Myel_bsub, verbose = FALSE)$RNA))
avg.Myel_b$gene <- rownames(avg.Myel_b)
avg.markers <- FindMarkers(Myel_bsub, ident.1 = "B_CIA", ident.2 = "C_JAK1i",
                           only.pos = F, min.pct = 0.1, logfc.threshold = 0.25, verbose = FALSE)
avg.markers <- subset(avg.markers, avg.markers$p_val_adj < 0.05)
avg.markers$gene <- rownames(avg.markers)
avg.markers <- dplyr :: arrange(avg.markers, desc(avg.markers$avg_log2FC))

markers.to.plot <- c("Inhba","Cxcl3","Il1b","Cd14","Cxcr4","Ets2","Nfkb1","Lilr4b","Mmp13","Nfkb2","Cxcl5","Junb","Tnf","Fosb","Prg4","Jun","Nr4a1","Fos")

pdf("./Output_new4_30/Try5_Myel_b_Ctr_CIA_JAK1i/avg.Myel_c1_CIAvsJAK1i.pdf", useDingbats = F, height = 3.3, width = 3)
p1 <- ggplot(avg.Myel_b, aes(B_CIA, C_JAK1i)) + geom_point(size=0.1) +
  geom_point(data=avg.Myel_b[avg.markers$gene, ], aes(B_CIA, C_JAK1i), size=0.1, colour="red") +
  #geom_point(data=avg.Myel_b[markers.to.plot, ], aes(B_CIA, C_JAK1i), colour="red") +
  ggtitle("Myel_c1") + theme(text = element_text(size=15))
LabelPoints(plot = p1, points = markers.to.plot, repel = TRUE, color="blue")
dev.off()

#
cluster.averages <- AverageExpression(Myel_bsub, return.seurat = T)
gene_for_heatmap <- c(head(avg.markers$gene,25), tail(avg.markers$gene,25))
pdf("./Output_new4_30/Try5_Myel_b_Ctr_CIA_JAK1i/avg.Myel_c1_CIAvsJAK1i_heatmap.pdf", useDingbats = F, height = 10, width = 4)
DoHeatmap(subset(cluster.averages, idents = c("B_CIA","C_JAK1i")), features = gene_for_heatmap, size = 3, draw.lines = F, raster =F,
          group.colors = c("#F8766D","#D39200")) + scale_fill_gradientn(colors = c("blue","white","red"))
dev.off()

# RNAseq Myel_c2
Idents(Myel_b) <- "saved.idents"
Myel_bsub <- subset(Myel_b, idents = "Myel_c2")
Idents(Myel_bsub) <- "group"
avg.Myel_b <- data.frame(log1p(AverageExpression(Myel_bsub, verbose = FALSE)$RNA))
avg.Myel_b$gene <- rownames(avg.Myel_b)
avg.markers <- FindMarkers(Myel_bsub, ident.1 = "B_CIA", ident.2 = "C_JAK1i",
                           only.pos = F, min.pct = 0.1, logfc.threshold = 0.25, verbose = FALSE)
avg.markers <- subset(avg.markers, avg.markers$p_val_adj < 0.05)
avg.markers$gene <- rownames(avg.markers)
avg.markers <- dplyr :: arrange(avg.markers, desc(avg.markers$avg_log2FC))

markers.to.plot <- c("Il1b","Cd14","Ccrl2","Cd14","Cxcr4","Mmp19","Osm","Mmp3","Pdpn","Mmp13","Prg4","Cx3cr1","Egr1","Prg4")

pdf("./Output_new4_30/Try5_Myel_b_Ctr_CIA_JAK1i/avg.Myel_c2_CIAvsJAK1i.pdf", useDingbats = F, height = 3.3, width = 3)
p1 <- ggplot(avg.Myel_b, aes(B_CIA, C_JAK1i)) + geom_point(size=0.1) +
  geom_point(data=avg.Myel_b[avg.markers$gene, ], aes(B_CIA, C_JAK1i), size=0.1, colour="red") +
  #geom_point(data=avg.Myel_b[markers.to.plot, ], aes(B_CIA, C_JAK1i), colour="red") +
  ggtitle("Myel_c2") + theme(text = element_text(size=15))
LabelPoints(plot = p1, points = markers.to.plot, repel = TRUE, color="blue")
dev.off()

#
cluster.averages <- AverageExpression(Myel_bsub, return.seurat = T)
gene_for_heatmap <- c(head(avg.markers$gene,25), tail(avg.markers$gene,25))
pdf("./Output_new4_30/Try5_Myel_b_Ctr_CIA_JAK1i/avg.Myel_c2_CIAvsJAK1i_heatmap.pdf", useDingbats = F, height = 10, width = 4)
DoHeatmap(subset(cluster.averages, idents = c("B_CIA","C_JAK1i")), features = gene_for_heatmap, size = 3, draw.lines = F, raster =F,
          group.colors = c("#F8766D","#D39200")) + scale_fill_gradientn(colors = c("blue","white","red"))
dev.off()

# RNAseq Myel_c3
Idents(Myel_b) <- "saved.idents"
Myel_bsub <- subset(Myel_b, idents = "Myel_c3")
Idents(Myel_bsub) <- "group"
avg.Myel_b <- data.frame(log1p(AverageExpression(Myel_bsub, verbose = FALSE)$RNA))
avg.Myel_b$gene <- rownames(avg.Myel_b)
avg.markers <- FindMarkers(Myel_bsub, ident.1 = "B_CIA", ident.2 = "C_JAK1i",
                           only.pos = F, min.pct = 0.1, logfc.threshold = 0.25, verbose = FALSE)
avg.markers <- subset(avg.markers, avg.markers$p_val_adj < 0.05)
avg.markers$gene <- rownames(avg.markers)
avg.markers <- dplyr :: arrange(avg.markers, desc(avg.markers$avg_log2FC))

markers.to.plot <- c("Ets2","Spp1","Cxcl3","Acp5","C1qa","C1qb")

pdf("./Output_new4_30/Try5_Myel_b_Ctr_CIA_JAK1i/avg.Myel_c3_CIAvsJAK1i.pdf", useDingbats = F, height = 3.3, width = 3)
p1 <- ggplot(avg.Myel_b, aes(B_CIA, C_JAK1i)) + geom_point(size=0.1) +
  geom_point(data=avg.Myel_b[avg.markers$gene, ], aes(B_CIA, C_JAK1i), size=0.1, colour="red") +
  geom_point(data=avg.Myel_b[markers.to.plot, ], aes(B_CIA, C_JAK1i), size=0.1, colour="red") +
  ggtitle("Myel_c3") + theme(text = element_text(size=15))
LabelPoints(plot = p1, points = markers.to.plot, repel = TRUE, color="blue")
dev.off()

#
cluster.averages <- AverageExpression(Myel_bsub, return.seurat = T)
gene_for_heatmap <- c(head(avg.markers$gene,25), tail(avg.markers$gene,8))
pdf("./Output_new4_30/Try5_Myel_b_Ctr_CIA_JAK1i/avg.Myel_c3_CIAvsJAK1i_heatmap.pdf", useDingbats = F, height = 7, width = 4)
DoHeatmap(subset(cluster.averages, idents = c("B_CIA","C_JAK1i")), features = gene_for_heatmap, size = 3, draw.lines = F, raster =F,
          group.colors = c("#F8766D","#D39200")) + scale_fill_gradientn(colors = c("blue","white","red"))
dev.off()

# RNAseq Myel_c4
Idents(Myel_b) <- "saved.idents"
Myel_bsub <- subset(Myel_b, idents = "Myel_c4")
Idents(Myel_bsub) <- "group"
avg.Myel_b <- data.frame(log1p(AverageExpression(Myel_bsub, verbose = FALSE)$RNA))
avg.Myel_b$gene <- rownames(avg.Myel_b)
avg.markers <- FindMarkers(Myel_bsub, ident.1 = "B_CIA", ident.2 = "C_JAK1i",
                           only.pos = F, min.pct = 0.1, logfc.threshold = 0.25, verbose = FALSE)
avg.markers <- subset(avg.markers, avg.markers$p_val_adj < 0.05)
avg.markers$gene <- rownames(avg.markers)
avg.markers <- dplyr :: arrange(avg.markers, desc(avg.markers$avg_log2FC))

markers.to.plot <- c("Ccl8","Ifitm1","Ccrl2","Cxcr4","Vcam1","Mmp19","Cd14","Il1b","Osm","Mmp13","Lyz2","Prg4","Cx3cr1")

pdf("./Output_new4_30/Try5_Myel_b_Ctr_CIA_JAK1i/avg.Myel_c4_CIAvsJAK1i.pdf", useDingbats = F, height = 3.3, width = 3)
p1 <- ggplot(avg.Myel_b, aes(B_CIA, C_JAK1i)) + geom_point(size=0.1) +
  geom_point(data=avg.Myel_b[avg.markers$gene, ], aes(B_CIA, C_JAK1i), size=0.1, colour="red") +
 # geom_point(data=avg.Myel_b[markers.to.plot, ], aes(B_CIA, C_JAK1i), colour="red") +
  ggtitle("Myel_c4") + theme(text = element_text(size=15))
LabelPoints(plot = p1, points = markers.to.plot, repel = TRUE, color="blue")
dev.off()

#
cluster.averages <- AverageExpression(Myel_bsub, return.seurat = T)
gene_for_heatmap <- c(head(avg.markers$gene,25), tail(avg.markers$gene,25))
pdf("./Output_new4_30/Try5_Myel_b_Ctr_CIA_JAK1i/avg.Myel_c4_CIAvsJAK1i_heatmap.pdf", useDingbats = F, height = 10, width = 4)
DoHeatmap(subset(cluster.averages, idents = c("B_CIA","C_JAK1i")), features = gene_for_heatmap, size = 3, draw.lines = F, raster =F,
          group.colors = c("#F8766D","#D39200")) + scale_fill_gradientn(colors = c("blue","white","red"))
dev.off()

# RNAseq Myel_c5
Idents(Myel_b) <- "saved.idents"
Myel_bsub <- subset(Myel_b, idents = "Myel_c5")
Idents(Myel_bsub) <- "group"
avg.Myel_b <- data.frame(log1p(AverageExpression(Myel_bsub, verbose = FALSE)$RNA))
avg.Myel_b$gene <- rownames(avg.Myel_b)
avg.markers <- FindMarkers(Myel_bsub, ident.1 = "B_CIA", ident.2 = "C_JAK1i",
                           only.pos = F, min.pct = 0.1, logfc.threshold = 0.25, verbose = FALSE)
avg.markers <- subset(avg.markers, avg.markers$p_val_adj < 0.05)
avg.markers$gene <- rownames(avg.markers)
avg.markers <- dplyr :: arrange(avg.markers, desc(avg.markers$avg_log2FC))

markers.to.plot <- c("Nme2")

pdf("./Output_new4_30/Try5_Myel_b_Ctr_CIA_JAK1i/avg.Myel_c5_CIAvsJAK1i.pdf", useDingbats = F, height = 3.3, width = 3)
p1 <- ggplot(avg.Myel_b, aes(B_CIA, C_JAK1i)) + geom_point(size=0.1) +
  geom_point(data=avg.Myel_b[avg.markers$gene, ], aes(B_CIA, C_JAK1i), size=0.1, colour="red") +
 # geom_point(data=avg.Myel_b[markers.to.plot, ], aes(B_CIA, C_JAK1i), colour="red") +
  ggtitle("Myel_c5") + theme(text = element_text(size=15))
LabelPoints(plot = p1, points = markers.to.plot, repel = TRUE, color="blue")
dev.off()
