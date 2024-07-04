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
library(ggplot2)
library(multtest)
library(Seurat)
library(R.utils)
library(rhdf5)
library(Ecfun)
library(sctransform)
#library(scater)
#library(scran)
library(SingleCellExperiment)
library(gplots)
library(ggbeeswarm)
library(ggthemes)
library(slingshot)
library(gam)
library(viridis)
library(autoimage)
#library(clusterExperiment)

set.seed(1234567)

# Read Ctr
C_JAK1i<- Read10X(data.dir = "./Data/Other/01_analysis/cellranger_count/NK003-C/filtered_feature_bc_matrix")
seusetB_CIA <- CreateSeuratObject(counts = C_JAK1i, min.cells = 3, min.features = 200)
seusetB_CIA$group <- "A_Ctr"

# Read CIA
CIA<- Read10X(data.dir = "./Data/Other/01_analysis/cellranger_count/NK004-ART/filtered_feature_bc_matrix")
seuset_CIA <- CreateSeuratObject(counts = CIA, min.cells = 3, min.features = 200)
seuset_CIA$group <- "B_CIA"

# Read JAK1i
JAK1i<- Read10X(data.dir = "./Data/Other/01_analysis/cellranger_count/L5matrix_30/outs/filtered_feature_bc_matrix")
seuset_JAK1i <- CreateSeuratObject(counts = JAK1i, min.cells = 3, min.features = 200)
seuset_JAK1i$group <- "C_JAK1i"

# Merge
merge <- merge(x= seusetB_CIA, y=c(seuset_CIA,seuset_JAK1i))
merge

# QC and selecting cells for further analysis
merge <- PercentageFeatureSet(merge, pattern = "^mt-", col.name = "percent.mt")
#Expression QC
VlnPlot(object = merge, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
FeatureScatter(object = merge, feature1 = "nCount_RNA", feature2 = "percent.mt")
FeatureScatter(object = merge, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

# We filter out cells
merge <- subset(x = merge, subset = nFeature_RNA < 7500 & nCount_RNA<75000 & percent.mt < 5)
merge

# SCT
merge.list <- SplitObject(object = merge, split.by = "group")
for (i in 1:length(x = merge.list)) {
  merge.list[[i]] <- SCTransform(object = merge.list[[i]], verbose = FALSE
                                 , vars.to.regress = "percent.mt"
  )
}
merge.features <- SelectIntegrationFeatures(object.list = merge.list, nfeatures = 3000)
merge.list <- PrepSCTIntegration(object.list = merge.list, anchor.features = merge.features, verbose = FALSE)
merge.anchors <- FindIntegrationAnchors(object.list = merge.list, normalization.method = "SCT", 
                                        anchor.features = merge.features, verbose = FALSE
                                        )
merge.integrated <- IntegrateData(anchorset = merge.anchors, normalization.method = "SCT", verbose = FALSE)

# PCA, UMAP
merge.integrated <- RunPCA(object = merge.integrated, verbose = FALSE)
merge.integrated <- RunUMAP(object = merge.integrated, dims = 1:30)
merge.integrated <- FindNeighbors(object = merge.integrated, dims = 1:30)
merge.integrated_1 <- FindClusters(merge.integrated, verbose = FALSE, resolution = 1.5)
DimPlot(merge.integrated_1, label = TRUE, label.size = 2)

# DefaultAssay
DefaultAssay(object = merge.integrated_1) <- "RNA"
# Normalize RNA data for visualization purposes
merge.integrated_1 <- NormalizeData(merge.integrated_1, verbose = FALSE)
#
markers.to.plot <- c("Ptprc","Cd3g","Cd19","Cd14","Pdpn","Fap","Ngp","Camp","Pecam1","Mcam","Notch3")
DotPlot(merge.integrated_1, features = markers.to.plot, dot.scale = 8) + RotatedAxis()

# Remove minor cells (Cycling_cell, Neural_cell, Erythrocyte #17, 19, 21)
merge.integrated_2 <- subset(merge.integrated_1, idents = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,18,20,22,23))

merge.integrated_2
merge.integrated_2$saved.idents <- Idents(object = merge.integrated_2)

# DefaultAssay
DefaultAssay(object = merge.integrated_2) <- "RNA"
merge.list <- SplitObject(object = merge.integrated_2, split.by = "group")

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
merge.integrated <- RunPCA(object = merge.integrated, verbose = FALSE)
merge.integrated <- RunUMAP(object = merge.integrated, dims = 1:30)
merge.integrated <- FindNeighbors(object = merge.integrated, dims = 1:30)

merge.integrated_3 <- FindClusters(merge.integrated, verbose = FALSE, resolution = 0.1)
table(merge.integrated_3@active.ident)

# DefaultAssay
DefaultAssay(object = merge.integrated_3) <- "RNA"
# Normalize RNA data for visualization purposes
merge.integrated_3 <- NormalizeData(merge.integrated_3, verbose = FALSE)

merge.integrated_3 <- RenameIdents(object=merge.integrated_3,"8"="Mural_cell")
merge.integrated_3 <- RenameIdents(object=merge.integrated_3,"4"="Endothelial")
merge.integrated_3 <- RenameIdents(object=merge.integrated_3,"3"="Neutrophil")

merge.integrated_3 <- RenameIdents(object=merge.integrated_3,"1"="Fibroblast")

merge.integrated_3 <- RenameIdents(object=merge.integrated_3,"0"="Myeloid_c")
merge.integrated_3 <- RenameIdents(object=merge.integrated_3,"5"="Myeloid_b")
merge.integrated_3 <- RenameIdents(object=merge.integrated_3,"2"="Myeloid_a")

merge.integrated_3 <- RenameIdents(object=merge.integrated_3,"7"="B_cell")
merge.integrated_3 <- RenameIdents(object=merge.integrated_3,"6"="T_cell")

# plot
pdf("./Output_new4_30/Try13_merge_Ctr_CIA_JAK1i/dimplot_0.1.pdf", useDingbats = F, height = 4, width = 6)
DimPlot(merge.integrated_3, label = TRUE, label.size = 4) + NoLegend()
DimPlot(merge.integrated_3, label = TRUE, label.size = 0,
        order=c("A_Ctr", "B_CIA", "C_JAK1i"), group.by = "group")
dev.off()

pdf("./Output_new4_30/Try13_merge_Ctr_CIA_JAK1i/dimplot_0.1b.pdf", useDingbats = F, height = 4, width = 14)
DimPlot(merge.integrated_3, label = TRUE, label.size = 0, split.by = "group")
dev.off()

# Markers
markers <- FindAllMarkers(object = merge.integrated_3, only.pos = F, min.pct = 0.1, logfc.threshold = 0.1)
markers <- subset(markers, markers$p_val_adj < 0.05)
markers <- dplyr ::arrange(markers, desc(markers$avg_log2FC))
markers <- dplyr ::arrange(markers, markers$cluster)
write.table(markers,"./Output_new4_30/Try13_merge_Ctr_CIA_JAK1i/markers_0.1.txt", sep="\t", quote=FALSE, row.names=FALSE)

#Plots
markers.to.plot <- c("Ptprc","Cd3g","Cd19","Cd14","Ccl4","Il1b","Tnf","Ccr2","Ifitm3","C1qb","Lgmn","Pdpn","Fap","Ngp","Camp","Pecam1","Mcam","Notch3")
pdf("./Output_new4_30/Try13_merge_Ctr_CIA_JAK1i/dotplot_0.1.pdf", useDingbats = F, height = 4, width = 8)
DotPlot(merge.integrated_3, features = markers.to.plot, dot.scale = 8) + RotatedAxis()
dev.off()

pdf("./Output_new4_30/Try13_merge_Ctr_CIA_JAK1i/dotplot_0.1b.pdf", useDingbats = F, height = 8, width = 9)
DotPlot(merge.integrated_3, features = markers.to.plot, cols = c("#f2756d", "#058992", "#c78d32"), dot.scale = 8,
        split.by = "group") + RotatedAxis()

dev.off()

#
DimPlot(merge.integrated_3, label = TRUE, label.size = 2) + NoLegend()
merge.integrated_3$clusters.group <- paste(Idents(merge.integrated_3), merge.integrated_3$group, sep = "_")
merge.integrated_3$clusters <- Idents(merge.integrated_3)
clusters <- as.data.frame(table(merge.integrated_3$clusters))
## RNAseq all cells
RNAseq <- merge.integrated_3
Idents(RNAseq) <- "group"
avg <- data.frame(log1p(AverageExpression(RNAseq, verbose = FALSE)$RNA))
avg$gene <- rownames(avg)

# CIA vs JAK1i
avg.markers <- FindMarkers(RNAseq, ident.1 = "B_CIA", ident.2 = "C_JAK1i",
                           only.pos = F, min.pct = 0.1, logfc.threshold = 0.25, verbose = FALSE)
avg.markers <- subset(avg.markers, avg.markers$p_val_adj < 0.05)
avg.markers$gene <- rownames(avg.markers)
avg.markers <- dplyr :: arrange(avg.markers, desc(avg.markers$avg_log2FC))
write.table(avg.markers,"./Output_new4_30/Try13_merge_Ctr_CIA_JAK1i/markers_CIAvsJAK1i.txt", sep="\t", quote=FALSE, row.names=F)

markers.to.plot <- c("Ccrl2","Cd14","Ccl3","Ccl4","Cxcl3","Ets2","Tnfrsf1b","Cxcl2","Ccr1","Cxcr4","Cxcl5","Mmp3",
                     "Mmp13","Ifitm1","Nfkb1","Il1a","Il6","Stat3",
                     "Ly6a","Cx3cr1","Cxcl12","Ctsk","Ly6c1","Ecm1","Clic5","C1qb","Lyz1","Hbegf","Dcn","Prg4")

match(markers.to.plot, avg.markers$gene)

pdf("./Output_new4_30/Try13_merge_Ctr_CIA_JAK1i/avg.CIAvsJAK1i.pdf", useDingbats = F, height = 3, width = 3)
p1 <- ggplot(avg, aes(B_CIA, C_JAK1i)) + geom_point(size=0.1) +
  geom_point(data=avg[avg.markers$gene, ], aes(B_CIA, C_JAK1i), size=0.1, colour="red") +
#  geom_point(data=avg[markers.to.plot, ], aes(B_CIA, C_JAK1i), colour="red") +
  #ggtitle("Synovial cells") + 
  theme(text = element_text(size=15))
LabelPoints(plot = p1, points = markers.to.plot, repel = TRUE, color="blue")
dev.off()

## Assign cluster
merge.integrated_4 <- readRDS("./Save_new/seuset_merge_Ctr_CIA_JAK1i_30_2.RDS")
DimPlot(merge.integrated_4, label = TRUE, label.size = 2) + NoLegend()

Idents(merge.integrated_4, cells= WhichCells(merge.integrated_4, idents = c("Mural_cell"))) <- "Mural_cell"
Idents(merge.integrated_4, cells= WhichCells(merge.integrated_4, idents = c("Endothelial"))) <- "Endothelial"
Idents(merge.integrated_4, cells= WhichCells(merge.integrated_4, idents = c("Neutrophil"))) <- "Neutrophil"

# Fibroblast
Fibro <- readRDS("./Save_new/Fibro_Ctr_CIA_JAK1i_30_2.RDS")
Idents(merge.integrated_4, cells= WhichCells(Fibro, idents = c("Chondrocyte"))) <- "Chondrocyte"
Idents(merge.integrated_4, cells= WhichCells(Fibro, idents = c("Fibro_4"))) <- "Fibro_4"
Idents(merge.integrated_4, cells= WhichCells(Fibro, idents = c("Fibro_3"))) <- "Fibro_3"
Idents(merge.integrated_4, cells= WhichCells(Fibro, idents = c("Fibro_2"))) <- "Fibro_2"
Idents(merge.integrated_4, cells= WhichCells(Fibro, idents = c("Fibro_1"))) <- "Fibro_1"

# Myeloid b
Myel_b <- readRDS("./Save_new/Myel_b_aggr_Ctr_CIA_30_2.RDS")
Idents(merge.integrated_4, cells= WhichCells(Myel_b, idents = c("Myel_c5"))) <- "Myel_c5"
Idents(merge.integrated_4, cells= WhichCells(Myel_b, idents = c("Myel_c4"))) <- "Myel_c4"
Idents(merge.integrated_4, cells= WhichCells(Myel_b, idents = c("Myel_c3"))) <- "Myel_c3"
Idents(merge.integrated_4, cells= WhichCells(Myel_b, idents = c("Myel_c2"))) <- "Myel_c2"
Idents(merge.integrated_4, cells= WhichCells(Myel_b, idents = c("Myel_c1"))) <- "Myel_c1"
Idents(merge.integrated_4, cells= WhichCells(Myel_b, idents = c("Myel_b"))) <- "Myel_b"

# Myeloid a
Myel_a <- readRDS("./Save_new/Myel_a_Ctr_CIA_JAK1i_30_2.RDS")
Idents(merge.integrated_4, cells= WhichCells(Myel_a, idents = c("Myel_a3"))) <- "Myel_a3"
Idents(merge.integrated_4, cells= WhichCells(Myel_a, idents = c("Myel_a2"))) <- "Myel_a2"
Idents(merge.integrated_4, cells= WhichCells(Myel_a, idents = c("Myel_a1"))) <- "Myel_a1"

Idents(merge.integrated_4, cells= WhichCells(merge.integrated_4, idents = c("B_cell"))) <- "B_cell"
Idents(merge.integrated_4, cells= WhichCells(merge.integrated_4, idents = c("T_cell"))) <- "T_cell"

# plot
pdf("./Output_new4_30/Try13_merge_Ctr_CIA_JAK1i/dimplot_celltype.pdf", useDingbats = F, height = 4, width = 6)
DimPlot(merge.integrated_4, label = TRUE, label.size = 2) + NoLegend()
DimPlot(merge.integrated_4, label = TRUE, label.size = 0, group.by = "group")
dev.off()

pdf("./Output_new4_30/Try13_merge_Ctr_CIA_JAK1i/dimplot_celltypeb.pdf", useDingbats = F, height = 4.5, width = 14)
DimPlot(merge.integrated_4, label = TRUE, label.size = 2, split.by = "group") + NoLegend()
dev.off()

#
pdf("./Output_new4_30/Try13_merge_Ctr_CIA_JAK1i/dotplot_celltypeb.pdf", useDingbats = F, height = 15, width = 15)
DotPlot(merge.integrated_4, features = markers.to.plot, cols = c("#f2756d", "#058992", "#c78d32"), dot.scale = 8,
        split.by = "group") + RotatedAxis()
dev.off()
markers.to.plot <- c("Il17a","Il1b","Tnf","Il6","Osm","Lif","Tnfsf11",'Tnfrsf11a',"Foxm1","Mmp2","Mmp3","Mmp11","Mmp13","Mmp14","Mmp19",
                      "Il2ra","Il2rb","Il2rg","Il4ra","Il7r","Il15ra","Il21r","Il6ra","Il6st","Osmr","Lifr",
                     "Ifnar1","Ifnar2","Il10ra","Il10rb","Il22ra2","Ifngr1","Ifngr2","Csf3r","Epor","Ghr",
                     "Jak1","Jak2","Jak3","Tyk2","Stat1","Stat2","Stat3","Stat5a")
pdf("./Output_new4_30/Try13_merge_Ctr_CIA_JAK1i/dotplot_celltypec.pdf", useDingbats = F, height = 6, width = 18)
DotPlot(merge.integrated_4, features = markers.to.plot, dot.scale = 8) + RotatedAxis()
dev.off()

#
markers.to.plot2 <- c("Il17a","Il1b","Tnf","Osm","Il6","Tnfsf11","Mmp3","Mmp13")
pdf("./Output_new4_30/Try13_merge_Ctr_CIA_JAK1i/dotplot_celltypeg.pdf", useDingbats = F, height = 16, width = 6)
DotPlot(merge.integrated_4,
        features = markers.to.plot2, cols = c("blue", "blue", "blue"), dot.scale = 8,
        split.by = "group") + RotatedAxis()
dev.off()

#
markers.to.plot <- c("Il6","Osm",
                     "Il2ra","Il2rb","Il2rg","Il4ra","Il7r","Il15ra","Il21r","Il6ra","Il6st","Osmr","Lifr",
                     "Ifnar1","Ifnar2","Il10ra","Il10rb","Il22ra2","Ifngr1","Ifngr2","Csf3r","Epor","Ghr")
pdf("./Output_new4_30/Try13_merge_Ctr_CIA_JAK1i/dotplot_celltypeh.pdf", useDingbats = F, height = 4, width = 9)
DotPlot(merge.integrated_4, features = markers.to.plot, dot.scale = 8, idents = c("Myel_b","Myel_c1","Myel_c2","Myel_c3","Myel_c4",
                                                                                  "Fibro_1","Fibro_2","Fibro_3")) + RotatedAxis()
dev.off()

#
markers.to.plot <- c("Itgam","Fcgr1","Adgre1","Mertk","Mrc1","Fcgr3","Tnf","Il1b","Csf1r","Ccr2","Trem2","Cx3cr1","Tnfrsf11a","Arg1","Folr2","Lyve1","Id2",
                     "H2-Aa","Itgax","Retnla")

#
pdf("./Output_new4_30/Try13_merge_Ctr_CIA_JAK1i/dotplot_celltypei.pdf", useDingbats = F, height = 4, width = 8)
DotPlot(merge.integrated_4, features = markers.to.plot, dot.scale = 8, idents = c("Myel_a1","Myel_a2","Myel_a3","Myel_b","Myel_c1",
                                                                                  "Myel_c2","Myel_c3","Myel_c4","Myel_c5")) + RotatedAxis()
dev.off()

#
markers.to.plot1 <- c("Il17a","Il1b","Tnf","Il6","Osm","Lif","Tnfsf11",'Tnfrsf11a',"Mmp3","Mmp13","Mmp14",
                      "Il2ra","Il2rb","Il2rg","Il4ra","Il7r","Il15ra","Il21r","Il6ra","Il6st","Osmr","Lifr",
                      "Ifnar1","Ifnar2","Il10ra","Il10rb","Il22ra2","Ifngr1","Ifngr2","Csf3r","Epor","Ghr",
                      "Jak1","Jak2","Jak3","Tyk2","Stat1","Stat2","Stat3","Stat5a")

pdf("./Output_new4_30/Try13_merge_Ctr_CIA_JAK1i/DotPlot_JAKSTAT.pdf", useDingbats = F, height = 12, width = 14)
DotPlot(merge.integrated_4, idents = c("T_cell","Myel_a1","Myel_a2","Myel_a3","Myel_b","Myel_c1",
                                       "Myel_c2","Myel_c3","Myel_c4","Myel_c5", "Fibro_1", "Fibro_2", "Fibro_3", "Fibro_4"), features = markers.to.plot1, cols = c("blue","blue","blue"), dot.scale = 8,
        split.by = "group") + RotatedAxis()
dev.off()

#
pdf("./Output_new4_30/Try13_merge_Ctr_CIA_JAK1i/VlnPlot_celltype3.pdf", useDingbats = F, height = 3, width = 7)
markers.to.plot <- c("Osm")
VlnPlot(merge.integrated_4, features = markers.to.plot
                 , split.by = "group" , idents = c("Myel_a1","Myel_a2","Myel_a3","Myel_b","Myel_c1",
                                                  "Myel_c2","Myel_c3","Myel_c4","Myel_c5"),
        cols = c("#078992","#F8766D","#D39200"),
                 pt.size = 0.1, combine = FALSE, split.plot = F)
dev.off()

# Markers
markers <- FindAllMarkers(object = merge.integrated_4, only.pos = F, min.pct = 0.1, logfc.threshold = 0.1)
markers <- subset(markers, markers$p_val_adj < 0.05)
markers <- dplyr ::arrange(markers, desc(markers$avg_log2FC))
markers <- dplyr ::arrange(markers, markers$cluster)
write.table(markers,"./Output_new4_30/Try13_merge_Ctr_CIA_JAK1i/markers_celltype.txt", sep="\t", quote=FALSE, row.names=FALSE)

## DE of each cluster
DimPlot(merge.integrated_4, label = TRUE, label.size = 2) + NoLegend()
merge.integrated_4$clusters.group <- paste(Idents(merge.integrated_4), merge.integrated_4$group, sep = "_")
merge.integrated_4$clusters <- Idents(merge.integrated_4)
clusters <- as.data.frame(table(merge.integrated_4$clusters))

# Each cluster Ctr_vs_CIA
plotlist <- list()
pdf("./Output_new4_30/Try13_merge_Ctr_CIA_JAK1i/marker_eachcluster_CtrvsCIA.pdf", useDingbats = F, height = 5, width = 5)
for (i in c(1:2,6:16,18:19)){ # No Myel_a1_2_3, Neutrophil
  Idents(merge.integrated_4) <- "clusters.group"
  
  marker <- FindMarkers(merge.integrated_4, ident.1 = paste0(paste(clusters$Var1[i]),"_A_Ctr"), ident.2 = paste0(paste(clusters$Var1[i]),"_B_CIA"), 
                        only.pos = F, min.pct = 0.1, logfc.threshold = 0.25, verbose = FALSE)
  marker <- subset(marker, marker$p_val_adj < 0.05)
  marker$gene <- rownames(marker)
  marker <- dplyr ::arrange(marker, desc(marker$avg_log2FC))
  write.table(marker,paste0("./Output_new4_30/Try13_merge_Ctr_CIA_JAK1i/marker_",paste(clusters$Var1[i]),"_CtrvsCIA.txt"), sep="\t", quote=FALSE, row.names=FALSE)  
  
  Idents(merge.integrated_4) <- "clusters"
  merge.integrated_4_sub <- subset(merge.integrated_4, idents = paste(clusters$Var1[i]))
  Idents(merge.integrated_4_sub) <- "group"
  avg <- data.frame(log1p(AverageExpression(merge.integrated_4_sub, verbose = FALSE)$RNA))
  avg$gene <- rownames(avg)
  
  plotlist[[i]] <- ggplot(avg, aes(A_Ctr, B_CIA)) + geom_point(size=0.1) +
    geom_point(data=avg[marker$gene, ], aes(A_Ctr, B_CIA), size=0.1, colour="red") +
    ggtitle(paste(clusters$Var1[i]))
  
}
print(plotlist)
dev.off()

# Each cluster Ctr_vs_JAK1i
plotlist <- list()
pdf("./Output_new4_30/Try13_merge_Ctr_CIA_JAK1i/marker_eachcluster_CtrvsJAK1i.pdf", useDingbats = F, height = 5, width = 5)
for (i in c(1:2,6:16,18:19)){ # No Myel_a1_2_3, Neutrophil
  Idents(merge.integrated_4) <- "clusters.group"
  
  marker <- FindMarkers(merge.integrated_4, ident.1 = paste0(paste(clusters$Var1[i]),"_A_Ctr"), ident.2 = paste0(paste(clusters$Var1[i]),"_C_JAK1i"), 
                        only.pos = F, min.pct = 0.1, logfc.threshold = 0.25, verbose = FALSE)
  marker <- subset(marker, marker$p_val_adj < 0.05)
  marker$gene <- rownames(marker)
  marker <- dplyr ::arrange(marker, desc(marker$avg_log2FC))
  write.table(marker,paste0("./Output_new4_30/Try13_merge_Ctr_CIA_JAK1i/marker_",paste(clusters$Var1[i]),"_CtrvsJAK1i.txt"), sep="\t", quote=FALSE, row.names=FALSE)  
  
  Idents(merge.integrated_4) <- "clusters"
  merge.integrated_4_sub <- subset(merge.integrated_4, idents = paste(clusters$Var1[i]))
  Idents(merge.integrated_4_sub) <- "group"
  avg <- data.frame(log1p(AverageExpression(merge.integrated_4_sub, verbose = FALSE)$RNA))
  avg$gene <- rownames(avg)
  
  plotlist[[i]] <- ggplot(avg, aes(A_Ctr, C_JAK1i)) + geom_point(size=0.1) +
    geom_point(data=avg[marker$gene, ], aes(A_Ctr, C_JAK1i), size=0.1, colour="red") +
    ggtitle(paste(clusters$Var1[i]))
  
}
print(plotlist)
dev.off()

# Each cluster CIA vs JAK1i
plotlist <- list()
pdf("./Output_new4_30/Try13_merge_Ctr_CIA_JAK1i/marker_eachcluster_CIAvsJAK1i.pdf", useDingbats = F, height = 5, width = 5)
for (i in 1:length(levels(merge.integrated_4$clusters))){
  Idents(merge.integrated_4) <- "clusters.group"
  
  marker <- FindMarkers(merge.integrated_4, ident.1 = paste0(paste(clusters$Var1[i]),"_B_CIA"), ident.2 = paste0(paste(clusters$Var1[i]),"_C_JAK1i"), 
                        only.pos = F, min.pct = 0.1, logfc.threshold = 0.25, verbose = FALSE)
  marker <- subset(marker, marker$p_val_adj < 0.05)
  marker$gene <- rownames(marker)
  marker <- dplyr ::arrange(marker, desc(marker$avg_log2FC))
  write.table(marker,paste0("./Output_new4_30/Try13_merge_Ctr_CIA_JAK1i/marker_",paste(clusters$Var1[i]),"_CIAvsJAK1i.txt"), sep="\t", quote=FALSE, row.names=FALSE)  
  
  Idents(merge.integrated_4) <- "clusters"
  merge.integrated_4_sub <- subset(merge.integrated_4, idents = paste(clusters$Var1[i]))
  Idents(merge.integrated_4_sub) <- "group"
  avg <- data.frame(log1p(AverageExpression(merge.integrated_4_sub, verbose = FALSE)$RNA))
  avg$gene <- rownames(avg)
  
  plotlist[[i]] <- ggplot(avg, aes(B_CIA, C_JAK1i)) + geom_point(size=0.1) +
    geom_point(data=avg[marker$gene, ], aes(B_CIA, C_JAK1i), size=0.1, colour="red") +
    ggtitle(paste(clusters$Var1[i]))
  }
print(plotlist)
dev.off()

