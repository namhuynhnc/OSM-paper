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
library(readxl)
library(ggrepel)
set.seed(1234567)

options(ggrepel.max.overlaps = Inf)

# Load data
merge.integrated_3 <- readRDS("./Save_new/seuset_merge_Ctr_CIA_JAK1i_30_2.RDS")

# Select Fibroblast
Fibro <- subset(merge.integrated_3, idents = c("Fibroblast"))
Fibro

# DefaultAssay
DefaultAssay(object = Fibro) <- "RNA"
merge.list <- SplitObject(object = Fibro, split.by = "group")

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
merge.integrated <- FindNeighbors(object = merge.integrated, dims = 1:30)

# Save
saveRDS(merge.integrated, "./Save_new/Fibro_Ctr_CIA_JAK1i_30_1.RDS")
merge.integrated <- readRDS("./Save_new/Fibro_Ctr_CIA_JAK1i_30_1.RDS")

# Clustering
Fibro <- FindClusters(merge.integrated, verbose = FALSE, resolution = 0.3)
table(Fibro@active.ident)

Fibro <- RenameIdents(object=Fibro,"3"="Chondrocyte")
Fibro <- RenameIdents(object=Fibro,"2"="Fibro_4")
Fibro <- RenameIdents(object=Fibro,"1"="Fibro_3")
Fibro <- RenameIdents(object=Fibro,"0"="Fibro_2")
Fibro <- RenameIdents(object=Fibro,"4"="Fibro_1")

Fibro$saved.idents <- Idents(object = Fibro)

# DefaultAssay
DefaultAssay(object = Fibro) <- "RNA"
Fibro <- NormalizeData(Fibro, verbose = FALSE)

# Save
saveRDS(Fibro, "./Save_new/Fibro_Ctr_CIA_JAK1i_30_2.RDS")
Fibro <- readRDS("./Save_new/Fibro_Ctr_CIA_JAK1i_30_2.RDS")

# Plot
pdf("./Output_new4_30/Try14_Fibro_Ctr_CIA_JAK1i/dimplot_0.3.pdf", useDingbats = F, height = 4, width = 7)
DimPlot(Fibro, label = TRUE)
DimPlot(Fibro, label = TRUE, label.size = 0, group.by = "group")
dev.off()

#
markers.to.plot <- c('Pdpn','Thy1','Fap','Cd55',
                     'Cd34','Il6','Tnfsf11','Tnfrsf11b',
                     "Hapln1","Prg4","Ets1","Osm",
                     "Cdh11","Vcam1","Cd276","Sdc1")

pdf("./Output_new4_30/Try14_Fibro_Ctr_CIA_JAK1i/featureplot_1.pdf", useDingbats = F, height = 10, width = 16)
FeaturePlot(Fibro, markers.to.plot, pt.size = 0.1,
            cols = c("lightgrey", "#005493")
            , min.cutoff = "q10", max.cutoff = "q90"
)
dev.off()

#
markers.to.plot <- c('Pdpn','Fap','Thy1','Cd34',
                     "Vcam1","Cdh11",
                     "Dkk3",
                     "Cd276","Sdc1",
                     'Cd55',"Prg4","Clic5","Tspan15",
                     'Tnfrsf11b',"Hapln1")

pdf("./Output_new4_30/Try14_Fibro_Ctr_CIA_JAK1i/DotPlot_0.3.pdf", useDingbats = F, height = 4, width = 7.5)
DotPlot(Fibro, features = markers.to.plot, dot.scale = 8) + RotatedAxis()
dev.off()

#
markers.to.plot1 <- c("Il17a","Il1b","Tnf","Il6","Osm","Lif","Tnfsf11",'Tnfrsf11a',"Mmp3","Mmp13","Mmp14",
                      "Il2ra","Il2rb","Il2rg","Il4ra","Il7r","Il15ra","Il21r","Il6ra","Il6st","Osmr","Lifr",
                      "Ifnar1","Ifnar2","Il10ra","Il10rb","Il22ra2","Ifngr1","Ifngr2","Csf3r","Epor","Ghr",
                      "Jak1","Jak2","Jak3","Tyk2","Stat1","Stat2","Stat3","Stat5a")

pdf("./Output_new4_30/Try14_Fibro_Ctr_CIA_JAK1i/DotPlot.pdf", useDingbats = F, height = 5, width = 14)
DotPlot(Fibro, idents = c("Fibro_1", "Fibro_2", "Fibro_3", "Fibro_4"), features = markers.to.plot1, dot.scale = 8, cols = c("blue","blue","blue"),
        split.by = "group") + RotatedAxis()
dev.off()

#
markers.to.plot2 <- c("Il17a","Il1b","Tnf","Il6","Osm","Tnfsf11","Mmp3","Mmp13")
pdf("./Output_new4_30/Try14_Fibro_Ctr_CIA_JAK1i/DotPlot2.pdf", useDingbats = F, height = 5, width = 6)
DotPlot(Fibro, idents = c("Fibro_1", "Fibro_2", "Fibro_3", "Fibro_4"), features = markers.to.plot2, dot.scale = 8, cols = c("blue","blue","blue"),
        split.by = "group") + RotatedAxis()
dev.off()

#
pdf("./Output_new4_30/Try14_Fibro_Ctr_CIA_JAK1i/VlnPlot_1.pdf", useDingbats = F, height = 25, width = 10)
markers.to.plot <- c("Thy1" ,"Vcam1","Cdh11", 'Il6',"Osm",'Tnfsf11','Tnfrsf11b',"Ets1", "Ets2", "Cd276", "Ccr1")
plots <- VlnPlot(Fibro, features = markers.to.plot
                 , split.by = "group" , cols = c("#078992","#F8766D","#D39200"),
                 pt.size = 0.1, combine = FALSE, split.plot = F)
CombinePlots(plots = plots, ncol = 1)

markers.to.plot <- c("Il2ra","Il2rb","Il2rg","Il4ra","Il7r","Il15ra","Il21r","Il6ra","Il6st","Osmr","Lifr")
plots <- VlnPlot(Fibro, features = markers.to.plot
                 , split.by = "group" , cols = c("#078992","#F8766D","#D39200"),
                 pt.size = 0.1, combine = FALSE, split.plot = F)
CombinePlots(plots = plots, ncol = 1)
dev.off()

#
pdf("./Output_new4_30/Try14_Fibro_Ctr_CIA_JAK1i/DotPlot_0.3g.pdf", useDingbats = F, height = 3, width = 6)
VlnPlot(Fibro, features = "Ets1", idents = c("Fibro_1", "Fibro_2", "Fibro_3", "Fibro_4")
                 , split.by = "group" , cols = c("#078992","#F8766D","#D39200"),
                 pt.size = 0.1, combine = FALSE, split.plot = F)
dev.off()

#
Fibro$clusters <- Idents(Fibro)
Idents(Fibro) <- "group"
CIA <- subset(Fibro, idents= c("A_Ctr","B_CIA"))
Idents(Fibro) <- "clusters"
Idents(CIA) <- "clusters"

pdf("./Output_new4_30/Try14_Fibro_Ctr_CIA_JAK1i/VlnPlot_2.pdf", useDingbats = F, height = 25, width = 7)
markers.to.plot <- c("Thy1" ,"Vcam1","Cdh11", 'Il6',"Osm",'Tnfsf11','Tnfrsf11b',"Ets1", "Ets2", "Cd276", "Ccr1")
plots <- VlnPlot(Fibro, features = markers.to.plot
                 , split.by = "group" , cols = c("#078992","#F8766D"),
                 pt.size = 0.1, combine = FALSE, split.plot = F)
CombinePlots(plots = plots, ncol = 1)

markers.to.plot <- c("Il2ra","Il2rb","Il2rg","Il4ra","Il7r","Il15ra","Il21r","Il6ra","Il6st","Osmr","Lifr")
plots <- VlnPlot(Fibro, features = markers.to.plot
                 , split.by = "group" , cols = c("#078992","#F8766D"),
                 pt.size = 0.1, combine = FALSE, split.plot = F)
CombinePlots(plots = plots, ncol = 1)

markers.to.plot <- c("Csf3r","Ifnar1","Ifnar2","Il10ra","Il10rb","Il22ra2","Ifngr1","Ifngr2","Stat1","Stat3","Stat5a")
plots <- VlnPlot(Fibro, features = markers.to.plot
                 , split.by = "group" , cols = c("#078992","#F8766D"),
                 pt.size = 0.1, combine = FALSE, split.plot = F)
CombinePlots(plots = plots, ncol = 1)
dev.off()

pdf("./Output_new4_30/Try14_Fibro_Ctr_CIA_JAK1i/VlnPlot_3.pdf", useDingbats = F, height = 3, width = 3.5)
VlnPlot(CIA, features = c("Il6"), pt.size = 0.1, cols = c("#078992","#F8766D"), log = T, split.by = "group")
VlnPlot(CIA, features = c("Tnfsf11"), pt.size = 0.1, cols = c("#078992","#F8766D"), log = T, split.by = "group")
dev.off()

# Markers
markers <- FindAllMarkers(object = Fibro, only.pos = F, min.pct = 0.1, logfc.threshold = 0.1)
markers <- subset(markers, markers$p_val_adj < 0.05)
markers <- dplyr ::arrange(markers, desc(markers$avg_log2FC))
markers <- dplyr ::arrange(markers, markers$cluster)
write.table(markers,"./Output_new4_30/Try14_Fibro_Ctr_CIA_JAK1i/markers_0.3.txt", sep="\t", quote=FALSE, row.names=FALSE)

# Markers CIA
Idents(Fibro) <- "group"
Fibro_CIA <- subset(Fibro, idents = "B_CIA")
Idents(Fibro_CIA) <- "saved.idents"
markers <- FindAllMarkers(object = Fibro_CIA, only.pos = T, min.pct = 0.1, logfc.threshold = 0.1)
markers <- subset(markers, markers$p_val_adj < 0.05)
markers <- dplyr ::arrange(markers, desc(markers$avg_log2FC))
markers <- dplyr ::arrange(markers, markers$cluster)
write.table(markers,"./Output_new4_30/Try14_Fibro_Ctr_CIA_JAK1i/markers_0.3_CIA.txt", sep="\t", quote=FALSE, row.names=FALSE)

# Check IL6 target genes in CIA fibroblast up-regulated genes
avg.markers <- read.table("./Output_new4_30/Try14_Fibro_Ctr_CIA_JAK1i/avg.markers_CIA_JAK1i.txt", sep="\t", header = T)
avg.markers_up <- subset(avg.markers, avg.markers$avg_log2FC>0)

#IL6 (mouse)
list <- read.table("/Volumes/Nam_HDD/Ref data (sc_bulk RNAseq)/phs001529.v1.p1/Results/JAK1-cytokine/Ref/IL6/result.txt",sep="\t", header = T)
list <- subset(list, list$adj.P.Val<0.05)
list <- subset(list, list$logFC>0)
list[list==""]<-NA
listn_il6 <- na.omit(list)
listn_il6 <- listn_il6[!duplicated(listn_il6$Gene.symbol),] # remove duplicated genes

list.merge_il6 <- merge(x=avg.markers_up, y=listn_il6, by.x="gene", by.y="Gene.symbol")
list.merge_il6 <- dplyr :: arrange(list.merge_il6, desc(list.merge_il6$avg_log2FC))
write.table(list.merge_il6,"./Output_new4_30/Try14_Fibro_Ctr_CIA_JAK1i/avg.Fibro_CIAvsJAK1i_heatmap_IL6_target.txt", sep="\t", quote=FALSE, row.names=F)

Idents(Fibro) <- "group"
cluster.averages <- AverageExpression(Fibro, return.seurat = T)

pdf("./Output_new4_30/Try14_Fibro_Ctr_CIA_JAK1i/avg.Fibro_CIAvsJAK1i_heatmap_IL6_target.pdf", useDingbats = F, height = 15, width = 4)
DoHeatmap(subset(cluster.averages, idents = c("B_CIA","C_JAK1i")), features = list.merge_il6$gene, size = 1.7, draw.lines = F,raster =F,
          group.colors = c("#F8766D","#D39200")) + scale_fill_gradientn(colors = c("blue","white","red"))

dev.off()

#OSM (mouse)
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

list.merge_osm <- merge(x=avg.markers_up, y=listn_osm, by.x="gene", by.y="Gene.symbol")
list.merge_osm <- dplyr :: arrange(list.merge_osm, desc(list.merge_osm$avg_log2FC))
write.table(list.merge_osm,"./Output_new4_30/Try14_Fibro_Ctr_CIA_JAK1i/avg.Fibro_CIAvsJAK1i_heatmap_OSM_target.txt", sep="\t", quote=FALSE, row.names=F)

pdf("./Output_new4_30/Try14_Fibro_Ctr_CIA_JAK1i/avg.Fibro_CIAvsJAK1i_heatmap_OSM_target.pdf", useDingbats = F, height = 15, width = 4)
DoHeatmap(subset(cluster.averages, idents = c("B_CIA","C_JAK1i")), features = list.merge_osm$gene, size = 1.6, draw.lines = F,raster =F,
          group.colors = c("#F8766D","#D39200")) + scale_fill_gradientn(colors = c("blue","white","red"))
dev.off()

# Il6 vs OsM
list.il6_osm <- merge(x=listn_il6 , y=listn_osm, by.x="Gene.symbol", by.y="Gene.symbol")
list.merge_il6_osm <- merge(x=avg.markers_up , y=list.il6_osm, by.x="gene", by.y="Gene.symbol")
pdf("./Output_new4_30/Try14_Fibro_Ctr_CIA_JAK1i/avg.Fibro_CIAvsJAK1i_heatmap_IL6_OSM_target.pdf", useDingbats = F, height = 8, width = 4)
DoHeatmap(subset(cluster.averages, idents = c("B_CIA","C_JAK1i")), features = list.merge_il6_osm$gene, size = 2, draw.lines = F,raster =F,
          group.colors = c("#F8766D","#D39200")) + scale_fill_gradientn(colors = c("blue","white","red"))

dev.off()

