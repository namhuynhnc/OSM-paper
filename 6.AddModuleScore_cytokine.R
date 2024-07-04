# Set working directory
setwd("/Volumes/Nam_HDD/Ref data (sc_bulk RNAseq)/CIA")

library(scater)
library(scran)
library(SingleCellExperiment)
library(Biobase)
library(irlba)
library(ggplot2)
library(gplots)
library(tidyverse)
library(jsonlite)
library(XVector)
library(Seurat)
library(pheatmap)
library(mclust)
library(edgeR)
library(dplyr)
library(sctransform)
library(RColorBrewer)
library(pheatmap)
library(readxl)
library(gdata)
library(dplyr)
library(taRifx)
set.seed(1234567)

# All cell types
seuset <- readRDS("./Save_new/seuset_merge_Ctr_CIA_JAK1i_30_3.RDS")
DimPlot(seuset, label = TRUE, label.size = 2) + NoLegend()

seuset$clusters.group <- paste(Idents(seuset), seuset$group, sep = "_")
seuset$clusters <- Idents(seuset)
clusters <- as.data.frame(table(seuset$clusters))

#DotPlot
DotPlot(seuset, features = c("Ifng","Il22","Il10","Ifnb1","Ifna2","Csf2","Il6","Il21","Il15","Il9","Il7","Il4","Il2")) + RotatedAxis()
DotPlot(seuset, features = c("Tyk2","Jak3","Jak2","Jak1")) + RotatedAxis()
FeaturePlot(
  seuset, 
  c("Jak1"),
  cols = c("lightgrey", "#005493")
  # , min.cutoff = "q10", max.cutoff = "q90"
)

#Match mouse and human genes
load("/Volumes/Nam_HDD/Ref data (sc_bulk RNAseq)/Ref/Single cell RNAseq Immunity 2019/Data/7098575/bm.query_mouse_allevo.RData")
genelist <- data.frame(bm.query_mouse_allevo_list[[1]])
genelist2 <- data.frame(bm.query_mouse_allevo_list[[5]]) # match human gene
genelist3 <- merge(genelist, genelist2, by.x=1, by.y=1)
genelist3[genelist3==""]<-NA
genelist3 <- na.omit(genelist3)
genelist3 <-  genelist3[,c(1:4)]

#Match gene accession ID
accession <- read_excel("/Volumes/Nam_HDD/Ref data (sc_bulk RNAseq)/phs001529.v1.p1/Results/JAK1-cytokine/Ref/IL10/gene_accession.xlsx",sheet=1, col_names = T)

##IL2 (mouse)
list <- read_excel("/Volumes/Nam_HDD/Ref data (sc_bulk RNAseq)/phs001529.v1.p1/Results/JAK1-cytokine/Ref/IL15_21/40425_2014_28_MOESM1_ESM.xlsx",sheet=2, col_names = T)
list <- list[,c(9,10)]
list <- subset(list, list[[2]]>1)
list <- dplyr ::arrange(list, desc(list$IL2FC))

#IL2 target list
list <- data.frame(list[1:50,1])
colnames(list) <- "gene"
Il2 <- list$gene 

##IL4 (mouse)
list <- read_excel("/Volumes/Nam_HDD/Ref data (sc_bulk RNAseq)/phs001529.v1.p1/Results/JAK1-cytokine/Ref/IL15_21/40425_2014_28_MOESM1_ESM.xlsx",sheet=2, col_names = T)
list <- list[,c(9,11)]
list <- subset(list, list[[2]]>1)
list <- dplyr ::arrange(list, desc(list$IL4FC))

#IL4 target list
list <- data.frame(list[1:50,1])
colnames(list) <- "gene"
Il4 <- list$gene 

##IL7 (mouse)
list <- read_excel("/Volumes/Nam_HDD/Ref data (sc_bulk RNAseq)/phs001529.v1.p1/Results/JAK1-cytokine/Ref/IL7/IL7_IFNg.xls",sheet=1, col_names = T)
list <- subset(list, list$log2FC>0)
list <- dplyr ::arrange(list, desc(list$log2FC))

#IL7 target list
list <- data.frame(list[1:50,1])
colnames(list) <- "gene"
Il7 <- list$gene 

##IFNg (human)
list <- read_excel("/Volumes/Nam_HDD/Ref data (sc_bulk RNAseq)/phs001529.v1.p1/Results/JAK1-cytokine/Ref/IFN/interferome.org.2171-2169-2170.xlsx",sheet=1, col_names = T)
list <- list[,c(2,5)]
list <- data.frame(list[!duplicated(list[[2]]),])
colnames(list)[[2]] <- c("gene")

#Match mouse gene
list <- merge(list, genelist3, by.x="gene", by.y=4)
list <- dplyr ::arrange(list, desc(list$Fold.Change))
list <- list[!duplicated(list$external_gene_name),]

#IFNg target list
list <- data.frame(list[1:50,4])
colnames(list) <- "gene"
Ifng <- list$gene 

#IFNa (human)
list <- read_excel("/Volumes/Nam_HDD/Ref data (sc_bulk RNAseq)/phs001529.v1.p1/Results/JAK1-cytokine/Ref/IFN/interferome.org.2171-2169-2170.xlsx",sheet=2, col_names = T)
list <- list[,c(2,5)]
list <- data.frame(list[!duplicated(list[[2]]),])
colnames(list)[[2]] <- c("gene")

#Match mouse gene
list <- merge(list, genelist3, by.x="gene", by.y=4)
list <- dplyr ::arrange(list, desc(list$Fold.Change))
list <- list[!duplicated(list$external_gene_name),]

#IFNa target list
list <- data.frame(list[1:50,4])
colnames(list) <- "gene"
Ifna <- list$gene

#IFNb (human)
list <- read_excel("/Volumes/Nam_HDD/Ref data (sc_bulk RNAseq)/phs001529.v1.p1/Results/JAK1-cytokine/Ref/IFN/interferome.org.2171-2169-2170.xlsx",sheet=3, col_names = T)
list <- list[,c(2,5)]
list <- data.frame(list[!duplicated(list[[2]]),])
colnames(list)[[2]] <- c("gene")

#Match mouse gene
list <- merge(list, genelist3, by.x="gene", by.y=4)
list <- dplyr ::arrange(list, desc(list$Fold.Change))
list <- list[!duplicated(list$external_gene_name),]

#IFNb target list
list <- data.frame(list[1:50,4])
colnames(list) <- "gene"
Ifnb <- list$gene

#IL10 (mouse)
list <- read_excel("/Volumes/Nam_HDD/Ref data (sc_bulk RNAseq)/phs001529.v1.p1/Results/JAK1-cytokine/Ref/IL10/IL10.xlsx",sheet=1, col_names = T)
#Match accession ID
list <- merge(list, accession, by.x=1, by.y=2)
list <- na.omit(list)
list <- subset(list, list$FC>1)
list <- dplyr ::arrange(list, desc(list$FC))
#list[32,3] <- "Mtlrp"
list <- list[!duplicated(list$Gene),]

#IL10 target list
list <- data.frame(list[1:50,3])
colnames(list) <- "gene"
Il10 <- list$gene

#OSM (mouse)
list1 <- read_excel("/Volumes/Nam_HDD/Ref data (sc_bulk RNAseq)/phs001529.v1.p1/Results/JAK1-cytokine/Ref/Osm/Table1.xlsx",sheet=1, col_names = T)
list1 <- list1[-1,]

list2 <- read_excel("/Volumes/Nam_HDD/Ref data (sc_bulk RNAseq)/phs001529.v1.p1/Results/JAK1-cytokine/Ref/Osm/jbc.M116.748483-1.xlsx",sheet=1, col_names = F)
colnames(list2) <- list2[2,]
list2 <- list2[-c(1:2),]

list3 <- read_excel("/Volumes/Nam_HDD/Ref data (sc_bulk RNAseq)/phs001529.v1.p1/Results/JAK1-cytokine/Ref/Osm/jbc.M116.748483-1.xlsx",sheet=2, col_names = F)
colnames(list3) <- list3[2,]
list3 <- list3[-c(1:2),]
list3 <- subset(list3, list3$logFC>0)

#OSM target list
list <- data.frame(c(list1$Rank, list2$Symbol, list3$Symbol))
colnames(list) <- "Gene.symbol"
list[list=="NA"]<-NA
list <- na.omit(list)
list <- data.frame(list[!duplicated(list$Gene.symbol),]) # remove duplicated genes
colnames(list) <- "gene"
Osm <- list$gene[1:50]

#IL6 (mouse)
list <- read.table("/Volumes/Nam_HDD/Ref data (sc_bulk RNAseq)/phs001529.v1.p1/Results/JAK1-cytokine/Ref/IL6/result.txt",sep="\t", header = T)
list <- subset(list, list$adj.P.Val<0.05)
list[list==""]<-NA
list <- na.omit(list)
list <- dplyr ::arrange(list, desc(list$logFC))
#list[61, 4] <- "Gm4070"

list <- list[!duplicated(list$Gene.symbol),]

#IL6 target list
list <- data.frame(list[1:50,4])
colnames(list) <- "gene"
Il6 <- list$gene

#IL22 (mouse)
list <- read_excel("/Volumes/Nam_HDD/Ref data (sc_bulk RNAseq)/phs001529.v1.p1/Results/JAK1-cytokine/Ref/IL22/journal.pbio.3000540.s006.xlsx",sheet=1, col_names = T)
colnames(list) <- list[1,]
list <- list[-1,c(1,7,9)]
list$logFC <- as.numeric(list$logFC)
list$adj.P.Val <- as.numeric(list$adj.P.Val)
list <- subset(list, list$adj.P.Val <0.05 & list$logFC>0)

list <- dplyr ::arrange(list, desc(list$logFC))
list <- list[!duplicated(list$ID),]

#IL22 target list
list <- data.frame(list[1:50,1])
colnames(list) <- "gene"
Il22 <- list$gene

#IL15 (mouse)
list <- read_excel("/Volumes/Nam_HDD/Ref data (sc_bulk RNAseq)/phs001529.v1.p1/Results/JAK1-cytokine/Ref/IL15_21/40425_2014_28_MOESM1_ESM.xlsx",sheet=2, col_names = T)
list <- list[,c(9,13)]
list <- subset(list, list[[2]]>1)
list <- dplyr ::arrange(list, desc(list$IL15FC))

#IL15 target list
list <- data.frame(list[1:50,1])
colnames(list) <- "gene"
Il15 <- list$gene 

#IL21 (mouse)
list <- read_excel("/Volumes/Nam_HDD/Ref data (sc_bulk RNAseq)/phs001529.v1.p1/Results/JAK1-cytokine/Ref/IL15_21/40425_2014_28_MOESM1_ESM.xlsx",sheet=2, col_names = T)
list <- list[,c(9,14)]
list <- subset(list, list[[2]]>1)
list <- dplyr ::arrange(list, desc(list$IL21FC))

#IL21 target list
list <- data.frame(list[1:50,1])
colnames(list) <- "gene"
Il21 <- list$gene 

#GM-CSF (human mouse)
list <- read_excel("/Volumes/Nam_HDD/Ref data (sc_bulk RNAseq)/phs001529.v1.p1/Results/JAK1-cytokine/Ref/GM-CSF/GMCSF.xlsx",sheet=1, col_names = T)

#Match mouse gene
listb <- merge(list[1:100,], genelist3, by.x="gene", by.y=4)
listc <- listb[,c(4,2)]
colnames(listc)[1] <- "gene"
listc <- rbind.data.frame(listc, list[101:200,])
listc <- subset(listc, listc$FC>0)
list <- listc[!duplicated(listc$gene),] # remove duplicated genes
list <- dplyr ::arrange(list, desc(list$FC))

#GM-CSF target list
list <- data.frame(list[1:50,1])
colnames(list) <- "gene"
Csf2 <- list$gene 

# Epo (human)
list <- read.table("/Volumes/Nam_HDD/Ref data (sc_bulk RNAseq)/phs001529.v1.p1/Results/JAK1-cytokine/Ref/EPO/GSE6260.top.table.tsv", sep="\t", header = T)
list <- list[!(list$Gene.symbol == ""), ]
list <- list[!duplicated(list$Gene.symbol),] # remove duplicated genes
list <- subset(list, list$P.Val<0.05 & list$logFC<0)

#Match mouse gene
list <- merge(list, genelist3, by.x="Gene.symbol", by.y=4)
list <- dplyr ::arrange(list, list$logFC) #down in control is up in Epo
list <- list[!duplicated(list$external_gene_name),]

#EPO target list
list <- data.frame(list[1:50,10])
colnames(list) <- "gene"
Epo <- list$gene 

# Gh (mouse)
list <- read.table("/Volumes/Nam_HDD/Ref data (sc_bulk RNAseq)/phs001529.v1.p1/Results/JAK1-cytokine/Ref/GH/GSE2120.top.table.tsv", sep="\t", header = T)
list <- subset(list, list$P.Val<0.05 & list$logFC<0) #down in control
list <- list[!(list$Gene.symbol == ""), ]
list <- list[!duplicated(list$Gene.symbol),] # remove duplicated genes

#GH target list
list <- data.frame(list[1:50,7])
colnames(list) <- "gene"
Gh <- list$gene 

# target genes
target_data <- cbind(Il2, Il4, Il7, Il15, Il21, Il6, Osm, Ifna, Ifnb, Il10, Il22, Ifng, Csf2, Epo, Gh)
write.table(target_data,"./Output_new4_30/Try15_AddmoduleScore_target/cytokine_targets.txt", sep="\t", quote=FALSE, row.names=FALSE)

# AddModuleScore
DefaultAssay(object = seuset) <- "SCT"
seuset_AddModuleScore <- AddModuleScore(
  object = seuset,
  features = list(Il2, Il4, Il7, Il15, Il21, Il6, Osm, Ifna, Ifnb, Il10, Il22, Ifng, Csf2, Epo, Gh)
  , ctrl = 100, nbin = 24,
  name = 'target')

# Heatmap all
target <- c("target1","target2","target3","target4",
            "target5","target6","target7","target8",
            "target9","target10","target11","target12",
            "target13","target14","target15")

seuset_AddModuleScore[['module']] <- CreateAssayObject(data = t(x = FetchData(object = seuset_AddModuleScore, vars = target)))
DoHeatmap(object = seuset_AddModuleScore, features = target, assay = 'module', slot = 'data')

# Heatmap group
seuset_AddModuleScore$clusters <- Idents(seuset_AddModuleScore)
Idents(seuset_AddModuleScore) <- "group"
seuset_CIA <- subset(seuset_AddModuleScore, idents= "B_CIA")
Idents(seuset_CIA) <- "clusters"

seuset_JAK1i <- subset(seuset_AddModuleScore, idents= "C_JAK1i")
Idents(seuset_JAK1i) <- "clusters"

## Fibro_2
Idents(seuset_AddModuleScore) <- "clusters"
fibro_2 <- subset(seuset_AddModuleScore, idents= "Fibro_2")
Idents(fibro_2) <- "group"
fibro_2 <- AverageExpression(fibro_2, return.seurat = T)

#
pdf("./Output_new4_30/Try15_AddmoduleScore_target/Heatmap_fibro2.pdf", useDingbats = F, height = 3, width = 3)
DoHeatmap(subset(fibro_2, idents = c("B_CIA","C_JAK1i")), assay = 'module', slot = 'data', raster=F,
          features = target, size = 2, draw.lines = F, group.colors = c("#F8766D","#D39200")) + scale_fill_gradientn(colors = c("blue","white", "red"))
dev.off()

#
Idents(seuset_AddModuleScore) <- "group"
CIA_JAK1i <- subset(seuset_AddModuleScore, idents = c("B_CIA","C_JAK1i"))
Idents(CIA_JAK1i) <- "clusters"
pdf("./Output_new4_30/Try15_AddmoduleScore_target/VlnPlot1.pdf", useDingbats = F, height = 3, width = 6)
VlnPlot(CIA_JAK1i, features = c("target1"), split.by = "group", log = F, split.plot=T, cols = c("#F8766D","#D39200"))
VlnPlot(CIA_JAK1i, features = c("target2"), split.by = "group", log = F, split.plot=T, cols = c("#F8766D","#D39200"))
VlnPlot(CIA_JAK1i, features = c("target3"), split.by = "group", log = F, split.plot=T, cols = c("#F8766D","#D39200"))
VlnPlot(CIA_JAK1i, features = c("target4"), split.by = "group", log = F, split.plot=T, cols = c("#F8766D","#D39200"))
VlnPlot(CIA_JAK1i, features = c("target5"), split.by = "group", log = F, split.plot=T, cols = c("#F8766D","#D39200"))
VlnPlot(CIA_JAK1i, features = c("target6"), split.by = "group", log = F, split.plot=T, cols = c("#F8766D","#D39200"))
VlnPlot(CIA_JAK1i, features = c("target7"), split.by = "group", log = F, split.plot=T, cols = c("#F8766D","#D39200"))
VlnPlot(CIA_JAK1i, features = c("target8"), split.by = "group", log = F, split.plot=T, cols = c("#F8766D","#D39200"))
VlnPlot(CIA_JAK1i, features = c("target9"), split.by = "group", log = F, split.plot=T, cols = c("#F8766D","#D39200"))
VlnPlot(CIA_JAK1i, features = c("target10"), split.by = "group", log = F, split.plot=T, cols = c("#F8766D","#D39200"))
VlnPlot(CIA_JAK1i, features = c("target11"), split.by = "group", log = F, split.plot=T, cols = c("#F8766D","#D39200"))
VlnPlot(CIA_JAK1i, features = c("target12"), split.by = "group", log = F, split.plot=T, cols = c("#F8766D","#D39200"))
VlnPlot(CIA_JAK1i, features = c("target13"), split.by = "group", log = F, split.plot=T, cols = c("#F8766D","#D39200"))
VlnPlot(CIA_JAK1i, features = c("target14"), split.by = "group", log = F, split.plot=T, cols = c("#F8766D","#D39200"))
VlnPlot(CIA_JAK1i, features = c("target15"), split.by = "group", log = F, split.plot=T, cols = c("#F8766D","#D39200"))
dev.off()

#
pdf("./Output_new4_30/Try15_AddmoduleScore_target/Heatmap1.pdf", useDingbats = F, height = 5, width = 12)
DoHeatmap(object = seuset_CIA, features = target, assay = 'module', slot = 'data', size = 3, lines.width= 1,raster=F, group.colors = c("#F8766D","#D39200")) + scale_fill_gradientn(colors = c("blue","white", "red"))
DoHeatmap(object = seuset_JAK1i, features = target, assay = 'module', slot = 'data', size = 3, lines.width= 1,raster=F, group.colors = c("#F8766D","#D39200")) + scale_fill_gradientn(colors = c("blue","white", "red"))
dev.off()

#
pdf("./Output_new4_30/Try15_AddmoduleScore_target/Heatmap2.pdf", useDingbats = F, height = 3, width = 6)
Idents(seuset_AddModuleScore) <- "clusters"
seuset_avg <- AverageExpression(seuset_AddModuleScore, return.seurat = T)
DoHeatmap(object = seuset_avg, features = target, assay = 'module', slot = 'data', draw.lines = F, raster=F, size = 2, group.colors = c("#F8766D","#D39200")) + scale_fill_gradientn(colors = c("blue","white", "red"))
seuset_CIA <- AverageExpression(seuset_CIA, return.seurat = T)
DoHeatmap(object = seuset_CIA, features = target, assay = 'module', slot = 'data', draw.lines = F, raster=F, size = 2, group.colors = c("#F8766D","#D39200")) + scale_fill_gradientn(colors = c("blue","white", "red"))
seuset_JAK1i <- AverageExpression(seuset_JAK1i, return.seurat = T)
DoHeatmap(object = seuset_JAK1i, features = target, assay = 'module', slot = 'data', draw.lines = F, raster=F,size = 2, group.colors = c("#F8766D","#D39200")) + scale_fill_gradientn(colors = c("blue","white", "red"))
dev.off()

pdf("./Output_new4_30/Try15_AddmoduleScore_target/Heatmap3_50.pdf", useDingbats = F, height = 3, width = 12)
Idents(seuset_AddModuleScore) <- "group"
seuset_CIA_JAK1i <- subset(seuset_AddModuleScore, idents= c("B_CIA", "C_JAK1i"))
Idents(seuset_CIA_JAK1i) <- "clusters.group"
Idents(seuset_CIA_JAK1i) <- factor(x = Idents(seuset_CIA_JAK1i), levels = c("T_cell_B_CIA","B_cell_B_CIA",
                                                                            "Myel_a1_B_CIA","Myel_a2_B_CIA",
                                                                            "Myel_a3_B_CIA","Myel_b_B_CIA",
                                                                            "Myel_c1_B_CIA","Myel_c2_B_CIA",
                                                                            "Myel_c3_B_CIA","Myel_c4_B_CIA",
                                                                            "Myel_c5_B_CIA","Fibro_1_B_CIA",
                                                                            "Fibro_2_B_CIA","Fibro_3_B_CIA",
                                                                            "Fibro_4_B_CIA","Chondrocyte_B_CIA",
                                                                            "Neutrophil_B_CIA","Endothelial_B_CIA",
                                                                            "Mural_cell_B_CIA",
                                                                            "T_cell_C_JAK1i","B_cell_C_JAK1i",
                                                                            "Myel_a1_C_JAK1i","Myel_a2_C_JAK1i",
                                                                            "Myel_a3_C_JAK1i","Myel_b_C_JAK1i",
                                                                            "Myel_c1_C_JAK1i","Myel_c2_C_JAK1i",
                                                                            "Myel_c3_C_JAK1i","Myel_c4_C_JAK1i",
                                                                            "Myel_c5_C_JAK1i","Fibro_1_C_JAK1i",
                                                                            "Fibro_2_C_JAK1i","Fibro_3_C_JAK1i",
                                                                            "Fibro_4_C_JAK1i","Chondrocyte_C_JAK1i",
                                                                            "Neutrophil_C_JAK1i","Endothelial_C_JAK1i",
                                                                            "Mural_cell_C_JAK1i"))
seuset_avg <- AverageExpression(seuset_CIA_JAK1i, return.seurat = T)
DoHeatmap(object = seuset_avg, features = target, assay = 'module', slot = 'data', draw.lines = F, raster=F, size = 2, group.colors = c("#F8766D","#D39200")) + scale_fill_gradientn(colors = c("blue","white", "red"))
dev.off()

#
pdf("./Output_new4_30/Try15_AddmoduleScore_target/Dotplot.pdf", useDingbats = F, height = 6, width = 8)
DotPlot(seuset_AddModuleScore, features = target
        #, split.by = "group", cols = c("blue", "red", "green"),
        ,dot.scale = 8) + RotatedAxis()

dev.off()

pdf("./Output_new4_30/Try15_AddmoduleScore_target/FeaturePlot.pdf", useDingbats = F, height = 3, width = 12)
FeaturePlot(seuset_AddModuleScore, "target5", split.by = "group",
            cols = c("lightgrey", "#005493"), pt.size = 0.1
            , min.cutoff = "q10", max.cutoff = "q90")
FeaturePlot(seuset_AddModuleScore, "target6", split.by = "group",
            cols = c("lightgrey", "#005493"), pt.size = 0.1
            , min.cutoff = "q10", max.cutoff = "q90")
FeaturePlot(seuset_AddModuleScore, "target7", split.by = "group",
            cols = c("lightgrey", "#005493"), pt.size = 0.1
            , min.cutoff = "q30", max.cutoff = "q90"
            )
FeaturePlot(seuset_AddModuleScore, "target8", split.by = "group",
            cols = c("lightgrey", "#005493"), pt.size = 0.1
            , min.cutoff = "q10", max.cutoff = "q90")
FeaturePlot(seuset_AddModuleScore, "target9", split.by = "group",
            cols = c("lightgrey", "#005493"), pt.size = 0.1
            , min.cutoff = "q10", max.cutoff = "q90")
FeaturePlot(seuset_AddModuleScore, "target10", split.by = "group",
            cols = c("lightgrey", "#005493"), pt.size = 0.1
            , min.cutoff = "q10", max.cutoff = "q90")
FeaturePlot(seuset_AddModuleScore, "target11", split.by = "group",
            cols = c("lightgrey", "#005493"), pt.size = 0.1
            , min.cutoff = "q10", max.cutoff = "q90")
FeaturePlot(seuset_AddModuleScore, "target12", split.by = "group",
            cols = c("lightgrey", "#005493"), pt.size = 0.1
            , min.cutoff = "q10", max.cutoff = "q90")
dev.off()

#
type <- data.frame(
  seuset_AddModuleScore@meta.data$target1,
  seuset_AddModuleScore@meta.data$target2,
  seuset_AddModuleScore@meta.data$target3,
  seuset_AddModuleScore@meta.data$target4,
  seuset_AddModuleScore@meta.data$target5,
  seuset_AddModuleScore@meta.data$target6,
  seuset_AddModuleScore@meta.data$target7,
  seuset_AddModuleScore@meta.data$target8,
  seuset_AddModuleScore@meta.data$target9,
  seuset_AddModuleScore@meta.data$target10,
  seuset_AddModuleScore@meta.data$target11,
  seuset_AddModuleScore@meta.data$target12,
  seuset_AddModuleScore@meta.data$target13,
  seuset_AddModuleScore@meta.data$target14,
  seuset_AddModuleScore@meta.data$target15
  )

rownames(type) <- Cells(seuset_AddModuleScore)

colnames(type) <- c("target1",
                    "target2",
                    "target3",
                    "target4",
                    "target5",
                    "target6",
                    "target7",
                    "target8",
                    "target9",
                    "target10",
                    "target11",
                    "target12",
                    "target13",
                    "target14",
                    "target15")

type$max <- do.call(`pmax`, type)

j1 <- max.col(type, "first")
type$type <- c("target1",
               "target2",
               "target3",
               "target4",
               "target5",
               "target6",
               "target7",
               "target8",
               "target9",
               "target10",
               "target11",
               "target12",
               "target13",
               "target14",
               "target15")[j1]

type$type2 <- ifelse(type$max < quantile(type$max, 0.96), "other" ,type$type)
type$type3 <- ifelse(type$target7 < quantile(type$target7, 0.96), "other" ,"target7")

pdf("./Output_new4_30/Try15_AddmoduleScore_target/DimPlot.pdf", useDingbats = F, height = 4.5, width = 14)
seuset_AddModuleScore <- AddMetaData(seuset_AddModuleScore, metadata = type$type2, col.name = "subtype")
DimPlot(object = seuset_AddModuleScore, group.by = 'subtype', split.by = "group")
seuset_AddModuleScore <- AddMetaData(seuset_AddModuleScore, metadata = type$type3, col.name = "Osm_targets")
DimPlot(object = seuset_AddModuleScore, group.by = 'Osm_targets', split.by = "group")
dev.off()

# Il6 vs OsM
list.il6_osm <- merge(x=Osm , y=Il6, by.x=1, by.y=1)

#
markers_RA <- read.table("./Output_new4_30/Try3_AddModuleScore_fibro/mmu05323_KEGG_RA.txt", header=T, sep="\t")
# AddModuleScore
markers_RA_list <- list(markers_RA$RA, markers_RA$Inflammatory[1:8], markers_RA$Destructive[1:8])
seuset_AddModuleScore <- AddModuleScore(
  object = seuset,
  features = markers_RA_list, ctrl = 100,
  name = 'RA_score')

pdf("./Output_new4_30/Try15_AddmoduleScore_target/RA_score1.pdf", useDingbats = F, height = 3, width = 5)
FeaturePlot(seuset_AddModuleScore, features = c("RA_score1"),
            cols = c("lightgrey", "#005493"),
            min.cutoff = "q10", max.cutoff = "q90")
FeaturePlot(seuset_AddModuleScore, features = c("RA_score2"),
            cols = c("lightgrey", "#005493"),
            min.cutoff = "q10", max.cutoff = "q90")
FeaturePlot(seuset_AddModuleScore, features = c("RA_score3"),
            cols = c("lightgrey", "#005493"),
            min.cutoff = "q10", max.cutoff = "q90")
VlnPlot(seuset_AddModuleScore, features = c("RA_score1"),pt.size = 0.1, log = T) + NoLegend()
VlnPlot(seuset_AddModuleScore, features = c("RA_score2"),pt.size = 0.1, log = T) + NoLegend()
VlnPlot(seuset_AddModuleScore, features = c("RA_score3"),pt.size = 0.1, log = T) + NoLegend()
dev.off()

pdf("./Output_new4_30/Try15_AddmoduleScore_target/RA_score2.pdf", useDingbats = F, height = 3, width = 7)
VlnPlot(seuset_AddModuleScore, features = c("RA_score1"), pt.size = 0, cols = c("#078992","#F8766D","#D39200"),log = T, split.by = "group")
VlnPlot(seuset_AddModuleScore, features = c("RA_score2"), pt.size = 0, cols = c("#078992","#F8766D","#D39200"),log = T, split.by = "group")
VlnPlot(seuset_AddModuleScore, features = c("RA_score3"), pt.size = 0, cols = c("#078992","#F8766D","#D39200"),log = T, split.by = "group")
dev.off()

seuset_AddModuleScore$clusters <- Idents(seuset_AddModuleScore)
Idents(seuset_AddModuleScore) <- "group"
CIA <- subset(seuset_AddModuleScore, idents= c("A_Ctr","B_CIA"))
Idents(seuset_AddModuleScore) <- "clusters"
Idents(CIA) <- "clusters"

pdf("./Output_new4_30/Try15_AddmoduleScore_target/RA_score3.pdf", useDingbats = F, height = 3, width = 6)
VlnPlot(CIA, features = c("RA_score1"), pt.size = 0, cols = c("#078992","#F8766D"), log = T, split.by = "group")
VlnPlot(CIA, features = c("RA_score2"), pt.size = 0, cols = c("#078992","#F8766D"), log = T, split.by = "group")
VlnPlot(CIA, features = c("RA_score3"), pt.size = 0, cols = c("#078992","#F8766D"), log = T, split.by = "group")
dev.off()
