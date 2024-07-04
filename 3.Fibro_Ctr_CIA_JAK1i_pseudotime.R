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
library(annotables)
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
Fibro <- readRDS("./Save_new/Fibro_Ctr_CIA_JAK1i_30_2.RDS")

#Converting to SingleCellExperiment
sce <- as.SingleCellExperiment(Fibro)
Fibro_sub <- subset(Fibro, idents = c("Fibro_1","Fibro_2","Fibro_3","Fibro_4"))
sce_sub <- as.SingleCellExperiment(Fibro_sub)

# Slingshot
sce1 <- slingshot(sce_sub, clusterLabels = 'ident', reducedDim = "UMAP", start.clus="Fibro_1")

# Color
colors <- colorRampPalette(brewer.pal(5,'Spectral'))(100)
colors1 <- viridis(5)

# SlingshotDataSet
UMAP <- data.frame(reducedDims(sce)$UMAP)
PCA <- data.frame(reducedDims(sce)$PCA)

# Plot
pdf("./Output_new4_30/Try2_Fibro_Ctr_CIA_JAK1i_pseudotime/slingshot.pdf", useDingbats = F, height = 4, width = 4.5)
plot(reducedDims(sce)$UMAP, col = colors1[sce$ident], pch=20, asp = 1, frame = FALSE, cex.axis = 1.5, cex.lab = 1.5)
lines(SlingshotDataSet(sce1), lwd=2,type = 'lineages', col = "black")
par(las = 1, tck = -0.01)
dev.off()

pdf("./Output_new4_30/Try2_Fibro_Ctr_CIA_JAK1i_pseudotime/slingshot2.pdf", useDingbats = F, height = 3, width = 5)
ggplot(UMAP, aes(UMAP_1, UMAP_2)) +
  geom_point(col= "grey") +
  geom_point(data= data.frame(reducedDims(sce1)$UMAP), colour= colors[cut(sce1$slingPseudotime_1, breaks = 100)]) +
#  geom_point(data= data.frame(SlingshotDataSet(sce1)@curves[["curve1"]][["s"]]), colour= "black", size = 2) +
  theme(aspect.ratio = 1, panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        text = element_text(size=20), axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"))

dev.off()

# Identifying temporally expressed genes
t <- sce1$slingPseudotime_1
# Select makers of Pseudotime_1
markers_<- read.table("./Output_new4_30/Try14_Fibro_Ctr_CIA_JAK1i/markers_0.3.txt", header = T, sep = "\t")
diff_genes <- subset(markers_$gene, markers_$cluster %in% c("Fibro_1","Fibro_2","Fibro_3","Fibro_4"))
diff_genes <- unique(diff_genes)
X <- as.matrix(GetAssayData(object = Fibro_sub))
Y <- subset(X, rownames(X) %in% diff_genes)
# fit a GAM with a loess term for pseudotime
gam.pval <- apply(Y,1,function(z){
  d <- data.frame(z=z, t=t)
  tmp <- gam(z ~ lo(t), data=d)
  p <- summary(tmp)[4][[1]][1,5]
  p
})

differ_genes <- as.data.frame(gam.pval)
differ_genes$gene <- rownames(differ_genes)
differ_genes <- subset(differ_genes, differ_genes$gam.pval < 0.05)
differ_genes <- dplyr ::arrange(differ_genes, differ_genes$gam.pval)
write.table(differ_genes,"./Output_new4_30/Try2_Fibro_Ctr_CIA_JAK1i_pseudotime/differ_genes.txt", sep="\t", quote=FALSE, row.names=FALSE)  

# Loading required package: clusterExperiment
topgenes <- differ_genes$gene[1:100]
heatdata <- X[rownames(X) %in% topgenes, order(t, na.last = NA)]
heatclus <- sce1$ident[order(t, na.last = NA)]
ce <- ClusterExperiment(heatdata, heatclus)
# Plot
pdf("./Output_new4_30/Try2_Fibro_Ctr_CIA_JAK1i_pseudotime/Genes_trajectory_1.pdf", useDingbats = F, height = 10, width = 6)
plotHeatmap(ce, clusterSamplesData = "orderSamplesValue",visualizeData = 'transformed') + theme(aspect.ratio = 1)
dev.off()

# plotExpression
pdf("./Output_new4_30/Try2_Fibro_Ctr_CIA_JAK1i_pseudotime/Genes_trajectory_2.pdf", useDingbats = F, height = 6, width = 6)
plotE(sce1, c("Thy1","Cd34","Cd55","Prg4"), x = "slingPseudotime_1", 
      colour_by = "ident", show_violin = FALSE,
      show_smooth = TRUE) +
  labs(x="Pseudotime") +
  labs(fill='Cluster') +
  scale_fill_viridis(discrete=TRUE) +
  theme(aspect.ratio = 0.7, text = element_text(size=15), axis.text.x = element_text(size = 15),
        axis.text.y = element_text( size = 15), strip.background = element_rect(fill = NA))

plotE(sce1, c("Cdh11","Vcam1","Cd276","Sdc1"), x = "slingPseudotime_1", 
      colour_by = "ident", show_violin = FALSE,
      show_smooth = TRUE) +
  labs(x="Pseudotime") +
  labs(fill='Cluster') +
  scale_fill_viridis(discrete=TRUE) +
  theme(aspect.ratio = 0.7, text = element_text(size=15), axis.text.x = element_text(size = 15),
        axis.text.y = element_text( size = 15), strip.background = element_rect(fill = NA))

plotE(sce1, c("Ets1","Ets2","Runx1","Tnfsf11"), x = "slingPseudotime_1", 
      colour_by = "ident", show_violin = FALSE,
      show_smooth = TRUE) +
  labs(x="Pseudotime") +
  labs(fill='Cluster') +
  scale_fill_viridis(discrete=TRUE) +
  theme(aspect.ratio = 0.7, text = element_text(size=15), axis.text.x = element_text(size = 15),
        axis.text.y = element_text( size = 15), strip.background = element_rect(fill = NA))

dev.off()
