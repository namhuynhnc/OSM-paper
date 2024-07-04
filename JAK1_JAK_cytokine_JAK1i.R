# Set working directory
setwd("/Volumes/Nam_HDD/Ref data (sc_bulk RNAseq)/CIA")

library(scater)
library(scran)
library(SingleCellExperiment)
library(Biobase)
library(irlba)
library(ggplot2)
library(gplots)
library(annotables)
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

#DotPlot
DotPlot(seuset, features = rev(c("Ifng","Il22","Il10","Ifnb1","Ifna2","Osm","Il6","Il21","Il15","Il7","Il4","Il2","Csf2","Epo","Gh"))) + RotatedAxis()
DotPlot(seuset, features = rev(c("Tyk2","Jak3","Jak2","Jak1"))) + RotatedAxis()
FeaturePlot(
  seuset, 
  c("Jak1"),
  cols = c("lightgrey", "#005493")
  # , min.cutoff = "q10", max.cutoff = "q90"
)

all <- as.data.frame(matrix(0, nrow= length(levels(seuset@active.ident)), ncol=15))
colnames(all) <- c("Gh","Epo","Csf2","Ifng","Il22","Il10","Ifnb","Ifna","Osm","Il6","Il21","Il15","Il7","Il4","Il2")
rownames(all) <- 1: length(levels(seuset@active.ident))

#Load marker
#markers <- read.table("./Output_new4_30/Try13_merge_Ctr_CIA_JAK1i/markers_celltype_CIA.txt", sep="\t", header = T)
markers <- read.table("./Output_new4_30/Try13_merge_Ctr_CIA_JAK1i/markers_celltype_JAK1i.txt", sep="\t", header = T)
markers <- subset(markers, markers$avg_logFC>0 
                  #                  & markers$pct.1>0.25
)

cluster_name <- data.frame(summary(seuset@active.ident))
cluster_name$cluster <- rownames(cluster_name)
cluster_name$number <- c(1: length(levels(seuset@active.ident)))

markers$cluster_number <- cluster_name$number[match(unlist(markers$cluster), cluster_name$cluster)]
table(markers$cluster)

#Match mouse and human genes
load("/Volumes/Nam_HDD/Ref data (sc_bulk RNAseq)/Ref/Single cell RNAseq Immunity 2019/Data/7098575/bm.query_mouse_allevo.RData")
genelist <- data.frame(bm.query_mouse_allevo_list[[1]])
genelist2 <- data.frame(bm.query_mouse_allevo_list[[5]]) # match human gene
genelist3 <- merge(genelist, genelist2, by.x=1, by.y=1)
genelist3[genelist3==""]<-NA
genelist3 <- na.omit(genelist3)
genelist3 <-  genelist3[,c(1:4)]
#genelist4 <- genelist3[!duplicated(genelist3$hsapiens_homolog_ensembl_gene),] # remove duplicated genes

#Match gene accession ID
accession <- read_excel("/Volumes/Nam_HDD/Ref data (sc_bulk RNAseq)/phs001529.v1.p1/Results/JAK1-cytokine/Ref/IL10/gene_accession.xlsx",sheet=1, col_names = T)

#IL2 (human)
list <- read_excel("/Volumes/Nam_HDD/Ref data (sc_bulk RNAseq)/phs001529.v1.p1/Results/JAK1-cytokine/Ref/IL2/dxh283_Supplementary_Data/dxh283_SupplementalTable2.xls",sheet=1, col_names = TRUE)
list <- dplyr :: arrange(list, desc(list$`log10(+/-)`))
list <- data.frame(list$`log10(+/-)`, list$`Gene Symbol`)
list <- remove.factors(list)
list <- subset(list, list[,1]>0)
listn <- na.omit(list)
colnames(listn) <- c("log10FC","gene")

#Match mouse gene
listn <- merge(listn, genelist3, by.x="gene", by.y=4)
listn <- listn[!duplicated(listn$external_gene_name),] # remove duplicated genes

#IL2 (mouse)
listm <- read_excel("/Volumes/Nam_HDD/Ref data (sc_bulk RNAseq)/phs001529.v1.p1/Results/JAK1-cytokine/Ref/IL15_21/40425_2014_28_MOESM1_ESM.xlsx",sheet=2, col_names = T)
listm <- listm[,c(9,10)]
listm <- subset(listm, listm[[2]]>1)

#IL2 target list
list <- data.frame(c(listn$external_gene_name, listm$Gene...9))
colnames(list) <- "gene"
Il2_ <- list[!duplicated(list$gene),] # remove duplicated genes
list <- as.data.frame(list[!duplicated(list$gene),]) # remove duplicated genes
colnames(list) <- "gene"

#Save txt
gene_numbers <- list()
for(i in 1: length(levels(seuset@active.ident))){
  #match
  list.merge  <- merge(x=subset(markers, markers$cluster_number %in% i) , y=list, by.x="gene", by.y="gene")
  gene_numbers[[i]] <- nrow(list.merge)/length(Il2_)
  print(gene_numbers[[i]])
  write.table(list.merge, paste0("./Output_new4_30/Try9_JAK_JAK1i/IL2/IL2.",i,".txt"), sep="\t", quote=FALSE, row.names=FALSE)
}

all$Il2 <- gene_numbers

#IL4 (human)
list <- read_excel("/Volumes/Nam_HDD/Ref data (sc_bulk RNAseq)/phs001529.v1.p1/Results/JAK1-cytokine/Ref/IL4/ScienceDirect_files_09Jan2020_09-22-33.657/1-s2.0-S1074761310002141-mmc2.xls",sheet=2, col_names = F)
colnames(list) <- list[6,]
list <- list[-c(1:6),]
list <- subset(list, list$`FDR F-test`<0.05)
list <- data.frame(list$`Symbol (IPA)`)
list <- remove.factors(list)
list <- na.omit(list)
listn <- data.frame(list[!duplicated(list),])
colnames(listn) <- c("gene")

#Match mouse gene
listn <- merge(listn, genelist3, by.x="gene", by.y=4)
listn <- listn[!duplicated(listn$external_gene_name),] # remove duplicated genes

#IL4 (mouse)
listm <- read_excel("/Volumes/Nam_HDD/Ref data (sc_bulk RNAseq)/phs001529.v1.p1/Results/JAK1-cytokine/Ref/IL15_21/40425_2014_28_MOESM1_ESM.xlsx",sheet=2, col_names = T)
listm <- listm[,c(9,11)]
listm <- subset(listm, listm[[2]]>1)

#IL4 target list
list <- data.frame(c(listn$external_gene_name, listm$Gene...9))
colnames(list) <- "gene"
Il4_ <- list[!duplicated(list$gene),] # remove duplicated genes
list <- as.data.frame(list[!duplicated(list$gene),]) # remove duplicated genes
colnames(list) <- "gene"

#Save txt
gene_numbers <- list()
for(i in 1: length(levels(seuset@active.ident))){
  #match
  list.merge  <- merge(x=subset(markers, markers$cluster_number %in% i) , y=list, by.x="gene", by.y="gene")
  gene_numbers[[i]] <- nrow(list.merge)/length(Il4_)
  print(gene_numbers[[i]])
  write.table(list.merge, paste0("./Output_new4_30/Try9_JAK_JAK1i/IL4/IL4.",i,".txt"), sep="\t", quote=FALSE, row.names=FALSE)
}
all$Il4 <- gene_numbers

#IL7 (mouse)
listn <- read_excel("/Volumes/Nam_HDD/Ref data (sc_bulk RNAseq)/phs001529.v1.p1/Results/JAK1-cytokine/Ref/IL7/IL7_IFNg.xls",sheet=1, col_names = T)
listn <- subset(listn, listn$log2FC>0)

#IL7 (mouse2)
listm <- read_excel("/Volumes/Nam_HDD/Ref data (sc_bulk RNAseq)/phs001529.v1.p1/Results/JAK1-cytokine/Ref/IL15_21/40425_2014_28_MOESM1_ESM.xlsx",sheet=2, col_names = T)
listm <- listm[,c(9,12)]
listm <- subset(listm, listm[[2]]>1)

#IL7 target list
list <- data.frame(c(listn$`Gene ID`, listm$Gene...9))
colnames(list) <- "gene"
Il7_ <- list[!duplicated(list$gene),] # remove duplicated genes
list <- as.data.frame(list[!duplicated(list$gene),]) # remove duplicated genes
colnames(list) <- "gene"

#Save txt
gene_numbers <- list()
for(i in 1: length(levels(seuset@active.ident))){
  #match
  list.merge  <- merge(x=subset(markers, markers$cluster_number %in% i) , y=list, by.x="gene", by.y="gene")
  gene_numbers[[i]] <- nrow(list.merge)/length(Il7_)
  print(gene_numbers[[i]])
  write.table(list.merge, paste0("./Output_new4_30/Try9_JAK_JAK1i/IL7/IL7.",i,".txt"), sep="\t", quote=FALSE, row.names=FALSE)
}
all$Il7 <- gene_numbers

#IFNg (human)
list <- read_excel("/Volumes/Nam_HDD/Ref data (sc_bulk RNAseq)/phs001529.v1.p1/Results/JAK1-cytokine/Ref/IFN/interferome.org.2171-2169-2170.xlsx",sheet=1, col_names = T)
list <- list[,c(2,5)]
list <- subset(list, list$`Fold Change`>1)
list <- data.frame(list[!duplicated(list[[2]]),])
colnames(list)[[2]] <- c("gene")

#Match mouse gene
list <- merge(list, genelist3, by.x="gene", by.y=4)
list <-as.data.frame(list$external_gene_name)
colnames(list) <- "gene"
Ifng_ <- list[!duplicated(list$gene),] # remove duplicated genes
list <- as.data.frame(list[!duplicated(list$gene),]) # remove duplicated genes
colnames(list) <- "gene"

#Save txt
gene_numbers <- list()
for(i in 1: length(levels(seuset@active.ident))){
  #match
  list.merge  <- merge(x=subset(markers, markers$cluster_number %in% i) , y=list, by.x="gene", by.y="gene")
  gene_numbers[[i]] <- nrow(list.merge)/length(Ifng_)
  print(gene_numbers[[i]])
  write.table(list.merge, paste0("./Output_new4_30/Try9_JAK_JAK1i/IFNg/IFNg.",i,".txt"), sep="\t", quote=FALSE, row.names=FALSE)
}
all$Ifng <- gene_numbers

#IFNa (human)
list <- read_excel("/Volumes/Nam_HDD/Ref data (sc_bulk RNAseq)/phs001529.v1.p1/Results/JAK1-cytokine/Ref/IFN/interferome.org.2171-2169-2170.xlsx",sheet=2, col_names = T)
list <- list[,c(2,5)]
list <- subset(list, list$`Fold Change`>1)
list <- data.frame(list[!duplicated(list[[2]]),])
colnames(list)[[2]] <- c("gene")

#Match mouse gene
list <- merge(list, genelist3, by.x="gene", by.y=4)
list <-as.data.frame(list$external_gene_name)
colnames(list) <- "gene"
Ifna_ <- list[!duplicated(list$gene),] # remove duplicated genes
list <- as.data.frame(list[!duplicated(list$gene),]) # remove duplicated genes
colnames(list) <- "gene"

#Save txt
gene_numbers <- list()
for(i in 1: length(levels(seuset@active.ident))){
  #match
  list.merge  <- merge(x=subset(markers, markers$cluster_number %in% i) , y=list, by.x="gene", by.y="gene")
  gene_numbers[[i]] <- nrow(list.merge)/length(Ifna_)
  print(gene_numbers[[i]])
  write.table(list.merge, paste0("./Output_new4_30/Try9_JAK_JAK1i/IFNa/IFNa.",i,".txt"), sep="\t", quote=FALSE, row.names=FALSE)
}
all$Ifna <- gene_numbers

#IFNb (human)
list <- read_excel("/Volumes/Nam_HDD/Ref data (sc_bulk RNAseq)/phs001529.v1.p1/Results/JAK1-cytokine/Ref/IFN/interferome.org.2171-2169-2170.xlsx",sheet=3, col_names = T)
list <- list[,c(2,5)]
list <- subset(list, list$`Fold Change`>1)
list <- data.frame(list[!duplicated(list[[2]]),])
colnames(list)[[2]] <- c("gene")

#Match mouse gene
list <- merge(list, genelist3, by.x="gene", by.y=4)
list <-as.data.frame(list$external_gene_name)
colnames(list) <- "gene"
Ifnb_ <- list[!duplicated(list$gene),] # remove duplicated genes
list <- as.data.frame(list[!duplicated(list$gene),]) # remove duplicated genes
colnames(list) <- "gene"

#Save txt
gene_numbers <- list()
for(i in 1: length(levels(seuset@active.ident))){
  #match
  list.merge  <- merge(x=subset(markers, markers$cluster_number %in% i) , y=list, by.x="gene", by.y="gene")
  gene_numbers[[i]] <- nrow(list.merge)/length(Ifnb_)
  print(gene_numbers[[i]])
  write.table(list.merge, paste0("./Output_new4_30/Try9_JAK_JAK1i/IFNb/IFNb.",i,".txt"), sep="\t", quote=FALSE, row.names=FALSE)
}
all$Ifnb <- gene_numbers

#IL10 (mouse)
list <- read_excel("/Volumes/Nam_HDD/Ref data (sc_bulk RNAseq)/phs001529.v1.p1/Results/JAK1-cytokine/Ref/IL10/IL10.xlsx",sheet=1, col_names = T)
list <- subset(list, list$FC>1)
#Match accession ID
list <- merge(list, accession, by.x=1, by.y=2)
list <- na.omit(list)
list <- list[!duplicated(list$Gene),] # remove duplicated genes

list <- data.frame(c(list$Gene))
colnames(list) <- "gene"
Il10_ <- list$gene

#Save txt
gene_numbers <- list()
for(i in 1: length(levels(seuset@active.ident))){
  #match
  list.bind  <- merge(x=subset(markers, markers$cluster_number %in% i) , y=list, by.x="gene", by.y="gene")
  gene_numbers[[i]] <- nrow(list.bind)/length(Il10_)
  print(gene_numbers[[i]])
  write.table(list.bind, paste0("./Output_new4_30/Try9_JAK_JAK1i/IL10/IL10.",i,".txt"), sep="\t", quote=FALSE, row.names=FALSE)
}
all$Il10 <- gene_numbers

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

list <- data.frame(c(list1$Rank, list2$Symbol, list3$Symbol))
colnames(list) <- "Gene.symbol"
list[list=="NA"]<-NA
list <- na.omit(list)
list <- data.frame(list[!duplicated(list$Gene.symbol),]) # remove duplicated genes
colnames(list) <- "gene"
Osm_ <- list$gene

#Save txt
gene_numbers <- list()
for(i in 1: length(levels(seuset@active.ident))){
  #match
  list.merge <- merge(x=subset(markers, markers$cluster_number %in% i) , y=list, by.x="gene", by.y="gene")
  gene_numbers[[i]] <- nrow(list.merge)/length(Osm_)
  print(gene_numbers[[i]])
  write.table(list.merge, paste0("./Output_new4_30/Try9_JAK_JAK1i/OSM/OSM.",i,".txt"), sep="\t", quote=FALSE, row.names=FALSE)
}
all$Osm <- gene_numbers

#IL6 (mouse)
list <- read.table("/Volumes/Nam_HDD/Ref data (sc_bulk RNAseq)/phs001529.v1.p1/Results/JAK1-cytokine/Ref/IL6/result.txt",sep="\t", header = T)
list <- subset(list, list$adj.P.Val<0.05)
list <- subset(list, list$logFC>0)
list[list==""]<-NA
list <- na.omit(list)
list <- list[!duplicated(list$Gene.symbol),] # remove duplicated genes

list <- data.frame(c(list$Gene.symbol))
colnames(list) <- "gene"
Il6_ <- list$gene

#Save txt
gene_numbers <- list()
for(i in 1: length(levels(seuset@active.ident))){
  #match
  list.merge <- merge(x=subset(markers, markers$cluster_number %in% i) , y=list, by.x="gene", by.y="gene")
  gene_numbers[[i]] <- nrow(list.merge)/length(Il6_)
  print(gene_numbers[[i]])
  write.table(list.merge, paste0("./Output_new4_30/Try9_JAK_JAK1i/IL6/IL6.",i,".txt"), sep="\t", quote=FALSE, row.names=FALSE)
}
all$Il6 <- gene_numbers

#IL22 (mouse)
list <- read_excel("/Volumes/Nam_HDD/Ref data (sc_bulk RNAseq)/phs001529.v1.p1/Results/JAK1-cytokine/Ref/IL22/journal.pbio.3000540.s006.xlsx",sheet=1, col_names = T)
colnames(list) <- list[1,]
list <- list[-1,c(1,7,9)]
list$logFC <- as.numeric(list$logFC)
list$adj.P.Val <- as.numeric(list$adj.P.Val)
list <- subset(list, list$adj.P.Val <0.05 & list$logFC>0)

list <- data.frame(c(list$ID))
colnames(list) <- "gene"
Il22_ <- list$gene

#Save txt
gene_numbers <- list()
for(i in 1: length(levels(seuset@active.ident))){
  #match
  list.merge <- merge(x=subset(markers, markers$cluster_number %in% i) , y=list, by.x="gene", by.y="gene")
  gene_numbers[[i]] <- nrow(list.merge)/length(Il22_)
  print(gene_numbers[[i]])
  write.table(list.merge, paste0("./Output_new4_30/Try9_JAK_JAK1i/IL22/IL22.",i,".txt"), sep="\t", quote=FALSE, row.names=FALSE)
}
all$Il22 <- gene_numbers

#IL15 (mouse)
list <- read_excel("/Volumes/Nam_HDD/Ref data (sc_bulk RNAseq)/phs001529.v1.p1/Results/JAK1-cytokine/Ref/IL15_21/40425_2014_28_MOESM1_ESM.xlsx",sheet=2, col_names = T)
list <- list[,c(9,13)]
list <- subset(list, list[[2]]>1)

list <- data.frame(c(list$Gene...9))
colnames(list) <- "gene"
Il15_ <- list$gene

#Save txt
gene_numbers <- list()
for(i in 1: length(levels(seuset@active.ident))){
  #match
  list.merge <- merge(x=subset(markers, markers$cluster_number %in% i) , y=list, by.x="gene", by.y="gene")
  gene_numbers[[i]] <- nrow(list.merge)/length(Il15_)
  print(gene_numbers[[i]])
  write.table(list.merge, paste0("./Output_new4_30/Try9_JAK_JAK1i/IL15/IL15.",i,".txt"), sep="\t", quote=FALSE, row.names=FALSE)
}
all$Il15 <- gene_numbers

#IL21 (mouse)
list <- read_excel("/Volumes/Nam_HDD/Ref data (sc_bulk RNAseq)/phs001529.v1.p1/Results/JAK1-cytokine/Ref/IL15_21/40425_2014_28_MOESM1_ESM.xlsx",sheet=2, col_names = T)
list <- list[,c(9,14)]
list <- subset(list, list[[2]]>1)

list <- data.frame(c(list$Gene...9))
colnames(list) <- "gene"
Il21_ <- list$gene

#Save txt
gene_numbers <- list()
for(i in 1: length(levels(seuset@active.ident))){
  #match
  list.merge <- merge(x=subset(markers, markers$cluster_number %in% i) , y=list, by.x="gene", by.y="gene")
  gene_numbers[[i]] <- nrow(list.merge)/length(Il21_)
  print(gene_numbers[[i]])
  write.table(list.merge, paste0("./Output_new4_30/Try9_JAK_JAK1i/IL21/IL21.",i,".txt"), sep="\t", quote=FALSE, row.names=FALSE)
}
all$Il21 <- gene_numbers

#GM-CSF (human mouse)
list <- read_excel("/Volumes/Nam_HDD/Ref data (sc_bulk RNAseq)/phs001529.v1.p1/Results/JAK1-cytokine/Ref/GM-CSF/GMCSF.xlsx",sheet=1, col_names = T)

#Match mouse gene
listb <- merge(list[1:100,], genelist3, by.x="gene", by.y=4)
listc <- listb[,c(4,2)]
colnames(listc)[1] <- "gene"
listc <- rbind.data.frame(listc, list[101:200,])
listc <- subset(listc, listc$FC>1)
list <- listc[!duplicated(listc$gene),] # remove duplicated genes
Csf2_ <- list$gene

#Save txt
gene_numbers <- list()
for(i in 1: length(levels(seuset@active.ident))){
  #match
  list.merge <- merge(x=subset(markers, markers$cluster_number %in% i) , y=list, by.x="gene", by.y="gene")
  gene_numbers[[i]] <- nrow(list.merge)/length(Csf2_)
  print(gene_numbers[[i]])
  write.table(list.merge, paste0("./Output_new4_30/Try9_JAK_JAK1i/GM-CSF/GM-CSF.",i,".txt"), sep="\t", quote=FALSE, row.names=FALSE)
}
all$Csf2 <- gene_numbers

# Epo (human)
list <- read.table("/Volumes/Nam_HDD/Ref data (sc_bulk RNAseq)/phs001529.v1.p1/Results/JAK1-cytokine/Ref/EPO/GSE6260.top.table.tsv", sep="\t", header = T)
list <- list[!(list$Gene.symbol == ""), ]
list <- list[!duplicated(list$Gene.symbol),] # remove duplicated genes
list <- subset(list, list$P.Val<0.05 & list$logFC<0)

#Match mouse gene
list <- merge(list, genelist3, by.x="Gene.symbol", by.y=4)
list <- dplyr ::arrange(list, list$logFC) #down in control is up in Epo
list <- list[!duplicated(list$external_gene_name),]

list <- data.frame(list[,10])
colnames(list) <- "gene"
Epo_ <- list$gene 

#Save txt
gene_numbers <- list()
for(i in 1: length(levels(seuset@active.ident))){
  #match
  list.merge <- merge(x=subset(markers, markers$cluster_number %in% i) , y=list, by.x="gene", by.y="gene")
  gene_numbers[[i]] <- nrow(list.merge)/length(Epo_)
  print(gene_numbers[[i]])
  write.table(list.merge, paste0("./Output_new4_30/Try9_JAK_JAK1i/EPO/EPO.",i,".txt"), sep="\t", quote=FALSE, row.names=FALSE)
}
all$Epo <- gene_numbers

# Gh (mouse)
list <- read.table("/Volumes/Nam_HDD/Ref data (sc_bulk RNAseq)/phs001529.v1.p1/Results/JAK1-cytokine/Ref/GH/GSE2120.top.table.tsv", sep="\t", header = T)
list <- subset(list, list$P.Val<0.05 & list$logFC<0) #down in control
list <- list[!(list$Gene.symbol == ""), ]
list <- list[!duplicated(list$Gene.symbol),] # remove duplicated genes

list <-as.data.frame(list$Gene.symbol)
colnames(list) <- "gene"
Gh_ <- list$gene

#Save txt
gene_numbers <- list()
for(i in 1: length(levels(seuset@active.ident))){
  #match
  list.merge <- merge(x=subset(markers, markers$cluster_number %in% i) , y=list, by.x="gene", by.y="gene")
  gene_numbers[[i]] <- nrow(list.merge)/length(Gh_)
  print(gene_numbers[[i]])
  write.table(list.merge, paste0("./Output_new4_30/Try9_JAK_JAK1i/GH/GH.",i,".txt"), sep="\t", quote=FALSE, row.names=FALSE)
}
all$Gh <- gene_numbers

#Make all data
library(reshape2)
library(ggplot2)
library(tidyr)
library(dplyr)
library(circlize)

all <- as.matrix.data.frame(all)
write.table(all, "./Output_new4_30/Try9_JAK_JAK1i/All/all_scale.txt", sep="\t", quote=FALSE, row.names=T)

all <- read.table("./Output_new4_30/Try9_JAK_JAK1i/All/all_scale.txt", sep="\t", header = T)
all$cell <- rownames(all)

plotDat <- gather(all, key = "Cytokine", value = "Percent_target_gene", -cell)
plotDat$cell <- factor(plotDat$cell,levels = all$cell)
plotDat$Cytokine <- factor(plotDat$Cytokine,levels = rev(colnames(all)[1:15]))

pdf("./Output_new4_30/Try9_JAK_JAK1i/All/All_scale.pdf", height = 5, width = 5)
# Matrix heatmap
ggplot(plotDat, aes(Cytokine, cell, col = Percent_target_gene, fill = Percent_target_gene, label = Percent_target_gene)) +
  geom_tile() +
  #  geom_text(col = "black",size=2) +
  theme_minimal() +
  scale_fill_gradient2(low = "lightgrey", high = "red") +
  #  scale_color_gradient2(low = "lightgrey", high = "blue") +
  scale_y_discrete(breaks=cluster_name$number,
                   labels=cluster_name$cluster) +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10),aspect.ratio = 1)

# chordDiagram
#chordDiagram(as.data.frame(plotDat), annotationTrack = "grid", annotationTrackHeight = c(0.01, 0.01),
#             transparency = 0.5)
#circos.track(track.index = 1, panel.fun = function(x, y) {
#  circos.text(CELL_META$xcenter, CELL_META$ylim[1]+1.2, CELL_META$sector.index, 
#              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 0.6)
#}, bg.border = NA) # here set bg.border to NA is important
dev.off()
