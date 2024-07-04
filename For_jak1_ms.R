# Set working directory
setwd("/Volumes/Nam_HDD/Ref data (sc_bulk RNAseq)/CIA")

library(dplyr)
library(plyr)

# Fig2C
Fig2c_data <- read.table("./Output_new4_30/Try19_jak1_ms/Fig2c.txt", sep= "\t",header = TRUE)

Fig2c_data$Group <- factor(Fig2c_data$Group, levels = c("CIA+JAK1i","CIA"))
Fig2c_data$Subtype <- factor(Fig2c_data$Subtype, levels = rev(c("Fibro_1","Fibro_2","Fibro_3","Fibro_4","Chondrocyte","Myel_a1","Myel_a2","Myel_a3","Myel_b","Myel_c1","Myel_c2",
                                                            "Myel_c3","Myel_c4","Myel_c5")))

Fig2c_data <- ddply(Fig2c_data, "Subtype",
                   transform, label_ypos=cumsum(Significant_gene_number)- 0.5*Significant_gene_number)

head(Fig2c_data)

pdf("./Output_new4_30/Try19_jak1_ms/Fig2c.pdf", useDingbats = F, height = 8, width = 6)
ggplot(data=Fig2c_data, aes(x=Subtype, y=Significant_gene_number, fill=Group)) +
  geom_bar(stat="identity")+
  geom_text(aes(y=label_ypos, label=Significant_gene_number), vjust=0.5, 
            color="black", size=3.5)+
  scale_fill_brewer(palette="Paired")+
  theme_classic() + coord_flip()
dev.off()

# Fig3A
Fig3a_data <- read.table("./Output_new4_30/Try19_jak1_ms/Fig3a.txt", sep= "\t",header = TRUE)

Fig3a_data$Group <- factor(Fig3a_data$Group, levels = c("CIA+JAK1i","CIA"))
Fig3a_data$Subtype.Subtype <- factor(Fig3a_data$Subtype.Subtype, levels = rev(c('Myel_b/Fibro_2','Myel_b/Fibro_3','Myel_c1/Fibro_2',
'Myel_c1/Fibro_3','Myel_c4/Fibro_2','Myel_c4/Fibro_3')))

pdf("./Output_new4_30/Try19_jak1_ms/Fig3a.pdf", useDingbats = F, height = 2, width = 3.5)
ggplot(data=Fig3a_data, aes(x=Subtype.Subtype, y=Mean_number_of_interaction, fill=Group)) +
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_brewer(palette="Paired")+
  theme_classic() + coord_flip()
dev.off()

# Fig3D
merge.integrated_3 <- readRDS("./Save_new/seuset_merge_Ctr_CIA_JAK1i_30_2.RDS")
Extract <- DotPlot(merge.integrated_3, features = c("Il6ra", "Il6st", "Osmr"), idents = c("Myeloid_b","Myeloid_c","Fibroblast"))
Extract$data






# Fig4A
Fig4a_data <- read.table("./Output_new4_30/Try19_jak1_ms/Fig4a.txt", sep= "\t",header = TRUE)

Fig4a_data$Subtype.Signaling <- factor(Fig4a_data$Subtype.Signaling, levels = rev(c('Fibro_2/WP_ONCOSTATIN_M_SIGNALING_PATHWAY','Myel_c2/HALLMARK_IL6_JAK_STAT3_SIGNALING',
                                                                                    'Chondrocyte/WP_ONCOSTATIN_M_SIGNALING_PATHWAY','Fibro_3/WP_ONCOSTATIN_M_SIGNALING_PATHWAY',
                                                                                    'Fibro_1/WP_ONCOSTATIN_M_SIGNALING_PATHWAY','Fibro_4/WP_ONCOSTATIN_M_SIGNALING_PATHWAY',
                                                                                    'Chondrocyte/HALLMARK_INTERFERON_GAMMA_RESPONSE','Fibro_1/HALLMARK_INTERFERON_GAMMA_RESPONSE',
                                                                                    'Fibro_3/HALLMARK_IL6_JAK_STAT3_SIGNALING','Fibro_3/HALLMARK_INTERFERON_GAMMA_RESPONSE',
                                                                                    'Fibro_1/HALLMARK_IL6_JAK_STAT3_SIGNALING')))

pdf("./Output_new4_30/Try19_jak1_ms/Fig4a.pdf", useDingbats = F, height = 3, width = 6)
ggplot(data=Fig4a_data, aes(x=Subtype.Signaling, y=NES)) +
  geom_bar(stat="identity", fill="steelblue")+
  theme_classic() + coord_flip()
dev.off()

# FigS5
FigS5_data <- read.table("./Output_new4_30/Try19_jak1_ms/FigS5.txt", sep= "\t",header = TRUE)

FigS5_data$Group <- factor(FigS5_data$Group, levels = c("CIA+JAK1i","CIA"))
FigS5_data$Subtype <- factor(FigS5_data$Subtype, levels = rev(c('T_cell', 'B_cell', 'Myel_a', 'Myel_b', 'Myel_c', 'Fibroblast', 'Neutrophil', 'Endothelial', 'Mural_cell')))

FigS5_data <- ddply(FigS5_data, "Subtype",
                    transform, label_ypos=cumsum(Significant_gene_number)- 0.5*Significant_gene_number)

head(FigS5_data)

pdf("./Output_new4_30/Try19_jak1_ms/FigS5.pdf", useDingbats = F, height = 4, width = 6)
ggplot(data=FigS5_data, aes(x=Subtype, y=Significant_gene_number, fill=Group)) +
  geom_bar(stat="identity")+
  geom_text(aes(y=label_ypos, label=Significant_gene_number), vjust=0.5, 
            color="black", size=3.5)+
  scale_fill_brewer(palette="Paired")+
  theme_classic() + coord_flip()
dev.off()







