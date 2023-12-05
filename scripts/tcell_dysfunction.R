library(data.table)
library(ggplot2)
library(Matrix)
library(ComplexHeatmap)
library(circlize)
library(stringr)
library(GSVA)
library(BiocParallel)
library(doMC)
library(immunedeconv)
library(GetoptLong)
library(tidyverse)
library(ggpubr)
library(rstatix)
library(extrafontdb)
library(extrafont)
library(preprocessCore)
library(ImmuCellAI)
library(sva)
library(ggsignif)
loadfonts()
registerDoMC(cores=20)

setwd("~/Documents/Misc_Work/Other Work/Miller_Related/CONSTRU/")

source("scripts/all_functions.R")
load("Data/GSE82191_all_info.Rdata")
load("Data/GSE9891_all_info.Rdata")
load("Data/GSEOV3_all_info.Rdata")

final_gse82191_expr_df <- gse82191_out[[13]]
final_gse9891_expr_df <- gse9891_out[[13]]
final_gseov3_expr_df <- gseov3_out[[13]]

common_genes <- intersect(intersect(rownames(final_gse82191_expr_df),rownames(final_gseov3_expr_df)),rownames(final_gse9891_expr_df))
all_expr_df <- as.data.frame(cbind(final_gse82191_expr_df[common_genes,], final_gse9891_expr_df[common_genes,], final_gseov3_expr_df[common_genes,]))

#Need to use combat in Training phase as expression distributions are very different
final_all_expr_df = ComBat(dat=all_expr_df, batch=c(rep("A",ncol(final_gse82191_expr_df)),rep("B",ncol(final_gse9891_expr_df)),rep("C",ncol(final_gseov3_expr_df))), 
                           mod=NULL, par.prior=TRUE, prior.plots=FALSE)

#output <- ImmuCellAI(sample=sample_expression_dataset,data_type="rnaseq",group_tag=1,response_tag=0,customer = 0)

constru_tertiles <- c(gse82191_out[[3]],gse9891_out[[3]],gseov3_out[[3]])
cyt_tertiles <- c(gse82191_out[[2]], gse9891_out[[2]], gseov3_out[[2]])

#Order the ids by CONSTRU scores followed by Cyt scores
order_ids <- order(constru_tertiles, cyt_tertiles)
rev_constru_tertiles <- constru_tertiles[order_ids]
rev_cyt_tertiles <- cyt_tertiles[order_ids]
rev_final_all_expr_df <- final_all_expr_df[,order_ids]

#Get the different levels or combinations for constru and cyt tertiles
unique_constru_tertiles <- unique(rev_constru_tertiles)
unique_cyt_tertiles <- unique(rev_cyt_tertiles)
column_split_values <- NULL
k <- 0
for (i in unique_constru_tertiles)
{
  for (j in unique_cyt_tertiles)
  {
    ids <- which(rev_constru_tertiles==i & rev_cyt_tertiles==j)   
    column_split_values <- c(column_split_values,rep(paste0(k),length(ids)))
    k <- k+1
  }
}

#Get the abundance for each constru tertile with cyt groups
constru_tertile_1 <- which(rev_constru_tertiles==1)
revised_constru_tertile_1_df <- as.data.frame(rbind(rev_cyt_tertiles[constru_tertile_1],rev_final_all_expr_df[,constru_tertile_1]))
revised_constru_tertile_1_df[1,] <- as.factor(as.numeric(revised_constru_tertile_1_df[1,]))
output_tertile_1 <- ImmuCellAI(sample = revised_constru_tertile_1_df, data_type = "microarray", 
                               group_tag = 1, response_tag = 0, customer = 0)
group_fre_tertile_1 <- group_fre

constru_tertile_2 <- which(rev_constru_tertiles==2)
revised_constru_tertile_2_df <- as.data.frame(rbind(rev_cyt_tertiles[constru_tertile_2],rev_final_all_expr_df[,constru_tertile_2]))
revised_constru_tertile_2_df[1,] <- as.factor(as.numeric(revised_constru_tertile_2_df[1,]))
output_tertile_2 <- ImmuCellAI(sample = revised_constru_tertile_2_df, data_type = "microarray", 
                               group_tag = 1, response_tag = 0, customer = 0)
group_fre_tertile_2 <- group_fre

constru_tertile_3 <- which(rev_constru_tertiles==3)
revised_constru_tertile_3_df <- as.data.frame(rbind(rev_cyt_tertiles[constru_tertile_3],rev_final_all_expr_df[,constru_tertile_3]))
revised_constru_tertile_3_df[1,] <- as.factor(as.numeric(revised_constru_tertile_3_df[1,]))
output_tertile_3 <- ImmuCellAI(sample = revised_constru_tertile_3_df, data_type = "microarray", 
                               group_tag = 1, response_tag = 0, customer = 0)
group_fre_tertile_3 <- group_fre

abundance_matrix <- as.data.frame(cbind(t(group_fre_tertile_1)[c(1:24),c(1:3)],t(group_fre_tertile_2)[c(1:24),c(1:3)],t(group_fre_tertile_3)[c(1:24),c(1:3)]))
abundance_matrix <- as.matrix(abundance_matrix)

##Immune Concentrations
ha = HeatmapAnnotation(
  empty = anno_empty(border = FALSE, height = unit(8, "mm")),
  Cyt = anno_block(gp = gpar(fill = c("black","red","darkgreen","black","red","darkgreen","black","red","darkgreen")), 
                   labels = c("Cyt low","Cyt int","Cyt high","Cyt low","Cyt int","Cyt high","Cyt low","Cyt int","Cyt high")),
  simple_anno_size = unit(8, "mm"),
  annotation_name_gp = gpar(fontsize=11, family="sans",col="white"),
  annotation_legend_param = list(directon="horizontal"),
  show_legend=F)
rev_column_split_values <- c("0","1","2","3","4","5","6","7","8")
rev_column_split_values <- as.factor(rev_column_split_values)

col_fun2 = colorRamp2(c(0,0.2),c("white","red"))
ht_immune_conc = Heatmap(matrix=abundance_matrix, col_fun2, 
                         name = "Immune Cell Fractions", column_title = qq("Immune Cell Fractions across Constru Groups"),
                         width = unit(15, "cm"),
                         height = unit(15, "cm"),
                         cluster_columns = F,
                         cluster_rows = T,
                         row_names_centered = F,
                         show_row_names = T,
                         row_labels = rownames(abundance_matrix),
                         row_names_max_width = max_text_width(rownames(abundance_matrix)),
                         column_split = rev_column_split_values,
                         cluster_column_slices = F,
                         show_column_names = F,
                         raster_quality = 2,
                         top_annotation = ha,
                         column_title_rot = 0,
                         column_title_side = "bottom",
                         column_dend_side = "top",
                         column_title_gp = gpar(fontsize=11, family="sans"),
                         row_names_gp = gpar(fontsize=11, family="sans"),
                         use_raster = T,
                         column_gap = unit(1, "mm"),
                         border_gp = gpar(col = "black", lty = 1),
                         border = T,
                         heatmap_legend_param = list(direction = "horizontal"))
#pdf("results/Training_ImmuneCellAI_Fractions_Constru_Cyt_Tertiles.pdf", height = 8, width=8, fonts = "sans")
draw(ht_immune_conc, heatmap_legend_side="bottom")
group_block_anno = function(group, empty_anno, gp = gpar(), 
                            label = NULL, label_gp = gpar()) {
  
  seekViewport(qq("annotation_@{empty_anno}_@{min(group)}"))
  loc1 = deviceLoc(x = unit(0, "npc"), y = unit(0, "npc"))
  seekViewport(qq("annotation_@{empty_anno}_@{max(group)}"))
  loc2 = deviceLoc(x = unit(1, "npc"), y = unit(1, "npc"))
  
  seekViewport("global")
  grid.rect(loc1$x, loc1$y, width = loc2$x - loc1$x, height = loc2$y - loc1$y, 
            just = c("left", "bottom"), gp = gp)
  if(!is.null(label)) {
    grid.text(label, x = (loc1$x + loc2$x)*0.5, y = (loc1$y + loc2$y)*0.5, gp = label_gp)
  }
}
group_block_anno(1:3, "empty", gp = gpar(fill = "brown"), label = "Con low")
group_block_anno(4:6, "empty", gp = gpar(fill = "orange"), label = "Con int")
group_block_anno(7:9, "empty", gp = gpar(fill = "#FEDD00"), label = "Con high")
#dev.off()

#Get information about each constru tertile and compare cyt 1 vs cyt 3
constru_tertile_1_df <- as.data.frame(output_tertile_1$Sample_abundance)
constru_tertile_2_df <- as.data.frame(output_tertile_2$Sample_abundance)
constru_tertile_3_df <- as.data.frame(output_tertile_3$Sample_abundance)

#Remove ids from cyt tertile 2
cyt_tertile_2_in_constru_tertile_1 <- which(rev_cyt_tertiles[constru_tertile_1]==2)
cyt_tertile_2_in_constru_tertile_2 <- which(rev_cyt_tertiles[constru_tertile_2]==2)
cyt_tertile_2_in_constru_tertile_3 <- which(rev_cyt_tertiles[constru_tertile_3]==2)
rev_constru_tertile_1_df <- constru_tertile_1_df[-cyt_tertile_2_in_constru_tertile_1,]
rev_constru_tertile_2_df <- constru_tertile_2_df[-cyt_tertile_2_in_constru_tertile_2,]
rev_constru_tertile_3_df <- constru_tertile_3_df[-cyt_tertile_2_in_constru_tertile_3,]

#Make the plot for Constru tertile 1 
melted_constru_tertile_1_df <- as.data.frame(reshape2::melt(rev_constru_tertile_1_df))
melted_constru_tertile_1_df$category <- "Cyt Low"
unique_variables <- unique(melted_constru_tertile_1_df$variable)
for (variable in unique_variables)
{
  ids <- which(melted_constru_tertile_1_df$variable==variable)
  cyt_3_in_constru_1 <- setdiff(ids,ids[which(rev_cyt_tertiles[constru_tertile_1]==1)])
  melted_constru_tertile_1_df[cyt_3_in_constru_1,]$category <- "Cyt High"
}
infiltration_score_ids <- which(melted_constru_tertile_1_df$variable=="InfiltrationScore")
melted_constru_tertile_1_df <- melted_constru_tertile_1_df[-infiltration_score_ids,]
melted_constru_tertile_1_df$category <- factor(melted_constru_tertile_1_df$category,levels=c("Cyt Low","Cyt High"))

g_constru_tertile1 <- ggplot(data=melted_constru_tertile_1_df, aes(x=category, y= value))+ facet_wrap(~variable, nrow=3, ncol=9) +
  geom_boxplot(aes(fill=category))+ geom_point(aes(color=category),position="jitter",alpha=0.5) + xlab("Category")+ylab("Abundance")+
  theme_light()+ geom_signif(comparisons=list(c("Cyt High","Cyt Low")),map_signif_level = T)+
  theme(legend.text=element_text(size=16),
        legend.title=element_text(size=16))+ 
  scale_fill_manual(name="Category",values=c("Cyt Low"="blue","Cyt High"="red"))+
  scale_color_manual(name="Category",values=c("Cyt Low"="blue","Cyt High"="red"))+
  theme(axis.text.x = element_text(color = "grey20", size = 12, angle = 90, hjust = 0.95, vjust = .5, face = "plain"),
        strip.text = element_text(color = "black", size=12, angle=0, hjust = 0.5, vjust = 0.5, face = "plain"),
        strip.background=element_rect(colour="black",fill="white"),
        axis.text.y = element_text(color = "grey20", size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"),
        axis.title.x = element_text(color = "grey20", size = 16, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 16, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        title=element_text(color="grey20", size=16, face="plain"))
ggsave(filename="results/Training_ImmuneCellAI_Fractions_Cyt1_vs_Cyt3_Constru1.pdf",plot = g_constru_tertile1,
       device = pdf(), height = 12, width=15, units="in", dpi=300)
dev.off()

#Make the plot for cyt1 vs cyt3 for constru tertile 3
melted_constru_tertile_3_df <- as.data.frame(reshape2::melt(rev_constru_tertile_3_df))
melted_constru_tertile_3_df$category <- "Cyt Low"
unique_variables <- unique(melted_constru_tertile_3_df$variable)
for (variable in unique_variables)
{
  ids <- which(melted_constru_tertile_3_df$variable==variable)
  cyt_3_in_constru_3 <- setdiff(ids,ids[which(rev_cyt_tertiles[constru_tertile_3]==1)])
  melted_constru_tertile_3_df[cyt_3_in_constru_3,]$category <- "Cyt High"
}
infiltration_score_ids <- which(melted_constru_tertile_3_df$variable=="InfiltrationScore")
melted_constru_tertile_3_df <- melted_constru_tertile_3_df[-infiltration_score_ids,]
melted_constru_tertile_3_df$category <- factor(melted_constru_tertile_3_df$category,levels=c("Cyt Low","Cyt High"))

g_constru_tertile3 <- ggplot(data=melted_constru_tertile_3_df, aes(x=category, y= value))+ facet_wrap(~variable, nrow=3, ncol=9) +
  geom_boxplot(aes(fill=category))+ geom_point(aes(color=category),position="jitter",alpha=0.5) + xlab("Category")+ylab("Abundance")+
  theme_light()+ geom_signif(comparisons=list(c("Cyt High","Cyt Low")),map_signif_level = T)+
  theme(legend.text=element_text(size=16),
        legend.title=element_text(size=16))+ 
  scale_fill_manual(name="Category",values=c("Cyt Low"="blue","Cyt High"="red"))+
  scale_color_manual(name="Category",values=c("Cyt Low"="blue","Cyt High"="red"))+
  theme(axis.text.x = element_text(color = "grey20", size = 12, angle = 90, hjust = 0.95, vjust = .5, face = "plain"),
        strip.text = element_text(color = "black", size=12, angle=0, hjust = 0.5, vjust = 0.5, face = "plain"),
        strip.background=element_rect(colour="black",fill="white"),
        axis.text.y = element_text(color = "grey20", size = 12, angle = 0, hjust = 1, vjust = 0, face = "plain"),
        axis.title.x = element_text(color = "grey20", size = 16, angle = 0, hjust = .5, vjust = 0, face = "plain"),
        axis.title.y = element_text(color = "grey20", size = 16, angle = 90, hjust = .5, vjust = .5, face = "plain"),
        title=element_text(color="grey20", size=16, face="plain"))
ggsave(filename="results/Training_ImmuneCellAI_Fractions_Cyt1_vs_Cyt3_Constru3.pdf",plot = g_constru_tertile3,
       device = pdf(), height = 12, width=15, units="in", dpi=300)
dev.off()

