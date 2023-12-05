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
library(sva)
loadfonts()
registerDoMC(cores=20)

setwd("~/Documents/Misc_Work/Other Work/Miller_Related/CONSTRU/")

source("scripts/all_functions.R")

get_all_df_info <- function(df, mode)
{
  #Get sample names
  sample_ids <- as.character(as.vector(df[1,c(6:ncol(df))]))
  
  #get the CYT-Sig tertiles
  cyt_tertiles <- as.numeric(as.vector(df[2,c(6:ncol(df))]))
  
  #get the constru tertiles
  constru_tertiles <- as.numeric(as.vector(df[3,c(6:ncol(df))]))
  
  #Survival information
  os_event_information <- as.numeric(as.vector(df[4,c(6:ncol(df))]))
  os_time_information <- as.numeric(as.vector(df[5,c(6:ncol(df))]))
  
  #Get patient age
  age_information <- as.numeric(as.vector(df[6,c(6:ncol(df))]))
  
  #Get cancer stage
  cancer_stage <- as.character(as.vector(df[7,c(6:ncol(df))]))
  
  #Get cyt and constru avg score
  cyt_score <- as.numeric(as.vector(df[9,c(6:ncol(df))]))
  constru_score <- as.numeric(as.vector(df[10, c(6:ncol(df))]))
  
  #Get the metadata_info
  metadata_info <- as.character(as.vector(df[12,]))
  pearson_r_cyt_score <- as.numeric(as.vector(df[13:nrow(df),3]))
  pearson_r_constru_score <- as.numeric(as.vector(df[13:nrow(df),4]))
  
  #Get gene expression matrix and revise
  gene_expr_df <- df[c(13:nrow(df)),c(6:ncol(df))]
  for (i in 1:ncol(gene_expr_df))
  {
    gene_expr_df[,i] <- as.numeric(as.vector(gene_expr_df[,i]))
  }
  colnames(gene_expr_df) <- as.character(as.vector(metadata_info[c(6:length(metadata_info))]))
  
  #Get rid of genes with NA
  all_gene_names <- as.character(as.vector(df[13:nrow(df),2]))
  if (mode=="train")
  {
    all_gene_names <- unlist(lapply(strsplit(all_gene_names,split=" ///"),`[[`,1))
  }
  na_ids <- which(is.na(rowSums(gene_expr_df)))
  if (length(na_ids)>0)
  {
    rev_gene_expr_df <- gene_expr_df[-na_ids,]
    rev_all_gene_names <- all_gene_names[-na_ids]
  }else{
    rev_gene_expr_df <- gene_expr_df
    rev_all_gene_names <- all_gene_names
  }
  
  #Get rid of all genes with 0 or empty gene names
  other_ids <- which(rev_all_gene_names=="" | rev_all_gene_names=="0")
  if (length(other_ids)>0)
  {
    revised_gene_names <- rev_all_gene_names[-other_ids]
    rev_gene_expr_df <- rev_gene_expr_df[-other_ids,]
  }else{
    revised_gene_names <- rev_all_gene_names
    rev_gene_expr_df <- rev_gene_expr_df
  }
  
  unique_genes <- unique(revised_gene_names)
  final_expr_df <- NULL
  for (i in 1:length(unique_genes))
  {
    gene_name <- unique_genes[i]
    ids <- which(revised_gene_names==gene_name)
    if (length(ids)>1)
    {
      index <- as.numeric(which.max(rowSums(rev_gene_expr_df[ids,])))
      temp <- rev_gene_expr_df[ids[index],]
    }else{
      temp <- rev_gene_expr_df[ids,]
    }
    final_expr_df <- rbind(final_expr_df, temp)
  }
  final_expr_df <- as.data.frame(final_expr_df)
  rownames(final_expr_df) <- unique_genes
  
  return(list(sample_ids, cyt_tertiles, constru_tertiles, os_event_information, os_time_information, 
              age_information, cancer_stage, cyt_score, constru_score, metadata_info,
              pearson_r_cyt_score, pearson_r_constru_score, final_expr_df))
}

#Get the RNASeq data
gse82191_df <- fread("Data/GSE82191_OV431_22277_for Raghvendra_03-15-22.txt",header=F,sep="\t")
gse9891_df <- fread("Data/GSE9891_Tothill_OV227_for Raghvendra_03-30-22.txt",header=F,sep="\t")
gseov3_df <- fread("Data/OV3_221_22277_for Raghvendra_03-31-22.txt",header=F,sep="\t")
gse82191_df <- as.data.frame(gse82191_df)
gse9891_df <- as.data.frame(gse9891_df)
gseov3_df <- as.data.frame(gseov3_df)

gse82191_out <- get_all_df_info(gse82191_df, mode="train")
gse9891_out <- get_all_df_info(gse9891_df, mode="train")
gseov3_out <- get_all_df_info(gseov3_df, mode="train")

save(gse82191_out, file="Data/GSE82191_all_info.Rdata")
save(gse9891_out, file="Data/GSE9891_all_info.Rdata")
save(gseov3_out, file="Data/GSEOV3_all_info.Rdata")

#Load the final expression table
final_gse82191_expr_df <- gse82191_out[[13]]
final_gse9891_expr_df <- gse9891_out[[13]]
final_gseov3_expr_df <- gseov3_out[[13]]
common_genes <- intersect(intersect(rownames(final_gse82191_expr_df),rownames(final_gseov3_expr_df)),rownames(final_gse9891_expr_df))
all_expr_df <- as.data.frame(cbind(final_gse82191_expr_df[common_genes,], final_gse9891_expr_df[common_genes,], final_gseov3_expr_df[common_genes,]))
#Need to use combat in Training phase as expression distributions are very different
final_all_expr_df = ComBat(dat=all_expr_df, batch=c(rep("A",ncol(final_gse82191_expr_df)),rep("B",ncol(final_gse9891_expr_df)),rep("C",ncol(final_gseov3_expr_df))), 
                       mod=NULL, par.prior=TRUE, prior.plots=FALSE)
save(final_all_expr_df,file="Data/Training_Combined.Rdata")
mapping_df <- as.data.frame(cbind(c(rep("GSE82191",ncol(final_gse82191_expr_df)), rep("GSE9891",ncol(final_gse9891_expr_df)), rep("OV3",ncol(final_gseov3_expr_df))),
                                  c(colnames(final_gse82191_expr_df),colnames(final_gse9891_expr_df), colnames(final_gseov3_expr_df))))
colnames(mapping_df) <- c("Training_Data","Sample_Name")
write.table(mapping_df, "results/Sample_Name_Mapping_df.csv", row.names=F, col.names=T, sep=",", quote=F)

#Load the frequently dysregulated pathways in cancer
load("Data/Selected.pathways.3.4.RData")
pathway_activities <- gsva(as.matrix(final_all_expr_df), gset.idx.list = Selected.pathways, kcdf="Gaussian", method="gsva", BPPARAM=SerialParam(), parallel.sz=20)

#Load the gene lists for markers of immune cells
load("Data/immune.gene.lists.v3.Rdata")
immune_activities <- gsva(as.matrix(final_all_expr_df), gset.idx.list = Bindea_REV1, kcdf = "Gaussian", method="gsva", BPPARAM=SerialParam(), parallel.sz=20 )

#Perform immune deconvolution
immune_cell_type_conc <- immunedeconv::deconvolute(final_all_expr_df, 'xcell')
immune_cell_types <- immune_cell_type_conc$cell_type
immune_cell_type_concentrations <- as.data.frame(immune_cell_type_conc[,c(2:ncol(immune_cell_type_conc))])
rownames(immune_cell_type_concentrations) <- immune_cell_types
# conc_values <- as.numeric(as.vector(colSums(immune_cell_type_concentrations)))
# for (i in 1:length(conc_values))
# {
#   immune_cell_type_concentrations[,i] <- immune_cell_type_concentrations[,i]/conc_values[i]
# }

#Perform immune AI




#Build the complexheatmaps for the pathway activites, immune activities and immune deconvolutions
###############################################################################
col_fun1 <- colorRamp2(c(-1,0,1.0),c("blue","white","red"))
constru_tertiles <- c(gse82191_out[[3]],gse9891_out[[3]],gseov3_out[[3]])
cyt_tertiles <- c(gse82191_out[[2]], gse9891_out[[2]], gseov3_out[[2]])

#Order the ids by CONSTRU scores followed by Cyt scores
order_ids <- order(constru_tertiles, cyt_tertiles)
rev_constru_tertiles <- constru_tertiles[order_ids]
rev_cyt_tertiles <- cyt_tertiles[order_ids]
rev_pathway_activities <- pathway_activities[,order_ids]
rev_immune_cell_type_concentrations <- immune_cell_type_concentrations[,order_ids]
rev_immune_activities <- immune_activities[,order_ids]
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

column_split_values <- factor(column_split_values)
ha = HeatmapAnnotation(
  empty = anno_empty(border = FALSE, height = unit(8, "mm")),
  Cyt = anno_block(gp = gpar(fill = c("black","red","darkgreen","black","red","darkgreen","black","red","darkgreen")), 
                   labels = c("Cyt low","Cyt int","Cyt high","Cyt low","Cyt int","Cyt high","Cyt low","Cyt int","Cyt high")),
  simple_anno_size = unit(8, "mm"),
  annotation_name_gp = gpar(fontsize=11, family="sans",col="white"),
  annotation_legend_param = list(directon="horizontal"),
  show_legend=F)
ht_pathway = Heatmap(matrix=rev_pathway_activities, col_fun1, 
                  name = "Pathway Activities", column_title = qq("Pathway Activities across Constru Groups"),
                  width = unit(15, "cm"),
                  height = unit(15, "cm"),
                  cluster_columns = F,
                  cluster_rows = T,
                  row_names_centered = F,
                  show_row_names = T,
                  row_labels = rownames(rev_pathway_activities),
                  row_names_max_width = max_text_width(rownames(rev_pathway_activities)),
                  column_split = column_split_values,
                  cluster_column_slices=F,
                  show_column_names = F,
                  raster_quality = 2,
                  top_annotation = ha,
                  column_title_rot = 0,
                  column_title_side = "bottom",
                  column_dend_side = "top",
                  column_title_gp = gpar(fontsize=11, family="sans"),
                  row_names_gp = gpar(fontsize=9, family="sans"),
                  use_raster = T,
                  column_gap = unit(1, "mm"),
                  border_gp = gpar(col = "black", lty = 1),
                  border = T,
                  heatmap_legend_param = list(direction = "horizontal")
)
pdf("results/Training_Pathway_Activities_Constru_Cyt_Tertiles.pdf", height = 12, width=12, fonts = "sans")
draw(ht_pathway, heatmap_legend_side="bottom")
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
dev.off()

##Immune Concentrations
col_fun2 = colorRamp2(c(0,0.5),c("white","red"))
ht_immune_conc = Heatmap(matrix=as.matrix(rev_immune_cell_type_concentrations), col_fun2, 
                             name = "Immune Cell Fractions", column_title = qq("Immune Cell Fractions across Constru Groups"),
                             width = unit(15, "cm"),
                             height = unit(15, "cm"),
                             cluster_columns = F,
                             cluster_rows = T,
                             row_names_centered = F,
                             show_row_names = T,
                             row_labels = rownames(rev_immune_cell_type_concentrations),
                             row_names_max_width = max_text_width(rownames(rev_immune_cell_type_concentrations)),
                             column_split = column_split_values,
                             cluster_column_slices = T,
                             show_column_names = F,
                             raster_quality = 2,
                             top_annotation = ha,
                             column_title_rot = 0,
                             column_title_side = "bottom",
                             column_dend_side = "top",
                             column_title_gp = gpar(fontsize=11, family="sans"),
                             row_names_gp = gpar(fontsize=9, family="sans"),
                             use_raster = T,
                             column_gap = unit(1, "mm"),
                             border_gp = gpar(col = "black", lty = 1),
                             border = T,
                             heatmap_legend_param = list(direction = "horizontal"))
pdf("results/Training_Immune_Cell_Fractions_Constru_Cyt_Tertiles.pdf", height = 8, width=8, fonts = "sans")
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
dev.off()

#Perform additional analysis on pathways
################################################################################
#Melt the dataset to perform wilcox Training between the tertiles: Cyt 1 vs Cyt 3 for Constru tertiles 1,2,3
constru_tertiles_labels <- c(rep("Tertile 1",294),rep("Tertile 2",291),rep("Tertile 3",294))
mod_rev_pathway_activities <- as.data.frame(t(rbind(rbind(rev_pathway_activities,constru_tertiles_labels),rev_cyt_tertiles)))
mod_rev_pathway_activities_df <- reshape2::melt(mod_rev_pathway_activities,id.vars=c("constru_tertiles_labels","rev_cyt_tertiles"))
colnames(mod_rev_pathway_activities_df) <- c("Constru","Cyt","Pathway","Value")
mod_rev_pathway_activities_df$Cyt <- as.factor(as.vector(mod_rev_pathway_activities_df$Cyt))
mod_rev_pathway_activities_df$Constru <- as.factor(as.vector(mod_rev_pathway_activities_df$Constru))
mod_rev_pathway_activities_df$Value <- as.numeric(as.vector(mod_rev_pathway_activities_df$Value))
mod_rev_pathway_activities_df$Pathway <- as.character(as.vector(mod_rev_pathway_activities_df$Pathway))

average_pathway_activity_matrix <- matrix(0, nrow=nrow(rev_pathway_activities), ncol=9)
rownames(average_pathway_activity_matrix) <- rownames(rev_pathway_activities)
constru_high_vs_low_ratio_matrix <- matrix(0, nrow=nrow(rev_pathway_activities), ncol=3)
constru_high_vs_low_pval_matrix <- matrix(0, nrow=nrow(rev_pathway_activities), ncol=3)
cyt_constru1_high_vs_low_ratio_matrix <- matrix(0, nrow=nrow(rev_pathway_activities), ncol=3)
cyt_constru2_high_vs_low_ratio_matrix <- matrix(0, nrow=nrow(rev_pathway_activities), ncol=3)
cyt_constru3_high_vs_low_ratio_matrix <- matrix(0, nrow=nrow(rev_pathway_activities), ncol=3)
cyt_constru1_high_vs_low_pval_matrix <- matrix(0, nrow=nrow(rev_pathway_activities), ncol=3)
cyt_constru2_high_vs_low_pval_matrix <- matrix(0, nrow=nrow(rev_pathway_activities), ncol=3)
cyt_constru3_high_vs_low_pval_matrix <- matrix(0, nrow=nrow(rev_pathway_activities), ncol=3)
rownames(constru_high_vs_low_pval_matrix) <- rownames(rev_pathway_activities)
rownames(constru_high_vs_low_ratio_matrix) <- rownames(rev_pathway_activities)
rownames(cyt_constru1_high_vs_low_pval_matrix) <- rownames(rev_pathway_activities)
rownames(cyt_constru2_high_vs_low_pval_matrix) <- rownames(rev_pathway_activities)
rownames(cyt_constru3_high_vs_low_pval_matrix) <- rownames(rev_pathway_activities)
rownames(cyt_constru1_high_vs_low_ratio_matrix) <- rownames(rev_pathway_activities)
rownames(cyt_constru2_high_vs_low_ratio_matrix) <- rownames(rev_pathway_activities)
rownames(cyt_constru3_high_vs_low_ratio_matrix) <- rownames(rev_pathway_activities)
unique_pathways <- rownames(average_pathway_activity_matrix)
colnames(average_pathway_activity_matrix) <- NULL
unique_constru_tertiles <- unique(mod_rev_pathway_activities_df$Constru)
unique_cyt_tertiles <- unique(mod_rev_pathway_activities_df$Cyt)
for (k in 1:length(unique_pathways))
{
  pathway <- unique_pathways[k]
  fixed_pathway_df <- mod_rev_pathway_activities_df[mod_rev_pathway_activities_df$Pathway==pathway,]
  constru_high_activities <- fixed_pathway_df[fixed_pathway_df$Constru=="Tertile 3",]$Value
  constru_low_activities <- fixed_pathway_df[fixed_pathway_df$Constru=="Tertile 1",]$Value
  wx_test <- wilcox.test(constru_high_activities,constru_low_activities,paired=F)
  constru_high_vs_low_pval_matrix[pathway,2] <- wx_test$p.value
  constru_high_vs_low_ratio_matrix[pathway,2] <- mean(constru_high_activities)-mean(constru_low_activities)
  for (i in 1:length(unique_constru_tertiles))
  {
    constru_tertile <- unique_constru_tertiles[i]
    cyt_constru_high_activities <- fixed_pathway_df[fixed_pathway_df$Constru==constru_tertile & fixed_pathway_df$Cyt==3,]$Value
    cyt_constru_low_activities <- fixed_pathway_df[fixed_pathway_df$Constru==constru_tertile & fixed_pathway_df$Cyt==1,]$Value
    rev_wx_test <- wilcox.test(cyt_constru_high_activities,cyt_constru_low_activities,paired=F)

    if (constru_tertile=="Tertile 1")
    {
      cyt_constru1_high_vs_low_pval_matrix[pathway,2] <- rev_wx_test$p.value
      cyt_constru1_high_vs_low_ratio_matrix[pathway,2] <- mean(cyt_constru_high_activities)-mean(cyt_constru_low_activities)
    }else if (constru_tertile =="Tertile 2")
    {
      cyt_constru2_high_vs_low_pval_matrix[pathway,2] <- rev_wx_test$p.value
      cyt_constru2_high_vs_low_ratio_matrix[pathway,2] <- mean(cyt_constru_high_activities)-mean(cyt_constru_low_activities)
    }else{
      cyt_constru3_high_vs_low_pval_matrix[pathway,2] <- rev_wx_test$p.value
      cyt_constru3_high_vs_low_ratio_matrix[pathway,2] <- mean(cyt_constru_high_activities)-mean(cyt_constru_low_activities)
    }
    
    for (j in 1:length(unique_cyt_tertiles))
    {
      cyt_tertile <- unique_cyt_tertiles[j]
      average_pathway_activity_matrix[pathway,((i-1)*3+j)] <- mean(mod_rev_pathway_activities_df[mod_rev_pathway_activities_df$Pathway==pathway & mod_rev_pathway_activities_df$Constru==constru_tertile & mod_rev_pathway_activities_df$Cyt==cyt_tertile,]$Value)
    }
  }
}
constru_high_vs_low_pval_matrix[,2] <- p.adjust(constru_high_vs_low_pval_matrix[,2], method="fdr")
cyt_constru1_high_vs_low_pval_matrix[,2] <- p.adjust(cyt_constru1_high_vs_low_pval_matrix[,2], method="fdr")
cyt_constru2_high_vs_low_pval_matrix[,2] <- p.adjust(cyt_constru2_high_vs_low_pval_matrix[,2], method="fdr")
cyt_constru3_high_vs_low_pval_matrix[,2] <- p.adjust(cyt_constru3_high_vs_low_pval_matrix[,2], method="fdr")

#Make the average heatmap figure for pathway activities 
#############################################################################################################################################################
col_fun3 <- colorRamp2(c(-0.5,0,0.5),c("blue","white","red"))
ha = HeatmapAnnotation(
  empty = anno_empty(border = FALSE, height = unit(8, "mm")),
  Constru = anno_block(gp=gpar(fill=c("brown","orange","#FEDD00"),col="white"),
                       labels=c("Con low","Con int","Con high")),
  Cyt = c("black","red","darkgreen","black","red","darkgreen","black","red","darkgreen"),
  col=list(Cyt=c("black"="black","red"="red","darkgreen"="darkgreen")),
  simple_anno_size = unit(8, "mm"),
  annotation_name_gp = gpar(fontsize=11,col="white"),
  annotation_legend_param = list(directon="horizontal"),
  show_legend=F)
ht_rev_pathway = Heatmap(matrix=average_pathway_activity_matrix, col_fun3, 
                     name = "Pathway Activity", column_title = qq("Average Pathway Activity"),
                     width = unit(6, "cm"),
                     height = unit(15, "cm"),
                     cluster_columns = F,
                     cluster_rows = F,
                     row_names_centered = F,
                     show_row_names = T,
                     row_labels = rownames(average_pathway_activity_matrix),
                     row_names_max_width = max_text_width(rownames(average_pathway_activity_matrix)),
                     show_column_names = F,
                     column_split = c(1,1,1,2,2,2,3,3,3),
                     raster_quality = 2,
                     top_annotation = ha,
                     column_title_rot = 0,
                     column_title_side = "top",
                     column_dend_side = "bottom",
                     row_names_side = "left",
                     column_title_gp = gpar(fontsize=11, family="sans"),
                     row_names_gp = gpar(fontsize=9, family="sans"),
                     use_raster = T,
                     column_gap = unit(1, "mm"),
                     border_gp = gpar(col = "black", lty = 1),
                     border = F,
                     heatmap_legend_param = list(direction = "horizontal")
)

ha2 = HeatmapAnnotation(
  empty = anno_empty(border = FALSE, height = unit(8, "mm")),
  Constru = anno_block(gp=gpar(fill=c("brown","white","#FEDD00"),col="white"),
                       labels=c("low","vs","high")),
  Cyt = c("black","red","darkgreen"),
  col=list(Cyt=c("black"="white","red"="white","darkgreen"="white")),
  simple_anno_size = unit(8, "mm"),
  annotation_name_gp = gpar(fontsize=11,col="white"),
  annotation_legend_param = list(directon="horizontal"),
  show_legend=F)
ht_constru_mean_pathway = Heatmap(matrix=constru_high_vs_low_ratio_matrix, col_fun3, 
                         name="Diff Activity",
                         width = unit(3, "cm"),
                         height = unit(15, "cm"),
                         cell_fun = function(j, i, x, y, w, h, col) {
                           if (constru_high_vs_low_pval_matrix[i,j]>0.05)
                           {
                             r = min(unit.c(w, h))*0.0
                             grid.circle(x,y,r,gp=gpar(fill="white"))
                           }
                           else if(constru_high_vs_low_pval_matrix[i,j]<0.05 & constru_high_vs_low_pval_matrix[i,j]>0.01 & constru_high_vs_low_pval_matrix[i,j]>0) {
                             r = min(unit.c(w, h))*2
                             grid.circle(x, y, r, gp = gpar(fill = "red"))
                           } 
                           else if(constru_high_vs_low_pval_matrix[i,j]<0.01 & constru_high_vs_low_pval_matrix[i,j]>0.001 & constru_high_vs_low_pval_matrix[i,j]>0) {
                             r = min(unit.c(w, h))*4
                             grid.circle(x, y, r, gp = gpar(fill = "red"))
                           }
                           else if(constru_high_vs_low_pval_matrix[i,j]<0.001 & constru_high_vs_low_pval_matrix[i,j]>0) {
                             r = min(unit.c(w, h))*6
                             grid.circle(x, y, r, gp = gpar(fill = "red"))
                           }
                         },
                         cluster_columns = F,
                         cluster_rows = F,
                         row_names_centered = F,
                         show_row_names = F,
                         row_labels = rownames(constru_high_vs_low_ratio_matrix),
                         row_names_max_width = max_text_width(rownames(constru_high_vs_low_ratio_matrix)),
                         show_column_names = F,
                         column_split = c(1,2,3),
                         raster_quality = 2,
                         top_annotation = ha2,
                         column_title = "",
                         column_title_rot = 0,
                         column_title_side = "top",
                         column_dend_side = "bottom",
                         row_names_side = "left",
                         column_title_gp = gpar(fontsize=11, family="sans"),
                         row_names_gp = gpar(fontsize=9, family="sans"),
                         use_raster = T,
                         column_gap = unit(1, "mm"),
                         border_gp = gpar(col = "black", lty = 1),
                         border = F,
                         heatmap_legend_param = list(direction = "horizontal")
)

ht_constru_cyt1_pathway <- make_bubble_plot_heatmap(activity_matrix = cyt_constru1_high_vs_low_ratio_matrix, cyt_constru1_high_vs_low_pval_matrix, tertile=1, col_fun3)
ht_constru_cyt2_pathway <- make_bubble_plot_heatmap(activity_matrix = cyt_constru2_high_vs_low_ratio_matrix, cyt_constru2_high_vs_low_pval_matrix, tertile=2, col_fun3)
ht_constru_cyt3_pathway <- make_bubble_plot_heatmap(activity_matrix = cyt_constru3_high_vs_low_ratio_matrix, cyt_constru3_high_vs_low_pval_matrix, tertile=3, col_fun3)

ht_list <- ht_rev_pathway+ht_constru_mean_pathway+ht_constru_cyt1_pathway+ht_constru_cyt2_pathway+ht_constru_cyt3_pathway
pdf(file="results/Training_Pathway_Activities_Combined.pdf", height = 10, width = 10, fonts="sans")
draw(ht_list, heatmap_legend_side="bottom")
dev.off()

##Get the mean and se information for pathway activites
#out_pathways <- perform_stat_Training(input_df = rev_pathway_activities_df)
#stat.Training.pathways <- out_pathways[[1]]
#stat_table.pathways <- out_pathways[[2]]
#g_pathways <- make_mean_se_stat_plot(rev_pathway_activities_df, stat.Training.pathways, plot_title = "Pathway Activity for Constru Tertiles", xlab_label = "Pathways", ylab_label = "Mean Activity")
#ggsave(filename="results/Pathway_Activities_Mean_Constru_Cyt.pdf",plot = g_pathways, device = pdf(), height = 10, width=10, units="in", dpi=300)
#dev.off()
#write.table(stat.Training.pathways, file="results/Pathway_Activities_Stats.csv",row.names=F, col.names=T, quote=F)

#Perform additional analysis on immune cell type concentrations
################################################################################
#Melt the dataset to perform wilcox Training between the tertiles: Cyt 1 vs Cyt 3 for Constru tertiles 1,2,3
mod_rev_immune_conc <- as.data.frame(t(rbind(rbind(rev_immune_cell_type_concentrations,constru_tertiles_labels),rev_cyt_tertiles)))
colnames(mod_rev_immune_conc)[c(40:41)] <- c("constru_tertiles_labels","cyt_tertiles")
mod_rev_immune_conc_df <- reshape2::melt(mod_rev_immune_conc,id.vars=c("constru_tertiles_labels","cyt_tertiles"))
colnames(mod_rev_immune_conc_df) <- c("Constru","Cyt","Immune_Celltype","Value")
mod_rev_immune_conc_df$Cyt <- as.factor(as.vector(mod_rev_immune_conc_df$Cyt))
mod_rev_immune_conc_df$Constru <- as.factor(as.vector(mod_rev_immune_conc_df$Constru))
mod_rev_immune_conc_df$Value <- as.numeric(as.vector(mod_rev_immune_conc_df$Value))
mod_rev_immune_conc_df$Immune_Celltype <- as.character(as.vector(mod_rev_immune_conc_df$Immune_Celltype))


average_immune_conc_matrix <- matrix(0, nrow=nrow(rev_immune_cell_type_concentrations), ncol=9)
rownames(average_immune_conc_matrix) <- rownames(rev_immune_cell_type_concentrations)
immune_constru_high_vs_low_ratio_matrix <- matrix(1, nrow=nrow(rev_immune_cell_type_concentrations), ncol=3)
immune_constru_high_vs_low_pval_matrix <- matrix(1, nrow=nrow(rev_immune_cell_type_concentrations), ncol=3)
immune_cyt_constru1_high_vs_low_ratio_matrix <- matrix(1, nrow=nrow(rev_immune_cell_type_concentrations), ncol=3)
immune_cyt_constru2_high_vs_low_ratio_matrix <- matrix(1, nrow=nrow(rev_immune_cell_type_concentrations), ncol=3)
immune_cyt_constru3_high_vs_low_ratio_matrix <- matrix(1, nrow=nrow(rev_immune_cell_type_concentrations), ncol=3)
immune_cyt_constru1_high_vs_low_pval_matrix <- matrix(1, nrow=nrow(rev_immune_cell_type_concentrations), ncol=3)
immune_cyt_constru2_high_vs_low_pval_matrix <- matrix(1, nrow=nrow(rev_immune_cell_type_concentrations), ncol=3)
immune_cyt_constru3_high_vs_low_pval_matrix <- matrix(1, nrow=nrow(rev_immune_cell_type_concentrations), ncol=3)
rownames(immune_constru_high_vs_low_pval_matrix) <- rownames(rev_immune_cell_type_concentrations)
rownames(immune_constru_high_vs_low_ratio_matrix) <- rownames(rev_immune_cell_type_concentrations)
rownames(immune_cyt_constru1_high_vs_low_pval_matrix) <- rownames(rev_immune_cell_type_concentrations)
rownames(immune_cyt_constru2_high_vs_low_pval_matrix) <- rownames(rev_immune_cell_type_concentrations)
rownames(immune_cyt_constru3_high_vs_low_pval_matrix) <- rownames(rev_immune_cell_type_concentrations)
rownames(immune_cyt_constru1_high_vs_low_ratio_matrix) <- rownames(rev_immune_cell_type_concentrations)
rownames(immune_cyt_constru2_high_vs_low_ratio_matrix) <- rownames(rev_immune_cell_type_concentrations)
rownames(immune_cyt_constru3_high_vs_low_ratio_matrix) <- rownames(rev_immune_cell_type_concentrations)
unique_immune_conc <- rownames(average_immune_conc_matrix)
colnames(average_immune_conc_matrix) <- NULL

unique_constru_tertiles <- unique(mod_rev_immune_conc_df$Constru)
unique_cyt_tertiles <- unique(mod_rev_immune_conc_df$Cyt)
for (k in 1:length(unique_immune_conc))
{
  immune_type <- unique_immune_conc[k]
  fixed_immune_conc_df <- mod_rev_immune_conc_df[mod_rev_immune_conc_df$Immune_Celltype==immune_type,]
  constru_high_conc <- fixed_immune_conc_df[fixed_immune_conc_df$Constru=="Tertile 3",]$Value
  constru_low_conc <- fixed_immune_conc_df[fixed_immune_conc_df$Constru=="Tertile 1",]$Value
  wx_test <- wilcox.test(constru_high_conc,constru_low_conc,paired=F)
  immune_constru_high_vs_low_pval_matrix[immune_type,2] <- wx_test$p.value
  immune_constru_high_vs_low_ratio_matrix[immune_type,2] <- mean(constru_high_conc)/mean(constru_low_conc)
  for (i in 1:length(unique_constru_tertiles))
  {
    constru_tertile <- unique_constru_tertiles[i]
    cyt_constru_high_conc <- fixed_immune_conc_df[fixed_immune_conc_df$Constru==constru_tertile & fixed_immune_conc_df$Cyt==3,]$Value
    cyt_constru_low_conc <- fixed_immune_conc_df[fixed_immune_conc_df$Constru==constru_tertile & fixed_immune_conc_df$Cyt==1,]$Value
    rev_wx_test <- wilcox.test(cyt_constru_high_conc,cyt_constru_low_conc,paired=F,exact=F)
    
    if (constru_tertile=="Tertile 1")
    {
      immune_cyt_constru1_high_vs_low_pval_matrix[immune_type,2] <- rev_wx_test$p.value
      immune_cyt_constru1_high_vs_low_ratio_matrix[immune_type,2] <- mean(cyt_constru_high_conc)/mean(cyt_constru_low_conc)
    }else if (constru_tertile =="Tertile 2")
    {
      immune_cyt_constru2_high_vs_low_pval_matrix[immune_type,2] <- rev_wx_test$p.value
      immune_cyt_constru2_high_vs_low_ratio_matrix[immune_type,2] <- mean(cyt_constru_high_conc)/mean(cyt_constru_low_conc)
    }else{
      immune_cyt_constru3_high_vs_low_pval_matrix[immune_type,2] <- rev_wx_test$p.value
      immune_cyt_constru3_high_vs_low_ratio_matrix[immune_type,2] <- mean(cyt_constru_high_conc)/mean(cyt_constru_low_conc)
    }
    
    for (j in 1:length(unique_cyt_tertiles))
    {
      cyt_tertile <- unique_cyt_tertiles[j]
      average_immune_conc_matrix[immune_type,((i-1)*3+j)] <- mean(mod_rev_immune_conc_df[mod_rev_immune_conc_df$Immune_Celltype==immune_type & 
                                                                                     mod_rev_immune_conc_df$Constru==constru_tertile & 
                                                                                     mod_rev_immune_conc_df$Cyt==cyt_tertile,]$Value)
    }
  }
}
immune_constru_high_vs_low_pval_matrix[,2] <- p.adjust(immune_constru_high_vs_low_pval_matrix[,2], method="fdr")
immune_cyt_constru1_high_vs_low_pval_matrix[,2] <- p.adjust(immune_cyt_constru1_high_vs_low_pval_matrix[,2], method="fdr")
immune_cyt_constru2_high_vs_low_pval_matrix[,2] <- p.adjust(immune_cyt_constru2_high_vs_low_pval_matrix[,2], method="fdr")
immune_cyt_constru3_high_vs_low_pval_matrix[,2] <- p.adjust(immune_cyt_constru3_high_vs_low_pval_matrix[,2], method="fdr")

#Make the average heatmap figure for pathway activities 
#############################################################################################################################################################
col_fun3 <- colorRamp2(c(0,0.01,0.1),c("blue","white","red"))
ha = HeatmapAnnotation(
  empty = anno_empty(border = FALSE, height = unit(8, "mm")),
  Constru = anno_block(gp=gpar(fill=c("brown","orange","#FEDD00"),col="white"),
                       labels=c("Con low","Con int","Con high")),
  Cyt = c("black","red","darkgreen","black","red","darkgreen","black","red","darkgreen"),
  col=list(Cyt=c("black"="black","red"="red","darkgreen"="darkgreen")),
  simple_anno_size = unit(8, "mm"),
  annotation_name_gp = gpar(fontsize=11,col="white"),
  annotation_legend_param = list(directon="horizontal"),
  show_legend=F)
ht_rev_immune_conc = Heatmap(matrix=average_immune_conc_matrix, col_fun3, 
                         name = "Immune Celltype Fraction", column_title = qq("Average Immune Celltype Fraction"),
                         width = unit(6, "cm"),
                         height = unit(15, "cm"),
                         cluster_columns = F,
                         cluster_rows = F,
                         row_names_centered = F,
                         show_row_names = T,
                         row_labels = rownames(average_immune_conc_matrix),
                         row_names_max_width = max_text_width(rownames(average_immune_conc_matrix)),
                         show_column_names = F,
                         column_split = c(1,1,1,2,2,2,3,3,3),
                         raster_quality = 2,
                         top_annotation = ha,
                         column_title_rot = 0,
                         column_title_side = "top",
                         column_dend_side = "bottom",
                         row_names_side = "left",
                         column_title_gp = gpar(fontsize=11, family="sans"),
                         row_names_gp = gpar(fontsize=9, family="sans"),
                         use_raster = T,
                         column_gap = unit(1, "mm"),
                         border_gp = gpar(col = "black", lty = 1),
                         border = F,
                         heatmap_legend_param = list(direction = "horizontal")
)

ha2 = HeatmapAnnotation(
  empty = anno_empty(border = FALSE, height = unit(8, "mm")),
  Constru = anno_block(gp=gpar(fill=c("brown","white","#FEDD00"),col="white"),
                       labels=c("low","vs","high")),
  Cyt = c("black","red","darkgreen"),
  col=list(Cyt=c("black"="white","red"="white","darkgreen"="white")),
  simple_anno_size = unit(8, "mm"),
  annotation_name_gp = gpar(fontsize=11,col="white"),
  annotation_legend_param = list(directon="horizontal"),
  show_legend=F)
col_fun4 = colorRamp2(c(0,1,2),c("blue","white","red"))
ht_constru_mean_immune_conc = Heatmap(matrix=immune_constru_high_vs_low_ratio_matrix, col_fun4, 
                                  name="Ratio",
                                  width = unit(3, "cm"),
                                  height = unit(15, "cm"),
                                  cell_fun = function(j, i, x, y, w, h, col) {
                                    if (immune_constru_high_vs_low_pval_matrix[i,j]>0.05)
                                    {
                                      r = min(unit.c(w, h))*0.0
                                      grid.circle(x,y,r,gp=gpar(fill="white"))
                                    }
                                    else if(immune_constru_high_vs_low_pval_matrix[i,j]<0.05 & immune_constru_high_vs_low_pval_matrix[i,j]>0.01 & immune_constru_high_vs_low_pval_matrix[i,j]>0) {
                                      r = min(unit.c(w, h))*2
                                      grid.circle(x, y, r, gp = gpar(fill = "red"))
                                    } 
                                    else if(immune_constru_high_vs_low_pval_matrix[i,j]<0.01 & immune_constru_high_vs_low_pval_matrix[i,j]>0.001 & immune_constru_high_vs_low_pval_matrix[i,j]>0) {
                                      r = min(unit.c(w, h))*4
                                      grid.circle(x, y, r, gp = gpar(fill = "red"))
                                    }
                                    else if(immune_constru_high_vs_low_pval_matrix[i,j]<0.001 & immune_constru_high_vs_low_pval_matrix[i,j]>0) {
                                      r = min(unit.c(w, h))*6
                                      grid.circle(x, y, r, gp = gpar(fill = "red"))
                                    }
                                  },
                                  cluster_columns = F,
                                  cluster_rows = F,
                                  row_names_centered = F,
                                  show_row_names = F,
                                  row_labels = rownames(immune_constru_high_vs_low_ratio_matrix),
                                  row_names_max_width = max_text_width(rownames(immune_constru_high_vs_low_ratio_matrix)),
                                  show_column_names = F,
                                  column_split = c(1,2,3),
                                  raster_quality = 2,
                                  top_annotation = ha2,
                                  column_title = "",
                                  column_title_rot = 0,
                                  column_title_side = "top",
                                  column_dend_side = "bottom",
                                  row_names_side = "left",
                                  column_title_gp = gpar(fontsize=11, family="sans"),
                                  row_names_gp = gpar(fontsize=9, family="sans"),
                                  use_raster = T,
                                  column_gap = unit(1, "mm"),
                                  border_gp = gpar(col = "black", lty = 1),
                                  border = F,
                                  heatmap_legend_param = list(direction = "horizontal")
)

ht_constru_cyt1_immune_conc <- make_bubble_plot_heatmap(activity_matrix = immune_cyt_constru1_high_vs_low_ratio_matrix, immune_cyt_constru1_high_vs_low_pval_matrix, tertile=1, col_fun4)
ht_constru_cyt2_immune_conc <- make_bubble_plot_heatmap(activity_matrix = immune_cyt_constru2_high_vs_low_ratio_matrix, immune_cyt_constru2_high_vs_low_pval_matrix, tertile=2, col_fun4)
ht_constru_cyt3_immune_conc <- make_bubble_plot_heatmap(activity_matrix = immune_cyt_constru3_high_vs_low_ratio_matrix, immune_cyt_constru3_high_vs_low_pval_matrix, tertile=3, col_fun4)

ht_list <- ht_rev_immune_conc+ht_constru_mean_immune_conc+ht_constru_cyt1_immune_conc+ht_constru_cyt2_immune_conc+ht_constru_cyt3_immune_conc
pdf(file="results/Training_Immune_Conc_Combined.pdf", height = 10, width = 10, fonts="sans")
draw(ht_list, heatmap_legend_side="bottom")
dev.off()



##Get the mean and se information for pathway activites
#out_immune_conc <- perform_stat_Training(input_df = rev_immune_conc_df)
#stat.Training.immune_conc <- out_immune_conc[[1]]
#g_immune_conc <- make_mean_se_stat_plot(rev_immune_conc_df, stat.Training.immune_conc, plot_title = "Immune Fractions for Constru Tertiles", xlab_label = "Immune Celltypes", ylab_label = "Mean Fractions")
#ggsave(filename="results/Immune_Conc_Mean_Constru_Cyt.pdf",plot = g_immune_conc, device = pdf(), height = 10, width=10, units="in", dpi=300)
#dev.off()
#write.table(stat.Training.immune_conc, file="results/Immune_Conc_Stats.csv",row.names=F, col.names=T, quote=F)


