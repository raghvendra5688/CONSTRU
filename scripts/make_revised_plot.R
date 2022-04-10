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

setwd("~/QCRI_PostDoc/Raghav_Related/Lance_Miller_Related/CONSTRU/")

source("scripts/all_functions.R")

#Training data
load("Data/GSE82191_all_info.Rdata")
load("Data/GSE9891_all_info.Rdata")
load("Data/GSEOV3_all_info.Rdata")
new_gse82191_out <- gse82191_out
new_gse9891_out <- gse9891_out
new_gseov3_out <- gseov3_out

#Test data
load("Data/GSE140082_all_info.Rdata")
load("Data/GSE32062_all_info.Rdata")
load("Data/GSE53963_all_info.Rdata")
new_gse140082_out <- gse82191_out
new_gse32062_out <- gse9891_out
new_gse53963_out <- gseov3_out

#Load the final expression table
final_gse82191_expr_df <- new_gse82191_out[[13]]
final_gse9891_expr_df <- new_gse9891_out[[13]]
final_gseov3_expr_df <- new_gseov3_out[[13]]
common_genes <- intersect(intersect(rownames(final_gse82191_expr_df),rownames(final_gseov3_expr_df)),rownames(final_gse9891_expr_df))
all_expr_df <- as.data.frame(cbind(final_gse82191_expr_df[common_genes,], final_gse9891_expr_df[common_genes,], final_gseov3_expr_df[common_genes,]))

#Need to use combat in Training phase as expression distributions are very different
final_all_expr_df = ComBat(dat=all_expr_df, batch=c(rep("A",ncol(final_gse82191_expr_df)),rep("B",ncol(final_gse9891_expr_df)),rep("C",ncol(final_gseov3_expr_df))), 
                           mod=NULL, par.prior=TRUE, prior.plots=FALSE)

#Load the test expression table
final_gse140082_expr_df <- new_gse140082_out[[13]]
final_gse32062_expr_df <- new_gse32062_out[[13]]
final_gse53963_expr_df <- new_gse53963_out[[13]]
common_genes <- intersect(intersect(rownames(final_gse140082_expr_df),rownames(final_gse32062_expr_df)),rownames(final_gse53963_expr_df))
test_all_expr_df <- as.data.frame(cbind(final_gse140082_expr_df[common_genes,], final_gse32062_expr_df[common_genes,], final_gse53963_expr_df[common_genes,]))

#Need to use combat in Training phase as expression distributions are very different
final_test_all_expr_df = ComBat(dat=test_all_expr_df, batch=c(rep("A",ncol(final_gse140082_expr_df)),rep("B",ncol(final_gse32062_expr_df)),rep("C",ncol(final_gse53963_expr_df))), 
                           mod=NULL, par.prior=TRUE, prior.plots=FALSE)


#Load the frequently dysregulated pathways in cancer
load("Data/Selected.pathways.3.4.RData")
train_pathway_activities <- gsva(as.matrix(final_all_expr_df), gset.idx.list = Selected.pathways, kcdf="Gaussian", method="gsva", BPPARAM=SerialParam(), parallel.sz=10)

#Load the gene lists for markers of immune cells
load("Data/immune.gene.lists.v3.Rdata")
train_immune_activities <- gsva(as.matrix(final_all_expr_df), gset.idx.list = Bindea_REV1, kcdf = "Gaussian", method="gsva", BPPARAM=SerialParam(), parallel.sz=10 )

#Perform immune deconvolution
train_immune_cell_type_conc <- immunedeconv::deconvolute(final_all_expr_df, 'xcell')
train_immune_cell_types <- train_immune_cell_type_conc$cell_type
train_immune_cell_type_concentrations <- as.data.frame(train_immune_cell_type_conc[,c(2:ncol(train_immune_cell_type_conc))])
rownames(train_immune_cell_type_concentrations) <- train_immune_cell_types

test_pathway_activities <- gsva(as.matrix(final_test_all_expr_df), gset.idx.list = Selected.pathways, kcdf="Gaussian", method="gsva", BPPARAM=SerialParam(), parallel.sz=10)
test_immune_activities <-  gsva(as.matrix(final_test_all_expr_df), gset.idx.list = Bindea_REV1, kcdf = "Gaussian", method="gsva", BPPARAM=SerialParam(), parallel.sz=10 )
test_immune_cell_type_conc <- immunedeconv::deconvolute(final_test_all_expr_df, 'xcell')
test_immune_cell_types <- test_immune_cell_type_conc$cell_type
test_immune_cell_type_concentrations <- as.data.frame(test_immune_cell_type_conc[,c(2:ncol(test_immune_cell_type_conc))])
rownames(test_immune_cell_type_concentrations) <- test_immune_cell_types

train_out <- make_pathway_activity_heatmap(new_gse82191_out, new_gse9891_out, new_gseov3_out, final_all_expr_df = final_all_expr_df,
                              pathway_activities = train_pathway_activities, immune_cell_type_concentrations = train_immune_cell_type_concentrations,
                              immune_activities = train_immune_activities)
test_out <- make_pathway_activity_heatmap(new_gse140082_out, new_gse32062_out, new_gse53963_out, final_all_expr_df = final_test_all_expr_df,
                                                 pathway_activities = test_pathway_activities, immune_cell_type_concentrations = test_immune_cell_type_concentrations,
                                                 immune_activities = test_immune_activities)
ht_list <- train_out[[1]]+test_out[[1]]
pdf("results/Training_Test_Pathway_Activities_Constru_Cyt_Tertiles.pdf", height = 12, width=20, fonts = "sans")
draw(ht_list, heatmap_legend_side="bottom")
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
rev_train_immune_cell_type_concentrations <- train_out[[2]]
rev_test_immune_cell_type_concentrations <- test_out[[2]]
train_column_split_values <- train_out[[3]]
test_column_split_values <- test_out[[3]]
ha <- train_out[[4]]
train_ht_immune_conc <- make_immune_conc_heatmap(rev_immune_cell_type_concentrations = rev_train_immune_cell_type_concentrations, train_column_split_values, ha=ha)
test_ht_immune_conc <- make_immune_conc_heatmap(rev_immune_cell_type_concentrations = rev_test_immune_cell_type_concentrations, test_column_split_values, ha)

ht_immune_conc <- train_ht_immune_conc+test_ht_immune_conc
pdf("results/Training_Test_Immune_Cell_Fractions_Constru_Cyt_Tertiles.pdf", height = 8, width=14, fonts = "sans")
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
