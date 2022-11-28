library(data.table)
library(ggplot2)
library(survival)
library(glmnet)
library(plotmo)
library(doMC)
library(dplyr)
library(Matrix)
library(preprocessCore)
library(survcomp)
library(survivalAnalysis)
library(biomaRt)
library(randomForestSRC)
library(Hmisc)
library(extrafont)
library(Rttf2pt1)
library(sva)
library(corrplot)
library(gplots)
loadfonts()
registerDoMC(cores=10)

setwd("~/Documents/Misc_Work/Other Work/Miller_Related/CONSTRU/")

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
test_pathway_activities <- gsva(as.matrix(final_test_all_expr_df), gset.idx.list = Selected.pathways, kcdf = "Gaussian", method = "gsva", BPPARAM=SerialParam(), parallel.sz=10)
write.table(train_pathway_activities, file="results/Train_Pathway_Activities.csv", row.names=T, col.names=T, quote=F, sep="\t")
write.table(test_pathway_activities, file="results/Test_Pathway_Activities.csv", row.names=T, col.names=T, quote=F, sep="\t")


#train_cyt_score <- c(new_gse82191_out[[2]], new_gse9891_out[[2]], new_gseov3_out[[2]])
#test_cyt_score <- c(new_gse140082_out[[2]], new_gse32062_out[[2]], new_gse53963_out[[2]])
train_cyt_score <- as.numeric(as.vector(colMeans(final_all_expr_df[c("GZMA","PRF1"),])))
test_cyt_score <- as.numeric(as.vector(colMeans(final_test_all_expr_df[c("GZMA","PRF1"),])))
train_constru_score <- c(new_gse82191_out[[9]],new_gse9891_out[[9]], new_gseov3_out[[9]])
test_constru_score <- c(new_gse140082_out[[9]],new_gse32062_out[[9]], new_gse53963_out[[9]])
train_cyt_tertiles <- c(new_gse82191_out[[2]],new_gse9891_out[[2]], new_gseov3_out[[2]])
test_cyt_tertiles <- c(new_gse140082_out[[2]],new_gse32062_out[[2]], new_gse53963_out[[2]])
train_constru_tertiles <- c(new_gse82191_out[[3]], new_gse9891_out[[3]], new_gseov3_out[[3]])
test_constru_tertiles <- c(new_gse140082_out[[3]], new_gse32062_out[[3]], new_gse53963_out[[3]])
train_os_event <- c(new_gse82191_out[[4]], new_gse9891_out[[4]], new_gseov3_out[[4]])
train_os_time <- c(new_gse82191_out[[5]], new_gse9891_out[[5]], new_gseov3_out[[5]])
test_os_event <- c(new_gse140082_out[[4]], new_gse32062_out[[4]], new_gse53963_out[[4]])
test_os_time <- c(new_gse140082_out[[5]], new_gse32062_out[[5]], new_gse53963_out[[5]])

#Make the training and testing data frame
train_df <- as.data.frame(cbind(data.frame(time=train_os_time, status=train_os_event, cyt_score = train_cyt_score, constru_score = train_constru_score, cyt_tertiles = train_cyt_tertiles, constru_tertiles = train_constru_tertiles), t(train_pathway_activities)))
test_df <- as.data.frame(cbind(data.frame(time=test_os_time, status=test_os_event, cyt_score = test_cyt_score, constru_score = test_constru_score, cyt_tertiles = test_cyt_tertiles, constru_tertiles = test_constru_tertiles), t(test_pathway_activities)))
save(train_df, test_df,file="Data/Training_Testing_OS.Rdata")

pathway_names <- colnames(train_df)[c(6:ncol(train_df))]
train_cox_proportional_df <- NULL
train_tertiles  <- list()
for (k in 1:length(pathway_names))
{
  df <- data.frame(Score=train_df[,pathway_names[k]])
  #df %>% mutate(tertiles = ntile(Score, 3)) %>% mutate(tertiles = if_else(tertiles == 1, 'Low', if_else(tertiles == 2, 'Medium', 'High'))) -> out_pathway_df
  tertile_cutoffs <- as.numeric(quantile(df$Score,probs = c(0.333,0.666)))
  high_group_ids <- which(df$Score>=tertile_cutoffs[2])
  medium_group_ids <- which(df$Score>tertile_cutoffs[1] & df$Score<tertile_cutoffs[2])
  low_group_ids <- which(df$Score<=tertile_cutoffs[1])
  train_tertiles[[k]] <- tertile_cutoffs
  #high_group_ids <- which(out_pathway_df$tertiles=="High")
  #medium_group_ids <- which(out_pathway_df$tertiles=="Medium")
  #low_group_ids <- which(out_pathway_df$tertiles=="Low")
  id_list <- list(high_group_ids, medium_group_ids, low_group_ids)
  names(id_list) <- c("High","Medium","Low")
  for (i in 1:length(id_list))
  {
    sample_ids <- id_list[[i]]
    sample_df <- train_df[sample_ids,c(1:4)]
    sample_cox <- coxph(Surv(time,status)~cyt_score, data = sample_df)
    temp_cox <- get_cox_info(sample_cox)
    temp_cox <- c(pathway_names[k],names(id_list)[i],temp_cox)
    train_cox_proportional_df <- rbind(train_cox_proportional_df, temp_cox)
  }
}
train_cox_proportional_df <- as.data.frame(train_cox_proportional_df)
colnames(train_cox_proportional_df) <- c("Pathway", "Tertile","beta", "mean", "low", "upper", "HR", "wald.test", "p.value")
rownames(train_cox_proportional_df) <- NULL
train_cox_proportional_df$beta <- as.numeric(as.vector(train_cox_proportional_df$beta))
train_cox_proportional_df$mean <- as.numeric(as.vector(train_cox_proportional_df$mean))
train_cox_proportional_df$low <- as.numeric(as.vector(train_cox_proportional_df$low))
train_cox_proportional_df$upper <- as.numeric(as.vector(train_cox_proportional_df$upper))
train_cox_proportional_df$wald.test <- as.numeric(as.vector(train_cox_proportional_df$wald.test))
train_cox_proportional_df$p.value <- as.numeric(as.vector(train_cox_proportional_df$p.value))
train_cox_proportional_df$p.adjust <- signif(p.adjust(train_cox_proportional_df$p.value, method="fdr"),digits=2)
write.table(x=train_cox_proportional_df,file="results/Training_Cox_Proportional_Pathway_Cytscore.csv",quote = F, row.names=F, col.names=T, sep="\t")
names(train_tertiles) <- pathway_names

################################################################################
test_cox_proportional_df <- NULL
for (k in 1:length(pathway_names))
{
  df <- data.frame(Score=test_df[,pathway_names[k]])
  #df %>% mutate(tertiles = ntile(Score, 3)) %>% mutate(tertiles = if_else(tertiles == 1, 'Low', if_else(tertiles == 2, 'Medium', 'High'))) -> out_pathway_df
  #high_group_ids <- which(out_pathway_df$tertiles=="High")
  #medium_group_ids <- which(out_pathway_df$tertiles=="Medium")
  #low_group_ids <- which(out_pathway_df$tertiles=="Low")
  tertile_cutoffs <- train_tertiles[[k]]
  high_group_ids <- which(df$Score>tertile_cutoffs[2])
  medium_group_ids <- which(df$Score>tertile_cutoffs[1] & df$Score<=tertile_cutoffs[2])
  low_group_ids <- which(df$Score<=tertile_cutoffs[1])
  
  id_list <- list(high_group_ids, medium_group_ids, low_group_ids)
  names(id_list) <- c("High","Medium","Low")
  for (i in 1:length(id_list))
  {
    sample_ids <- id_list[[i]]
    sample_df <- test_df[sample_ids,c(1:4)]
    sample_cox <- coxph(Surv(time,status)~cyt_score, data = sample_df)
    temp_cox <- get_cox_info(sample_cox)
    temp_cox <- c(pathway_names[k],names(id_list)[i],temp_cox)
    test_cox_proportional_df <- rbind(test_cox_proportional_df, temp_cox)
  }
}
test_cox_proportional_df <- as.data.frame(test_cox_proportional_df)
colnames(test_cox_proportional_df) <- c("Pathway","Tertile","beta", "mean", "low", "upper", "HR", "wald.test", "p.value")
rownames(test_cox_proportional_df) <- NULL
test_cox_proportional_df$beta <- as.numeric(as.vector(test_cox_proportional_df$beta))
test_cox_proportional_df$mean <- as.numeric(as.vector(test_cox_proportional_df$mean))
test_cox_proportional_df$low <- as.numeric(as.vector(test_cox_proportional_df$low))
test_cox_proportional_df$upper <- as.numeric(as.vector(test_cox_proportional_df$upper))
test_cox_proportional_df$wald.test <- as.numeric(as.vector(test_cox_proportional_df$wald.test))
test_cox_proportional_df$p.value <- as.numeric(as.vector(test_cox_proportional_df$p.value))
test_cox_proportional_df$p.adjust <- signif(p.adjust(test_cox_proportional_df$p.value, method="fdr"),digits=2)
write.table(x=test_cox_proportional_df,file="results/Testing_Cox_Proportional_Pathway_Cytscore.csv",quote = F, row.names=F, col.names=T, sep="\t")

#Make the plot based on pathway tertiles and hazard ratio comparison for CYT score in each tertile
################################################################################
#Make the hazards matrix
column_split_values <- c("Tertile low","Tertile int","Tertile high")
train_pathway_hazards_matrix <- as.matrix(cbind(as.numeric(train_cox_proportional_df[train_cox_proportional_df$Tertile=="Low",]$mean), 
                                                as.numeric(train_cox_proportional_df[train_cox_proportional_df$Tertile=="Medium",]$mean), 
                                                as.numeric(train_cox_proportional_df[train_cox_proportional_df$Tertile=="High",]$mean)))
colnames(train_pathway_hazards_matrix) <- column_split_values
rownames(train_pathway_hazards_matrix) <- pathway_names
train_pathway_pval_matrix <- as.matrix(cbind(as.numeric(train_cox_proportional_df[train_cox_proportional_df$Tertile=="Low",]$p.value), 
                                             as.numeric(train_cox_proportional_df[train_cox_proportional_df$Tertile=="Medium",]$p.value), 
                                             as.numeric(train_cox_proportional_df[train_cox_proportional_df$Tertile=="High",]$p.value)))
colnames(train_pathway_pval_matrix) <- column_split_values
rownames(train_pathway_pval_matrix) <- pathway_names

test_pathway_hazards_matrix <- as.matrix(cbind(as.numeric(test_cox_proportional_df[test_cox_proportional_df$Tertile=="Low",]$mean), 
                                                as.numeric(test_cox_proportional_df[test_cox_proportional_df$Tertile=="Medium",]$mean), 
                                                as.numeric(test_cox_proportional_df[test_cox_proportional_df$Tertile=="High",]$mean)))
colnames(test_pathway_hazards_matrix) <- column_split_values
rownames(test_pathway_hazards_matrix) <- pathway_names
test_pathway_pval_matrix <- as.matrix(cbind(as.numeric(test_cox_proportional_df[test_cox_proportional_df$Tertile=="Low",]$p.value), 
                                             as.numeric(test_cox_proportional_df[test_cox_proportional_df$Tertile=="Medium",]$p.value), 
                                             as.numeric(test_cox_proportional_df[test_cox_proportional_df$Tertile=="High",]$p.value)))
colnames(test_pathway_pval_matrix) <- column_split_values
rownames(test_pathway_pval_matrix) <- pathway_names

#Make the heatmap
column_split_values <- factor(column_split_values, levels = c("Tertile low","Tertile int","Tertile high"))
ha = HeatmapAnnotation(
  #empty = anno_empty(border = FALSE, height = unit(8, "mm")),
  Con = anno_block(gp = gpar(fill = c("brown","orange","#FEDD00"))), 
  show_legend=F)
col_fun1 <- colorRamp2(c(0.5,1,1.5),c("blue","white","red"))
ht_pathway_train = Heatmap(matrix=train_pathway_hazards_matrix, col_fun1, 
                       name = "Hazards Ratio", column_title = qq("Hazards for Cyt score across Pathway Tertiles in Training Set"),
                       cell_fun = function(j, i, x, y, width, height, fill) {
                         grid.rect(x = x, y = y, width = width, height = height, 
                                   gp = gpar(col = "grey", fill = NA))
                         if (train_pathway_pval_matrix[i,j]<0.05 & train_pathway_pval_matrix[i,j]>0.01)
                         {
                           grid.text(sprintf("*"), x, y, gp = gpar(fontsize = 11))
                         }
                         else if (train_pathway_pval_matrix[i,j]<0.01 & train_pathway_pval_matrix[i,j]>0.0001)
                         {
                           grid.text(sprintf("**"),x,y, gp=gpar(fontsize=11))
                         }
                         else if (train_pathway_pval_matrix[i,j]<1e-4)
                         {
                           grid.text(sprintf("***"),x,y, gp=gpar(fontsize=11))
                         }
                       },
                       width = unit(6, "cm"),
                       height = unit(15, "cm"),
                       cluster_columns = F,
                       cluster_rows = F,
                       row_names_centered = F,
                       show_row_names = T,
                       row_labels = rownames(train_pathway_hazards_matrix),
                       row_names_max_width = max_text_width(rownames(train_pathway_hazards_matrix)),
                       column_split = column_split_values,
                       cluster_column_slices=F,
                       show_column_names = T,
                       raster_quality = 2,
                       top_annotation = ha,
                       column_title_rot = 0,
                       column_title_side = "top",
                       column_dend_side = "top",
                       column_names_side = "top",
                       column_title_gp = gpar(fontsize=11, family="sans"),
                       row_names_gp = gpar(fontsize=9, family="sans"),
                       use_raster = T,
                       column_gap = unit(1, "mm"),
                       border_gp = gpar(col = "black", lty = 1),
                       border = T,
                       heatmap_legend_param = list(direction = "horizontal")
)
pdf("results/Training_Cox_Proportional_Pathway_Tertiles_CYT_Score.pdf",height = 10, width=6)
draw(ht_pathway_train, heatmap_legend_side="bottom")
dev.off()

#Make the figure for the test set
ht_pathway_test = Heatmap(matrix=test_pathway_hazards_matrix, col_fun1, 
                      name = "Hazards Ratio", column_title = qq("Hazards for Cyt score across Pathway Tertiles in Test Set"),
                      cell_fun = function(j, i, x, y, width, height, fill) {
                        grid.rect(x = x, y = y, width = width, height = height, 
                                  gp = gpar(col = "grey", fill = NA))
                        if (test_pathway_pval_matrix[i,j]<0.05 & test_pathway_pval_matrix[i,j]>0.01)
                        {
                          grid.text(sprintf("*"), x, y, gp = gpar(fontsize = 11))
                        }
                        else if (test_pathway_pval_matrix[i,j]<0.01 & test_pathway_pval_matrix[i,j]>0.0001)
                        {
                          grid.text(sprintf("**"),x,y, gp=gpar(fontsize=11))
                        }
                        else if (test_pathway_pval_matrix[i,j]<1e-4)
                        {
                          grid.text(sprintf("***"),x,y, gp=gpar(fontsize=11))
                        }
                      },
                      width = unit(6, "cm"),
                      height = unit(15, "cm"),
                      cluster_columns = F,
                      cluster_rows = F,
                      row_names_centered = F,
                      show_row_names = T,
                      row_labels = rownames(test_pathway_hazards_matrix),
                      row_names_max_width = max_text_width(rownames(test_pathway_hazards_matrix)),
                      column_split = column_split_values,
                      cluster_column_slices=F,
                      show_column_names = T,
                      raster_quality = 2,
                      top_annotation = ha,
                      column_title_rot = 0,
                      column_title_side = "top",
                      column_dend_side = "top",
                      column_names_side = "top",
                      column_title_gp = gpar(fontsize=11, family="sans"),
                      row_names_gp = gpar(fontsize=9, family="sans"),
                      use_raster = T,
                      column_gap = unit(1, "mm"),
                      border_gp = gpar(col = "black", lty = 1),
                      border = T,
                      heatmap_legend_param = list(direction = "horizontal")
)
pdf("results/Test_Cox_Proportional_Pathway_Tertiles_CYT_Score.pdf",height = 10, width=6)
draw(ht_pathway_test, heatmap_legend_side="bottom")
dev.off()


##################################################################################
correlation_matrix <- matrix(data=0, nrow=nrow(train_pathway_activities), ncol=2)
pval_matrix <- matrix(data=0, nrow=nrow(train_pathway_activities), ncol=2)
rownames(correlation_matrix) <- rownames(train_pathway_activities)
colnames(correlation_matrix) <- c("Training","Test")
rownames(pval_matrix) <- rownames(train_pathway_activities)
colnames(pval_matrix) <- c("Training","Test")
for (i in 1:nrow(train_pathway_activities))
{
  train_cor_info <- cor.test(train_pathway_activities[i,],train_cyt_score)
  train_cor_p_value <- train_cor_info$p.value
  train_cor_value <- as.numeric(train_cor_info$estimate)
  correlation_matrix[i,1] <- train_cor_value
  pval_matrix[i,1] <- train_cor_p_value
  
  test_cor_info <- cor.test(test_pathway_activities[i,],test_cyt_score)
  test_cor_p_value <- test_cor_info$p.value
  test_cor_value <- as.numeric(test_cor_info$estimate)
  correlation_matrix[i,2] <- test_cor_value
  pval_matrix[i,2] <- test_cor_p_value
}
pval_matrix[,1] <- p.adjust(pval_matrix[,1],method = "fdr")
pval_matrix[,2] <- p.adjust(pval_matrix[,2], method="fdr")

pdf("results/Correlation_Plot.pdf",height=16, width = 8)
corrplot(corr=correlation_matrix, 
         p.mat = pval_matrix,
         type="full", insig="pch", sig.level =.05, pch.cex = 0.5, col=bluered(100), cl.pos = "b", col.lim = c(-0.7,0.7))
dev.off()
