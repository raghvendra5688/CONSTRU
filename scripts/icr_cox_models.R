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

load("Data/Training_Testing_OS.Rdata")
icr_genes <- c("CXCL9","CXCL10","CCL5","IFNG","IL12B","TBX21", "CD8A", "CD8B", "STAT1", "IRF1","GNLY", "PRF1", "GZMA", "GZMB", "GZMH", "CD274", "CTLA4", "FOXP3", "IDO1", "PDCD1")
subset_icr_genes <- intersect(icr_genes[icr_genes %in% rownames(final_all_expr_df)],icr_genes[icr_genes %in% rownames(final_test_all_expr_df)])

train_icr_genes_expr <- t(final_all_expr_df[subset_icr_genes,rownames(train_df)])
test_icr_genes_expr <- t(final_test_all_expr_df[subset_icr_genes, rownames(test_df)])

#Make the final training and test datasets
rev_train_df <- as.data.frame(cbind(train_df[,c(1:4)],train_icr_genes_expr))
rev_test_df <- as.data.frame(cbind(test_df[,c(1:4)],test_icr_genes_expr))

#Build the Univariate models for each of the three constru tertiles
################################################################################
survival_analysis_function <- function(survival_df, subset_icr_genes)
{
  univ_formulas <- sapply(subset_icr_genes,
                        function(x) as.formula(paste('Surv(time, status)~', x)))

  univ_models <- lapply(univ_formulas, function(x){coxph(x, data = survival_df)})
  univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         p.value<-signif(x$wald["pvalue"], digits=2)
                         wald.test<-signif(x$wald["test"], digits=2)
                         beta<-signif(x$coef[1], digits=2); #coeficient beta
                         HR.mean <-signif(x$coef[2], digits=2); #exp(beta)
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         HR <- paste0(HR.mean, " (", 
                                      HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(beta, HR.mean, HR.confint.lower, HR.confint.upper, HR, wald.test, p.value)
                         names(res)<-c("beta", "mean", "low", "upper", "HR", "wald.test", 
                                       "p.value")
                         return(res)
                       })
  res <- t(as.data.frame(univ_results, check.names = FALSE))
  univariate_cox_results <- as.data.frame(res)
  univariate_cox_results$p.ajust <- p.adjust(univariate_cox_results$p.value, method="fdr")
  return(univariate_cox_results)
}

train_constru_low_univariate_df <- survival_analysis_function(survival_df = rev_train_df[rev_train_df$constru_tertiles==1,], subset_icr_genes)
train_constru_int_univariate_df <- survival_analysis_function(survival_df = rev_train_df[rev_train_df$constru_tertiles==2,], subset_icr_genes)
train_constru_high_univariate_df <- survival_analysis_function(survival_df = rev_train_df[rev_train_df$constru_tertiles==3,], subset_icr_genes)

test_constru_low_univariate_df <- survival_analysis_function(survival_df = rev_test_df[rev_test_df$constru_tertiles==1,], subset_icr_genes)
test_constru_int_univariate_df <- survival_analysis_function(survival_df = rev_test_df[rev_test_df$constru_tertiles==2,], subset_icr_genes)
test_constru_high_univariate_df <- survival_analysis_function(survival_df = rev_test_df[rev_test_df$constru_tertiles==3,], subset_icr_genes)

#Make the hazards matrix
column_split_values <- c("Con low","Con int","Con high")
train_constru_hazards_matrix <- as.matrix(cbind(as.numeric(train_constru_low_univariate_df$mean), as.numeric(train_constru_int_univariate_df$mean), as.numeric(train_constru_high_univariate_df$mean)))
colnames(train_constru_hazards_matrix) <- column_split_values
rownames(train_constru_hazards_matrix) <- rownames(train_constru_low_univariate_df)
train_constru_pval_matrix <- as.matrix(cbind(as.numeric(train_constru_low_univariate_df$p.ajust), as.numeric(train_constru_int_univariate_df$p.ajust), as.numeric(train_constru_high_univariate_df$p.ajust)))
colnames(train_constru_pval_matrix) <- column_split_values
rownames(train_constru_hazards_matrix) <- rownames(train_constru_low_univariate_df)

test_constru_hazards_matrix <- as.matrix(cbind(as.numeric(test_constru_low_univariate_df$mean), as.numeric(test_constru_int_univariate_df$mean), as.numeric(test_constru_high_univariate_df$mean)))
colnames(test_constru_hazards_matrix) <- column_split_values
rownames(test_constru_hazards_matrix) <- rownames(test_constru_low_univariate_df)
test_constru_pval_matrix <- as.matrix(cbind(as.numeric(test_constru_low_univariate_df$p.ajust), as.numeric(test_constru_int_univariate_df$p.ajust), as.numeric(test_constru_high_univariate_df$p.ajust)))
colnames(test_constru_pval_matrix) <- column_split_values
rownames(test_constru_hazards_matrix) <- rownames(test_constru_low_univariate_df)


column_split_values <- c("Con low","Con int","Con high")
column_split_values <- factor(column_split_values, levels = c("Con low","Con int","Con high"))
ha = HeatmapAnnotation(
  empty = anno_empty(border = FALSE, height = unit(8, "mm")),
  Con = anno_block(gp = gpar(fill = c("brown","orange","#FEDD00"))), 
  show_legend=F)
col_fun1 <- colorRamp2(c(0.5,1,1.5),c("blue","white","red"))
ht_icr_train = Heatmap(matrix=train_constru_hazards_matrix, col_fun1, 
                     name = "Hazards Ratio", column_title = qq("Hazards for ICR genes across Constru Groups in Training Set"),
                     cell_fun = function(j, i, x, y, width, height, fill) {
                       grid.rect(x = x, y = y, width = width, height = height, 
                                 gp = gpar(col = "grey", fill = NA))
                       if (train_constru_pval_matrix[i,j]<0.05 & train_constru_pval_matrix[i,j]>0.01)
                       {
                         grid.text(sprintf("*"), x, y, gp = gpar(fontsize = 11))
                       }
                       else if (train_constru_pval_matrix[i,j]<0.01 & train_constru_pval_matrix[i,j]>0.0001)
                       {
                         grid.text(sprintf("**"),x,y, gp=gpar(fontsize=11))
                       }
                       else if (train_constru_pval_matrix[i,j]<1e-4)
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
                     row_labels = rownames(train_constru_hazards_matrix),
                     row_names_max_width = max_text_width(rownames(train_constru_hazards_matrix)),
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
pdf("results/Training_Cox_Proportional_ICR_Genes.pdf",height = 8, width=4)
draw(ht_icr_train, heatmap_legend_side="bottom")
dev.off()

#Make the figure for the test set
ht_icr_test = Heatmap(matrix=test_constru_hazards_matrix, col_fun1, 
                     name = "Hazards Ratio", column_title = qq("Hazards for ICR genes across Constru Groups in Test Set"),
                     cell_fun = function(j, i, x, y, width, height, fill) {
                       grid.rect(x = x, y = y, width = width, height = height, 
                                 gp = gpar(col = "grey", fill = NA))
                       if (test_constru_pval_matrix[i,j]<0.05 & test_constru_pval_matrix[i,j]>0.01)
                       {
                         grid.text(sprintf("*"), x, y, gp = gpar(fontsize = 11))
                       }
                       else if (test_constru_pval_matrix[i,j]<0.01 & test_constru_pval_matrix[i,j]>0.0001)
                       {
                         grid.text(sprintf("**"),x,y, gp=gpar(fontsize=11))
                       }
                       else if (test_constru_pval_matrix[i,j]<1e-4)
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
                     row_labels = rownames(test_constru_hazards_matrix),
                     row_names_max_width = max_text_width(rownames(test_constru_hazards_matrix)),
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
pdf("results/Test_Cox_Proportional_ICR_Genes.pdf",height = 8, width=4)
draw(ht_icr_test, heatmap_legend_side="bottom")
dev.off()
