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
