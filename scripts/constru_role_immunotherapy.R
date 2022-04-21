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
library(GEOquery)
loadfonts()
registerDoMC(cores=20)

setwd("~/Documents/Misc_Work/Other Work/Miller_Related/CONSTRU/")
source("scripts/all_functions.R")
load("Data/Master_Datasets.RData")

#Get the constru upper tertile gene expression - constru lower tertile gene expression
upper_tertile_genes <- c("RAF1","TPD52","DDX21","UQCRB","TUBGCP4","RTF1","NCAPD3","DNAJC9","ZNF250","MEGF6","CDC42EP4","UBP1","C8orf33")
lower_tertile_genes <- c("PEPD","CAV2","DNAJB4","GALC","TRPC1","HAS2","PLAGL1","BAG2","AMACR","STAM2","PCGF1","WNT7A","MFGE8","ALB","WDFY3","NPTXR","OSR2","HOXA9",
                         "BDH2","CAMK2N1","PID1","PCOLCE2","FAT4","POGLUT2","EVA1B","CDC14B","APBB2")

#Perform the analysis for gse67501 Ascierto
###############################################################################
z_gse67501 <- getGEO("GSE67501")[[1]]
gse67501_expr <- z_gse67501@assayData$exprs
gse65701_fData <- fData(z_gse67501)
gse67501_pData <- pData(z_gse67501)
all_genes <- gse65701_fData$Symbol

gse67501_CONSTRU_scores <- get_constru_scores(gse67501_expr = gse67501_expr, gse67501_pData = gse67501_pData, gse65701_fData = gse65701_fData,
                                              all_genes = all_genes, upper_tertile_genes, lower_tertile_genes)

#Get the constru and icr tertiles
Ascierto_GSE67501$CONSTRU <- gse67501_CONSTRU_scores
Ascierto_GSE67501$Response <- str_replace(Ascierto_GSE67501$Response," ","")
icr_tertiles <- quantile(Ascierto_GSE67501$ICR,probs=c(0.33,0.66))
gse67501_CONSTRU_tertiles <- quantile(gse67501_CONSTRU_scores, probs = c(0.33,0.66))

g_gse67501 <- make_icr_vs_constru_plot(data_df = Ascierto_GSE67501, icr_tertiles = icr_tertiles, constru_tertiles = gse67501_CONSTRU_tertiles, titlename = "Ascierto GSE65701 (RCC)")
ggsave(filename = "results/Ascierto_GSE67501_ICR_vs_CONSTRU.pdf", plot=g_gse67501, device=pdf(), height=6, width = 8, units="in", dpi = 300)
dev.off()

#Perform the analysis for gse79691
#################################################################################
z_gse79691 <- getGEO("GSE79691")[[1]]
gse79691_expr <- z_gse79691@assayData$exprs
gse79691_fData <- fData(z_gse79691)
gse79691_pData <- pData(z_gse79691)
all_genes <- gse79691_fData$Symbol

gse79691_CONSTRU_scores <- get_constru_scores(gse67501_expr = gse79691_expr, gse67501_pData = gse79691_pData, gse65701_fData = gse79691_fData, 
                                              all_genes = all_genes, upper_tertile_genes, lower_tertile_genes)

#Get the constru and icr tertiles
Ascierto_GSE79691$CONSTRU <- gse79691_CONSTRU_scores
Ascierto_GSE79691$Response2 <- "R"
Ascierto_GSE79691[Ascierto_GSE79691$response.to.anti.pd.1..nivolumab..immunotherapy..regression.or.progression.=="Progression",]$Response2 <- "NR"
Ascierto_GSE79691$Response <- "NR"
Ascierto_GSE79691[Ascierto_GSE79691$Response2=="R",]$Response <- "CR"
gse79691_icr_tertiles <- quantile(Ascierto_GSE79691$ICR,probs=c(0.33,0.66))
gse79691_constru_tertiles <- quantile(Ascierto_GSE79691$CONSTRU, probs = c(0.33, 0.66))

g_gse79761 <- make_icr_vs_constru_plot(data_df = Ascierto_GSE79691, icr_tertiles = gse79691_icr_tertiles, constru_tertiles = gse79691_constru_tertiles, titlename = "Ascierto GSE79691 (Melanoma)")
ggsave(filename="results/Ascierto_GSE79691_ICR_vs_CONSTRU.pdf",plot=g_gse79761, device=pdf(), height=6, width=8, units="in",dpi=300)
dev.off()


#Perform the analysis for GSE115821
################################################################################
z_gse115821 <- getGEO("GSE115821")[[1]]
gse115821_df <- fread("Data/GSE115821_MGH_counts.csv")
gse115821_df <- as.data.frame(gse115821_df)
gse115821_expr_df <- gse115821_df[,c(7:ncol(gse115821_df))]
gse115821_expr_mat <- normalize.quantiles(as.matrix(log2(gse115821_expr_df+1)))
colnames(gse115821_expr_mat) <- colnames(gse115821_expr_df)
all_genes <- gse115821_df$Geneid
gse115821_fData <- gse115821_df[,c(1:6)]
gse115821_pData <- Auslander_GSE115821[,c(1:10)]

gse115821_CONSTRU_scores <- get_constru_scores(gse67501_expr = gse115821_expr_mat, gse67501_pData = gse115821_pData, gse65701_fData = gse115821_fData,
                                               all_genes = all_genes, upper_tertile_genes, lower_tertile_genes)

#Get the constru and icr tertiles
Auslander_GSE115821$CONSTRU <- gse115821_CONSTRU_scores
Auslander_GSE115821$Response2 <- Auslander_GSE115821$response
Auslander_GSE115821$Response <- Auslander_GSE115821$Response2
Auslander_GSE115821[Auslander_GSE115821$Response2==" R",]$Response <- "CR"
gse115821_icr_tertiles <- quantile(Auslander_GSE115821$ICR, probs=c(0.33,0.66))
gse115821_constru_tertiles <- quantile(Auslander_GSE115821$CONSTRU, probs=c(0.33,0.66))

g_gse115821 <- make_icr_vs_constru_plot(data_df = Auslander_GSE115821, icr_tertiles = gse115821_icr_tertiles, constru_tertiles = gse115821_constru_tertiles, titlename="Auslander GSE115821 (Melanoma)")
ggsave(filename="results/Auslander_GSE115821_ICR_vs_CONSTRU.pdf", plot=g_gse115821, device = pdf(), height=6, width =8, units="in", dpi=300)
dev.off()

#Perform the analysis on Chen dataset
################################################################################
z_gse126044 <- getGEO("GSE126044")[[1]]
gse126044_df <- fread("Data/GSE126044_counts.txt")
gse126044_df <- as.data.frame(gse126044_df)
gse126044_expr_df <- gse126044_df[,c(2:ncol(gse126044_df))]
gse126044_expr_mat <- normalize.quantiles(as.matrix(log2(gse126044_expr_df+1)))
all_genes <- gse126044_df$V1
gse126044_fData <- all_genes
gse126044_pData <- Cho_GSE126044[,c(1:10)]

gse126044_CONSTRU_scores <- get_constru_scores(gse67501_expr = gse126044_expr_mat, gse67501_pData = gse126044_pData, gse65701_fData = gse126044_fData,
                                               all_genes = all_genes, upper_tertile_genes, lower_tertile_genes)

#Get the constru and icr tertiles
Cho_GSE126044$CONSTRU <- gse126044_CONSTRU_scores
Cho_GSE126044$Response2 <- "NR"
Cho_GSE126044[Cho_GSE126044$Response=="Response",]$Response2 <- "R"
Cho_GSE126044$Response <- "NR"
Cho_GSE126044[Cho_GSE126044$Response2=="R",]$Response <- "CR"
gse126044_icr_tertiles <- quantile(Cho_GSE126044$ICR, probs = c(0.33, 0.66))
gse126044_constru_tertiles <- quantile(Cho_GSE126044$CONSTRU, probs = c(0.33, 0.66))

g_gse126044 <- make_icr_vs_constru_plot(data_df = Cho_GSE126044, icr_tertiles = gse126044_icr_tertiles, constru_tertiles = gse126044_constru_tertiles, titlename="Cho GSE126044 (NSCLC)")
ggsave(filename="results/Cho_GSE126044_ICR_vs_CONSTRU.pdf", plot=g_gse126044, device = pdf(), height=6, width=8, units="in", dpi=300)
dev.off()

#Perform analysis on Gide dataset
################################################################################
gide_df <- fread("Data/DATASET-PRJEB23709_Pre_73samples_Gide.txt")
gide_df <- as.data.frame(gide_df)
gide_expr_df <- gide_df[,c(2:75)]
subset_Gide_PRJEB23709 <- Gide_PRJEB23709[colnames(gide_expr_df),]
all_genes <- gide_df$Symbol
gide_fData <- all_genes
gide_pData <- subset_Gide_PRJEB23709[,c(1:10)]

gide_CONSTRU_scores <- get_constru_scores(gse67501_expr = gide_expr_df, gse67501_pData = gide_pData, gse65701_fData = gide_fData,
                                          all_genes = all_genes, upper_tertile_genes, lower_tertile_genes)

#Get the constru and icr tertiles
subset_Gide_PRJEB23709$CONSTRU <- gide_CONSTRU_scores
subset_Gide_PRJEB23709 <- subset_Gide_PRJEB23709[!is.na(subset_Gide_PRJEB23709$Status),]
subset_Gide_PRJEB23709$Response2 <- "NR"
subset_Gide_PRJEB23709[subset_Gide_PRJEB23709$Status==1,]$Response2 <- "R"
subset_Gide_PRJEB23709$Response <- "NR"
subset_Gide_PRJEB23709[subset_Gide_PRJEB23709$Response2=="R",]$Response <- "CR"
gide_icr_tertiles <- quantile(subset_Gide_PRJEB23709$ICR, probs = c(0.33, 0.66))
gide_constru_tertiles <- quantile(subset_Gide_PRJEB23709$CONSTRU, probs = c(0.33, 0.66))

g_gide <- make_icr_vs_constru_plot(data_df = subset_Gide_PRJEB23709, icr_tertiles = gide_icr_tertiles, constru_tertiles = gide_constru_tertiles, titlename = "Gide (Melanoma)")
ggsave(filename="results/Gide_Melanoma_ICR_vs_CONSTRU.pdf", plot=g_gide, device=pdf(), height=6, width=8, units="in", dpi=300)
dev.off()


#Perform analysis on Riaz
################################################################################
load("Data/002.Riaz_Data_normalized.gene.counts.Rdata")
riaz_expr_df <- dataNorm
all_genes <- rownames(riaz_expr_df)
riaz_fData <- all_genes
riaz_pData <- Riaz_GSE91061[,c(1:10)]

riaz_expr_mat <- log2(riaz_expr_df+1)
riaz_CONSTRU_scores <- get_constru_scores(gse67501_expr = riaz_expr_mat, gse67501_pData = riaz_pData, gse65701_fData = riaz_fData,
                                          all_genes = all_genes, upper_tertile_genes, lower_tertile_genes)

#Get the constru and icr tertiles
Riaz_GSE91061$CONSTRU <- riaz_CONSTRU_scores
subset_Riaz_GSE91061 <- Riaz_GSE91061[!is.na(Riaz_GSE91061$Status),]
subset_Riaz_GSE91061$Response2 <- "NR"
subset_Riaz_GSE91061[subset_Riaz_GSE91061$Status==1,]$Response2 <- "R"
subset_Riaz_GSE91061$Response <- "NR"
subset_Riaz_GSE91061[subset_Riaz_GSE91061$Status==1,]$Response <- "CR"
riaz_icr_tertiles <- quantile(subset_Riaz_GSE91061$ICR, probs=c(0.33, 0.66))
riaz_constru_tertiles <- quantile(subset_Riaz_GSE91061$CONSTRU, probs = c(0.33, 0.66))

g_riaz <- make_icr_vs_constru_plot(data_df = subset_Riaz_GSE91061, icr_tertiles = riaz_icr_tertiles, constru_tertiles = riaz_constru_tertiles, titlename = "Riaz (Melanoma)")
ggsave(filename="results/Riaz_Melanoma_ICR_vs_CONSTRU.pdf", plot=g_riaz, device = pdf(), height=6, width=8, units="in", dpi=300)
dev.off()

#Perform analysis on RIbas
################################################################################
z_gse78220 <- getGEO("GSE78220")[[1]]
gse78220_df <- fread("Data/GSE78220_Ribas_FPKM.csv")
gse78220_df <- as.data.frame(gse78220_df)
gse78220_expr_df <- gse78220_df[,c(2:ncol(gse78220_df))]
all_genes <- gse78220_df$Gene
gse78220_pData <- Ribas_GSE78220[,c(1:10)]
gse78220_fData <- all_genes
gse78220_expr_mat <- normalize.quantiles(as.matrix(log2(gse78220_expr_df+1)))
colnames(gse78220_expr_mat) <- colnames(gse78220_expr_df)

ribas_CONSTRU_scores <- get_constru_scores(gse67501_expr = gse78220_expr_mat, gse67501_pData = gse78220_pData, gse65701_fData = gse78220_fData,
                                           all_genes = all_genes, upper_tertile_genes, lower_tertile_genes)

#Get the constru and icr tertiles
Ribas_GSE78220$CONSTRU <- ribas_CONSTRU_scores
Ribas_GSE78220$Response2 <- "NR"
Ribas_GSE78220[Ribas_GSE78220$Status==1,]$Response2 <- "R"
Ribas_GSE78220$Response <- "NR"
Ribas_GSE78220[Ribas_GSE78220$Status==1,]$Response <- "CR"
ribas_icr_tertiles <- quantile(Ribas_GSE78220$ICR, probs = c(0.33,0.66))
ribas_constru_tertiles <- quantile(Ribas_GSE78220$CONSTRU, probs = c(0.33,0.66))

g_ribas <- make_icr_vs_constru_plot(data_df = Ribas_GSE78220, icr_tertiles = ribas_icr_tertiles, constru_tertiles = ribas_constru_tertiles, title = "Ribas (Melanoma)")
ggsave(filename="results/Ribas_Melanoma_ICR_vs_CONSTRU.pdf",plot=g_ribas, device=pdf(), height=6, width=8, units="in", dpi =300)
dev.off()

#Perform analysis on UlloaMontoya
################################################################################
z_gse35640 <- getGEO("GSE35640")[[1]]
gse35640_expr_df <- z_gse35640@assayData$exprs
gse35640_pData <- pData(z_gse35640)
gse35640_fData <- fData(z_gse35640)
all_genes <- gse35640_fData$`Gene Symbol`

gse35640_CONSTRU_scores <- get_constru_scores(gse67501_expr = gse35640_expr_df, gse67501_pData = gse35640_pData, gse65701_fData = gse35640_fData,
                                              all_genes=all_genes, upper_tertile_genes, lower_tertile_genes)

#Get the constru and icr tertiles
UlloaMontoya_GSE35640$CONSTRU <- gse35640_CONSTRU_scores
UlloaMontoya_GSE35640 <-  UlloaMontoya_GSE35640[!is.na(UlloaMontoya_GSE35640$Status),]
UlloaMontoya_GSE35640$Response2 <- "NR"
UlloaMontoya_GSE35640[UlloaMontoya_GSE35640$Status==1,]$Response2 <- "R"
UlloaMontoya_GSE35640$Response <- "NR"
UlloaMontoya_GSE35640[UlloaMontoya_GSE35640$Status==1,]$Response <- "CR"
ulloamontoya_icr_tertiles <- quantile(UlloaMontoya_GSE35640$ICR, probs = c(0.33,0.66))
ulloamontoya_constru_tertiles <- quantile(UlloaMontoya_GSE35640$CONSTRU, probs = c(0.33, 0.66))

g_ulloa <- make_icr_vs_constru_plot(data_df = UlloaMontoya_GSE35640, icr_tertiles = ulloamontoya_icr_tertiles, constru_tertiles = ulloamontoya_constru_tertiles, title = "Ulloa Montoya (Melanoma)")
ggsave(filename="results/UlloaMontoya_Melanoma_ICR_vs_CONSTRU.pdf", plot = g_ulloa, device=pdf(), height=6, width=8, units="in", dpi=300)
dev.off()

#Perform analysis on Glioblastoma
###############################################################################
z_gse121810 <- getGEO("GSE121810")[[1]]
gse121810_df <- fread("Data/GSE121810_Prins.PD1NeoAdjv.Jul2018.HUGO.PtID.csv")
gse121810_df <- as.data.frame(gse121810_df)
gse121810_expr_df <- gse121810_df[,c(2:ncol(gse121810_df))]
all_genes <- gse121810_df$Genes
gse121810_fData <- all_genes
gse121810_pData <- pData(z_gse121810)
gse121810_expr_mat <- normalize.quantiles(as.matrix(log2(gse121810_expr_df+1)))

gse121810_CONSTRU_scores <- get_constru_scores(gse67501_expr = gse121810_expr_mat, gse67501_pData = gse121810_pData, gse65701_fData = gse121810_fData,
                                               all_genes, upper_tertile_genes, lower_tertile_genes)
gse121810_final_expr_df <- get_expression_df(gse67501_expr = gse121810_expr_mat, all_genes)
gse121810_final_expr_mat <- as.matrix(gse121810_final_expr_df)
colnames(gse121810_final_expr_mat) <- colnames(gse121810_final_expr_df)
rownames(gse121810_final_expr_mat) <- rownames(gse121810_final_expr_df)

#Get the constru and icr tertiles
icr_scores <- gsva(gse121810_final_expr_mat, gset.idx.list = list(c("CXCL10","CXCL9","CCL5","IFNG","IL12B","TBX21","CD8A","CD8B","STAT1","IRF1","GNLY",
                                                           "PRF1","GZMA","GZMB","GZMH","CD274","CTLA4","FOXP3","IDO1","PDCD1")), method="ssgsea", kcdf=c("Gaussian"),ssgsea.norm=T)
icr_scores <- as.numeric(as.vector(icr_scores))
gse121810_pData$CONSTRU <- gse121810_CONSTRU_scores
gse121810_pData$ICR <- icr_scores
gse121810_pData$Response2 <- "NR"
gse121810_pData[gse121810_pData$`therapy:ch1`=="neoadjuvant pembrolizumab",]$Response2 <- "R"
gse121810_pData$Response <- "NR"
gse121810_pData[gse121810_pData$`therapy:ch1`=="neoadjuvant pembrolizumab",]$Response <- "CR"
gse121810_icr_tertiles <- quantile(gse121810_pData$ICR, probs = c(0.33, 0.66))
gse121810_constru_tertiles <- quantile(gse121810_pData$CONSTRU, probs = c(0.33, 0.66))

g_gse121810 <- make_icr_vs_constru_plot(data_df = gse121810_pData, icr_tertiles = gse121810_icr_tertiles, constru_tertiles = gse121810_constru_tertiles, title = "GSE121810 (Glioblastoma)")
ggsave(filename = "results/NatMedicine_Glioblastoma_ICR_vs_CONSTRU.pdf", plot=g_gse121810, device = pdf(), height = 6, width=8, units="in", dpi =300)
dev.off()