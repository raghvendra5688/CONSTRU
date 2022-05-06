library(data.table)
library(ggplot2)
library(Matrix)
library(stringr)
library(BiocManager)
library(GSVA)
library(BiocParallel)
library(doMC)
library(extrafontdb)
library(extrafont)
library(preprocessCore)
library(sva)
library(GEOquery)
library(survminer)
library(survival)
library(survivalAnalysis)
loadfonts()
registerDoMC(cores=20)

setwd("~/Documents/Misc_Work/Other Work/Miller_Related/CONSTRU/")

load("Data/OV_gene_RNAseq_normalized_TP_filtered.Rdata")

make_survival_plot <- function(constru_df, tertile_value)
{
  #Constru Low survival plot
  survival_df <- constru_df[constru_df$CONSTRU_Tertile==tertile_value,]
  model_fit <- survfit(Surv(OS.time, OS) ~ CYT_Tertile, data = survival_df)
  g <- ggsurvplot(model_fit,                     # survfit object with calculated statistics.
                  data = survival_df,            # data used to fit survival curves. 
                  risk.table = F,                # show risk table.
                  pval = TRUE,                   # show p-value of log-rank test.
                  conf.int = FALSE,              # show confidence intervals for 
                  xlim = c(0,10),                # present narrower X axis, but not affect
                  # survival estimates.
                  break.time.by = 1,             # break X axis in time intervals by 500.
                  ggtheme = theme_classic(),     # customize plot and risk table with a theme.
                  risk.table.y.text.col = T,     # colour risk table text annotations.
                  risk.table.y.text = FALSE,     # show bars instead of names in text annotations
                  palette = c("black","red","green"),
                  font.main = c(11, "bold","black"),
                  font.x = c(11, "bold","black"),
                  font.y = c(11, "bold","black"),
                  font.tickslab = c(11, "bold","black")
  )
  return(g)
  
}

#Get the constru upper tertile gene expression - constru lower tertile gene expression
upper_tertile_genes <- c("RAF1","TPD52","DDX21","UQCRB","TUBGCP4","RTF1","NCAPD3","DNAJC9","ZNF250","MEGF6","CDC42EP4","UBP1","C8orf33")
lower_tertile_genes <- c("PEPD","CAV2","DNAJB4","GALC","TRPC1","HAS2","PLAGL1","BAG2","AMACR","STAM2","PCGF1","WNT7A","MFGE8","ALB","WDFY3","NPTXR","OSR2","HOXA9",
                         "BDH2","CAMK2N1","PID1","PCOLCE2","FAT4","POGLUT2","EVA1B","CDC14B","APBB2")

log2_normalized_dataset <- log2(filtered.norm.RNAseqData+1)
samples_with_rnaseq <- colnames(log2_normalized_dataset)
samples_of_interest <- read.table("Data/TCGA_OV_Samples_of_Interest.csv",header=F)

#Make the CYT scores for TCGA
tcga_cyt_score <- as.numeric(as.vector(colMeans(log2_normalized_dataset[c("GZMA","PRF1"),])))

#Make the CONSTRU scores for TCGA
tcga_constru_score <- as.numeric(as.vector(colMeans(log2_normalized_dataset[upper_tertile_genes,])))-as.numeric(as.vector(colMeans(log2_normalized_dataset[lower_tertile_genes[lower_tertile_genes %in% rownames(log2_normalized_dataset)],])))

phenotype_df <- data.frame(CYT_scores = tcga_cyt_score, CONSTRU_scores = tcga_constru_score, Samples = samples_with_rnaseq)

#Get survival information
survival_df <- fread("Data/survival_OV_survival.txt")
survival_df <- as.data.frame(survival_df)

#Make the revised survival df matching the sample ids in tcga with gene expression
revised_survival_df <- NULL
for (i in 1:nrow(survival_df))
{
  sample_id <- survival_df$sample[i]
  phenotype_df_id <- grep(sample_id, samples_with_rnaseq)
  samples_interest_id <- grep(sample_id, samples_of_interest$V1)
  
  if (length(phenotype_df_id)==1 & length(samples_interest_id)==1)
  {
    temp <- cbind(survival_df$sample[i],phenotype_df$Samples[phenotype_df_id],survival_df$OS[i],survival_df$OS.time[i],phenotype_df$CYT_scores[phenotype_df_id],phenotype_df$CONSTRU_scores[phenotype_df_id])
    revised_survival_df <- rbind(revised_survival_df, temp)
  }
}
revised_survival_df <- as.data.frame(revised_survival_df)
colnames(revised_survival_df) <- c("Sample","TCGA_Sample","OS","OS.time","CYT_scores","CONSTRU_scores")
for (i in 3:ncol(revised_survival_df))
{
  revised_survival_df[,i] <- as.numeric(as.vector(revised_survival_df[,i]))
}
revised_survival_df$OS.time <- revised_survival_df$OS.time/365


#Make the cyt tertiles
rev_tcga_cyt_score <- revised_survival_df$CYT_scores
rev_tcga_cyt_tertiles <- as.numeric(quantile(rev_tcga_cyt_score, probs = c(0.33,0.66)))
rev_tcga_cyt_tertile_values <- rep(1, length(rev_tcga_cyt_score))
rev_tcga_cyt_tertile_values[rev_tcga_cyt_score>rev_tcga_cyt_tertiles[1] & rev_tcga_cyt_score< rev_tcga_cyt_tertiles[2]] <- 2
rev_tcga_cyt_tertile_values[rev_tcga_cyt_score>rev_tcga_cyt_tertiles[2]] <- 3
revised_survival_df$CYT_Tertile <- rev_tcga_cyt_tertile_values

#Make the constru tertiles
rev_tcga_constru_score <- revised_survival_df$CONSTRU_scores
rev_tcga_constru_tertiles <- as.numeric(quantile(rev_tcga_constru_score, probs = c(0.33,0.66)))
rev_tcga_constru_tertile_values <- rep(1, length(rev_tcga_constru_score))
rev_tcga_constru_tertile_values[rev_tcga_constru_score>=rev_tcga_constru_tertiles[1] & rev_tcga_constru_score <=rev_tcga_constru_tertiles[2]] <- 2
rev_tcga_constru_tertile_values[rev_tcga_constru_score>rev_tcga_constru_tertiles[2]] <- 3
revised_survival_df$CONSTRU_Tertile <- rev_tcga_constru_tertile_values

#Convert the constru and cyt tertiles into factors
revised_survival_df$CONSTRU_Tertile <- as.factor(as.vector(revised_survival_df$CONSTRU_Tertile))
revised_survival_df$CYT_Tertile <- as.factor(as.vector(revised_survival_df$CYT_Tertile))

#Make the KM plot
g_constru_low <- make_survival_plot(constru_df=revised_survival_df, tertile_value = 1)
ggsave(filename="results/KM_Plot_TCGA_Constru_Low.pdf",plot=g_constru_low$plot, device = pdf(), height = 4, width = 6, units="in", dpi=300)
dev.off()

g_constru_mid <- make_survival_plot(constru_df=revised_survival_df, tertile_value = 2)
ggsave(filename="results/KM_Plot_TCGA_Constru_Mid.pdf",plot=g_constru_mid$plot, device=pdf(), height=4, width=6, units="in", dpi=300)
dev.off()

g_constru_high <- make_survival_plot(constru_df = revised_survival_df, tertile_value = 3)
ggsave(filename="results/KM_Plot_TCGA_Constru_High.pdf", plot=g_constru_high$plot, device=pdf(), height=4, width=6, units="in", dpi = 300)
dev.off()

#Get TCGA clinical matrix
tcga_clinical_df <- fread("Data/TCGA.OV.sampleMap_OV_clinicalMatrix")
tcga_clinical_df <- as.data.frame(tcga_clinical_df)
subset_tcga_clinical_df <- tcga_clinical_df[tcga_clinical_df$sampleID %in% revised_survival_df$Sample,]
subset_tcga_clinical_df <- subset_tcga_clinical_df[order(subset_tcga_clinical_df$sampleID),]

#Match the clinical data with the survival data
revised_survival_df <- revised_survival_df[order(revised_survival_df$Sample),]
revised_survival_df$age <- subset_tcga_clinical_df$age_at_initial_pathologic_diagnosis
sample_name_matching_tmb <- paste(unlist(lapply(strsplit(revised_survival_df$Sample,split="-"),`[[`,1)),
                                  unlist(lapply(strsplit(revised_survival_df$Sample,split="-"),`[[`,2)),
                                  unlist(lapply(strsplit(revised_survival_df$Sample,split="-"),`[[`,3)),sep = "-")
revised_survival_df$Subset_Sample <- sample_name_matching_tmb

#Get the TMB information
tmb_info_df <- fread("Data/Thorsson_TCGA.csv")
tmb_info_df <- as.data.frame(tmb_info_df)
ov_tmb_info_df <- tmb_info_df[tmb_info_df$`TCGA Study`=="OV",]
rev_ov_tmb_info_df <- ov_tmb_info_df[ov_tmb_info_df$`TCGA Participant Barcode` %in% revised_survival_df$Subset_Sample,]
rev_ov_tmb_info_df <- rev_ov_tmb_info_df[order(rev_ov_tmb_info_df$`TCGA Participant Barcode`),]

#Match the tmb data with survival data
final_survival_df <- revised_survival_df[revised_survival_df$Subset_Sample %in% rev_ov_tmb_info_df$`TCGA Participant Barcode`,]
final_survival_df <- final_survival_df[order(final_survival_df$Subset_Sample),]
final_survival_df$Aneuploidy <- rev_ov_tmb_info_df$`Aneuploidy Score`

#Make the coxph model
res.cox <- coxph(Surv(OS.time, OS) ~ CYT_Tertile+Aneuploidy, data = final_survival_df[final_survival_df$CONSTRU_Tertile==3,])
multivariate_cox <- summary(res.cox)
multivariate_cox_df <- cox_as_data_frame(multivariate_cox,unmangle_dict = NULL,factor_id_sep = ":",sort_by = NULL)

#Get the constru high info
constru_high_survival_df <- final_survival_df[final_survival_df$CONSTRU_Tertile==3,]
g_box <- ggplot(data=constru_high_survival_df,aes(x=CYT_Tertile,y=Aneuploidy))+geom_boxplot()

