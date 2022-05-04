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
                              data = survival_df,  # data used to fit survival curves. 
                              risk.table = F,       # show risk table.
                              pval = TRUE,             # show p-value of log-rank test.
                              conf.int = FALSE,         # show confidence intervals for 
                              xlim = c(0,10),        # present narrower X axis, but not affect
                              # survival estimates.
                              break.time.by = 1,     # break X axis in time intervals by 500.
                              ggtheme = theme_classic(), # customize plot and risk table with a theme.
                              risk.table.y.text.col = T, # colour risk table text annotations.
                              risk.table.y.text = FALSE, # show bars instead of names in text annotations
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

#Make the CYT scores for TCGA
tcga_cyt_score <- as.numeric(as.vector(colMeans(log2_normalized_dataset[c("GZMA","PRF1"),])))
tcga_cyt_tertiles <- as.numeric(quantile(tcga_cyt_score, probs = c(0.33,0.66)))
tcga_cyt_tertile_values <- rep(1, length(tcga_cyt_score))
tcga_cyt_tertile_values[tcga_cyt_score>tcga_cyt_tertiles[1] & tcga_cyt_score<tcga_cyt_tertiles[2]] <- 2
tcga_cyt_tertile_values[tcga_cyt_score>tcga_cyt_tertiles[2]] <- 3

#Make the CONSTRU scores for TCGA
tcga_constru_score <- as.numeric(as.vector(colMeans(log2_normalized_dataset[upper_tertile_genes,])))-as.numeric(as.vector(colMeans(log2_normalized_dataset[lower_tertile_genes[lower_tertile_genes %in% rownames(log2_normalized_dataset)],])))
tcga_constru_tertiles <- as.numeric(quantile(tcga_constru_score, probs = c(0.33,0.66)))
tcga_constru_tertile_values <- rep(1, length(tcga_constru_score))
tcga_constru_tertile_values[tcga_constru_score>=tcga_constru_tertiles[1] & tcga_constru_score <=tcga_constru_tertiles[2]] <- 2
tcga_constru_tertile_values[tcga_constru_score>tcga_constru_tertiles[2]] <- 3

phenotype_df <- data.frame(CYT_tertiles = tcga_cyt_tertile_values, CONSTRU_tertiles = tcga_constru_tertile_values, Samples = colnames(log2_normalized_dataset))
survival_df <- fread("Data/survival_OV_survival.txt")
survival_df <- as.data.frame(survival_df)

#Make the revised survival df matching the sample ids in tcga with gene expression
revised_survival_df <- NULL
for (i in 1:nrow(survival_df))
{
  sample_id <- survival_df$sample[i]
  phenotype_df_id <- grep(sample_id, phenotype_df$Samples)
  if (length(phenotype_df_id)==1)
  {
    temp <- cbind(survival_df$sample[i],phenotype_df$Samples[phenotype_df_id],survival_df$OS[i],survival_df$OS.time[i], phenotype_df$CYT_tertiles[phenotype_df_id], phenotype_df$CONSTRU_tertiles[phenotype_df_id])
    revised_survival_df <- rbind(revised_survival_df, temp)
  }
}
revised_survival_df <- as.data.frame(revised_survival_df)
colnames(revised_survival_df) <- c("Sample","TCGA_Sample","OS","OS.time","CYT_Tertile","CONSTRU_Tertile")
for (i in 3:ncol(revised_survival_df))
{
  revised_survival_df[,i] <- as.numeric(as.vector(revised_survival_df[,i]))
}
revised_survival_df$CYT_Tertile <- as.factor(as.vector(revised_survival_df$CYT_Tertile))
revised_survival_df$CONSTRU_Tertile <- as.factor(as.vector(revised_survival_df$CONSTRU_Tertile))
revised_survival_df$OS.time <- revised_survival_df$OS.time/365

g_constru_low <- make_survival_plot(constru_df=revised_survival_df, tertile_value = 1)
ggsave(filename="results/KM_Plot_TCGA_Constru_Low.pdf",plot=g_constru_low$plot, device = pdf(), height = 4, width = 6, units="in", dpi=300)
dev.off()

g_constru_mid <- make_survival_plot(constru_df=revised_survival_df, tertile_value = 2)
ggsave(filename="results/KM_Plot_TCGA_Constru_Mid.pdf",plot=g_constru_mid$plot, device=pdf(), height=4, width=6, units="in", dpi=300)
dev.off()

g_constru_high <- make_survival_plot(constru_df = revised_survival_df, tertile_value = 3)
ggsave(filename="results/KM_Plot_TCGA_Constru_High.pdf", plot=g_constru_high$plot, device=pdf(), height=4, width=6, units="in", dpi = 300)
dev.off()