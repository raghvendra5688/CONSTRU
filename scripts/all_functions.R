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
loadfonts()
registerDoMC(cores=20)

setwd(".")

#Make the heatmap for combined training and test set
make_pathway_activity_heatmap <- function(x1, x2, x3, final_all_expr_df, pathway_activities, immune_cell_type_concentrations, immune_activities){
  col_fun1 <- colorRamp2(c(-1,0,1.0),c("blue","white","red"))
  constru_tertiles <- c(x1[[3]],x2[[3]],x3[[3]])
  cyt_tertiles <- c(x1[[2]], x2[[2]], x3[[2]])
  
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
                       cluster_rows = F,
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
  return(list(ht_pathway, rev_immune_cell_type_concentrations, column_split_values, ha, rev_pathway_activities, rev_cyt_tertiles, rev_constru_tertiles))
}

#Make the immune concentration for training + test set
make_immune_conc_heatmap <- function(rev_immune_cell_type_concentrations, column_split_values, ha){
  col_fun2 = colorRamp2(c(0,0.5),c("white","red"))
  ht_immune_conc = Heatmap(matrix=as.matrix(rev_immune_cell_type_concentrations), col_fun2, 
                           name = "Immune Cell Fractions", column_title = qq("Immune Cell Fractions across Constru Groups"),
                           width = unit(15, "cm"),
                           height = unit(15, "cm"),
                           cluster_columns = F,
                           cluster_rows = F,
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
  return(ht_immune_conc)
}

#Perform the statistical test
perform_stat_test <- function(input_df)
{
  stat.test <- input_df %>%
    group_by(Constru, Pathway) %>%
    wilcox_test(Value~Cyt) %>%
    adjust_pvalue(method = "fdr") %>%
    add_significance("p.adj") %>% tibble() 
  stat.test <- as.data.frame(stat.test)
  out_table <- desc_statby(input_df, "Value", grps=c("Constru","Cyt","Pathway"))
  return(list(stat.test,out_table))
}

#Make the plot comparing mean and se of mean for all cell types w.r.t. phenotype
make_mean_se_stat_plot <- function(df, stat.test, plot_title, xlab_label, ylab_label)
{
  #Make the plot of a particular gene expression across cell types (mean + sd of mean) and pvalues for the disease status
  ###########################################################################################3
  bp <- ggbarplot(df, x = "Pathway", y = "Value", add = "mean_se", 
                  fill= "Cyt", palette = c("brown", "yellow"), facet.by = "Constru", nrow=3,
                  position = position_dodge(0.8), width = 0.80)
  
  # Add p-values onto the bar plots
  stat.test <- stat.test %>%
    add_xy_position(fun = "mean_se", x = "Pathway", group="Cyt", dodge = 0.8) 
  
  bp <- bp + stat_pvalue_manual(
    stat.test,  label = "p.adj.signif", tip.length = 0.01
  )
  
  # Move down the brackets using `bracket.nudge.y`
  bp + stat_pvalue_manual(
    stat.test,  label = "p.adj.signif", tip.length = 0,
    bracket.nudge.y = -2
  )
  
  bp <- bp+xlab(xlab_label) + ylab(ylab_label) + 
    theme_light()+
    theme(legend.text=element_text(size=16),
          legend.title=element_text(size=16))+ 
    ggtitle(plot_title)+
    theme(axis.text.x = element_text(color = "grey20", size = 10, angle = 90, hjust = 0.95, vjust = .5, face = "plain"),
          strip.text = element_text(color = "white", size=10, angle=0, hjust = 0.5, vjust = 0.5, face = "plain"),
          axis.text.y = element_text(color = "grey20", size = 10, angle = 0, hjust = 1, vjust = 0, face = "plain"),
          axis.title.x = element_text(color = "grey20", size = 16, angle = 0, hjust = .5, vjust = 0, face = "plain"),
          axis.title.y = element_text(color = "grey20", size = 16, angle = 90, hjust = .5, vjust = .5, face = "plain"),
          title=element_text(color="grey20", size=16, face="plain"))
  return(bp)
}

#Make the bubble plot for pathway activity
make_bubble_plot_heatmap <- function(activity_matrix, pval_matrix, tertile, col_fun3)
{
  if (tertile==1)
  {
    ha = HeatmapAnnotation(
      empty = anno_empty(border = FALSE, height = unit(8, "mm")),
      Constru = c(rep("brown",3)),
      Cyt = c("black","red","darkgreen"),
      col=list(Constru=c("brown"="brown"), 
               Cyt=c("black"="black","red"="white","darkgreen"="darkgreen")),
      simple_anno_size = unit(8, "mm"),
      annotation_name_gp = gpar(fontsize=11,col="white"),
      annotation_legend_param = list(directon="horizontal"),
      show_legend=F)
  }else if (tertile==2)
  {
    ha = HeatmapAnnotation(
      empty = anno_empty(border = FALSE, height = unit(8, "mm")),
      Constru = c(rep("orange",3)),
      Cyt = c("black","red","darkgreen"),
      col=list(Constru=c("orange"="orange"), 
               Cyt=c("black"="black","red"="white","darkgreen"="darkgreen")),
      simple_anno_size = unit(8, "mm"),
      annotation_name_gp = gpar(fontsize=11,col="white"),
      annotation_legend_param = list(directon="horizontal"),
      show_legend=F)
  }else if (tertile==3)
  {
    ha = HeatmapAnnotation(
      empty = anno_empty(border = FALSE, height = unit(8, "mm")),
      Constru = c(rep("#FEDD00",3)),
      Cyt = c("black","red","darkgreen"),
      col=list(Constru=c("#FEDD00"="#FEDD00"), 
               Cyt=c("black"="black","red"="white","darkgreen"="darkgreen")),
      simple_anno_size = unit(8, "mm"),
      annotation_name_gp = gpar(fontsize=11,col="white"),
      annotation_legend_param = list(directon="horizontal"),
      show_legend=F)
  }else if (tertile==4){
    ha = HeatmapAnnotation(
      empty = anno_empty(border = FALSE, height = unit(8, "mm")),
      Constru = anno_block(gp=gpar(fill=c("brown","white","#FEDD00"),col="white"),
                           labels=c("low","vs","high")),
      Cyt = c(rep("white",3)),
      col=list(Cyt=c("white"="white")),
      simple_anno_size = unit(8, "mm"),
      annotation_name_gp = gpar(fontsize=11,col="white"),
      annotation_legend_param = list(directon="horizontal"),
      show_legend=F)
  }else{
    ha = HeatmapAnnotation(
      empty = anno_empty(border = FALSE, height = unit(8, "mm")),
      Constru = anno_block(gp=gpar(fill=c("brown","white","#FEDD00"),col="white"),
                           labels=c("low","vs","high")),
      Cyt = c("black","red","darkgreen"),
      col=list(Cyt=c("black"="darkgreen","red"="white","darkgreen"="darkgreen")),
      simple_anno_size = unit(8, "mm"),
      annotation_name_gp = gpar(fontsize=11,col="white"),
      annotation_legend_param = list(directon="horizontal"),
      show_legend=F)
  }
  ht_pathway = Heatmap(matrix=activity_matrix, col_fun3, 
                       name="Diff Activity",
                       width = unit(3, "cm"),
                       height = unit(15, "cm"),
                       cell_fun = function(j, i, x, y, w, h, col) {
                         if (pval_matrix[i,j]>0.05)
                         {
                           r = min(unit.c(w, h))*0.0
                           grid.circle(x,y,r,gp=gpar(fill="white"))
                         }
                         else if(pval_matrix[i,j]<0.05 & pval_matrix[i,j]>0.01 & pval_matrix[i,j]>0) {
                           r = min(unit.c(w, h))*1
                           grid.circle(x, y, r, gp = gpar(fill = "red"))
                         } 
                         else if(pval_matrix[i,j]<0.01 & pval_matrix[i,j]>0.001 & pval_matrix[i,j]>0) {
                           r = min(unit.c(w, h))*2
                           grid.circle(x, y, r, gp = gpar(fill = "red"))
                         }
                         else if(pval_matrix[i,j]<0.001 & pval_matrix[i,j]>0) {
                           r = min(unit.c(w, h))*4
                           grid.circle(x, y, r, gp = gpar(fill = "red"))
                         }
                       },
                       cluster_columns = F,
                       cluster_rows = F,
                       row_names_centered = F,
                       show_row_names = F,
                       row_labels = rownames(activity_matrix),
                       row_names_max_width = max_text_width(rownames(activity_matrix)),
                       show_column_names = F,
                       column_split = c(1,2,3),
                       raster_quality = 2,
                       top_annotation = ha,
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
  return(ht_pathway)
}

#Get the data for making the heatmap of pathway activity
get_pathway_heatmap_with_comparisons <- function(train_out){
  
  rev_pathway_activities <- train_out[[5]]
  constru_tertiles <- train_out[[7]]
  cyt_tertiles <- train_out[[6]]
  
  mod_rev_pathway_activities <- as.data.frame(t(rbind(rbind(rev_pathway_activities,constru_tertiles),cyt_tertiles)))
  mod_rev_pathway_activities_df <- reshape2::melt(mod_rev_pathway_activities,id.vars=c("constru_tertiles","cyt_tertiles"))
  colnames(mod_rev_pathway_activities_df) <- c("Constru","Cyt","Pathway","Value")
  mod_rev_pathway_activities_df$Cyt <- as.factor(as.vector(mod_rev_pathway_activities_df$Cyt))
  mod_rev_pathway_activities_df$Constru <- as.factor(as.vector(mod_rev_pathway_activities_df$Constru))
  mod_rev_pathway_activities_df$Value <- as.numeric(as.vector(mod_rev_pathway_activities_df$Value))
  mod_rev_pathway_activities_df$Pathway <- as.character(as.vector(mod_rev_pathway_activities_df$Pathway))
  
  average_pathway_activity_matrix <- matrix(0, nrow=nrow(rev_pathway_activities), ncol=9)
  rownames(average_pathway_activity_matrix) <- rownames(rev_pathway_activities)
  constru_high_vs_low_ratio_matrix <- matrix(0, nrow=nrow(rev_pathway_activities), ncol=3)
  constru_high_vs_low_pval_matrix <- matrix(0, nrow=nrow(rev_pathway_activities), ncol=3)
  rownames(constru_high_vs_low_pval_matrix) <- rownames(rev_pathway_activities)
  rownames(constru_high_vs_low_ratio_matrix) <- rownames(rev_pathway_activities)
  
  cyt3_constru_high_vs_low_ratio_matrix <- matrix(0, nrow=nrow(rev_pathway_activities), ncol=3)
  cyt3_constru_high_vs_low_pval_matrix <- matrix(0, nrow=nrow(rev_pathway_activities), ncol=3)
  rownames(cyt3_constru_high_vs_low_ratio_matrix) <- rownames(rev_pathway_activities)
  rownames(cyt3_constru_high_vs_low_pval_matrix) <- rownames(rev_pathway_activities)
  
  unique_pathways <- rownames(average_pathway_activity_matrix)
  colnames(average_pathway_activity_matrix) <- NULL
  unique_constru_tertiles <- unique(mod_rev_pathway_activities_df$Constru)
  unique_cyt_tertiles <- unique(mod_rev_pathway_activities_df$Cyt)
  for (k in 1:length(unique_pathways))
  {
    pathway <- unique_pathways[k]
    fixed_pathway_df <- mod_rev_pathway_activities_df[mod_rev_pathway_activities_df$Pathway==pathway,]
    constru_high_activities <- fixed_pathway_df[fixed_pathway_df$Constru==3,]$Value
    constru_low_activities <- fixed_pathway_df[fixed_pathway_df$Constru==1,]$Value
    wx_test <- wilcox.test(constru_high_activities,constru_low_activities,paired=F)
    constru_high_vs_low_pval_matrix[pathway,2] <- wx_test$p.value
    constru_high_vs_low_ratio_matrix[pathway,2] <- mean(constru_high_activities)-mean(constru_low_activities)
    
    #Get the data comparing Constru High vs Low Cyt3 tertile
    cyt3_constru_high_activities <- fixed_pathway_df[fixed_pathway_df$Constru==3 & fixed_pathway_df$Cyt==3,]$Value
    cyt3_constru_low_activities <- fixed_pathway_df[fixed_pathway_df$Constru==1 & fixed_pathway_df$Cyt==3,]$Value
    rev_wx_test <- wilcox.test(cyt3_constru_high_activities,cyt3_constru_low_activities,paired=F)
    
    cyt3_constru_high_vs_low_pval_matrix[pathway,2] <- rev_wx_test$p.value
    cyt3_constru_high_vs_low_ratio_matrix[pathway,2] <- mean(cyt3_constru_high_activities)-mean(cyt3_constru_low_activities)
    for (i in 1:length(unique_constru_tertiles))
    {
      constru_tertile <- unique_constru_tertiles[i]
      
      for (j in 1:length(unique_cyt_tertiles))
      {
        cyt_tertile <- unique_cyt_tertiles[j]
        average_pathway_activity_matrix[pathway,((i-1)*3+j)] <- mean(mod_rev_pathway_activities_df[mod_rev_pathway_activities_df$Pathway==pathway & mod_rev_pathway_activities_df$Constru==constru_tertile 
                                                                                                   & mod_rev_pathway_activities_df$Cyt==cyt_tertile,]$Value)
      }
    }
  }
  constru_high_vs_low_pval_matrix[,2] <- p.adjust(constru_high_vs_low_pval_matrix[,2], method="fdr")
  cyt3_constru_high_vs_low_pval_matrix[,2] <- p.adjust(cyt3_constru_high_vs_low_pval_matrix[,2], method="fdr")
  
  return(list(average_pathway_activity_matrix, constru_high_vs_low_ratio_matrix, constru_high_vs_low_pval_matrix, 
              cyt3_constru_high_vs_low_ratio_matrix, cyt3_constru_high_vs_low_pval_matrix))
}

#Get the heatmap for immune concentrations
get_immune_conc_heatmap_with_comparisons <- function(train_out){
  
  rev_immune_cell_type_concentrations <- train_out[[2]]
  constru_tertiles <- train_out[[7]]
  cyt_tertiles <- train_out[[6]]
  
  mod_rev_immune_conc <- as.data.frame(t(rbind(rbind(rev_immune_cell_type_concentrations,constru_tertiles),cyt_tertiles)))
  colnames(mod_rev_immune_conc)[c(40:41)] <- c("constru_tertiles","cyt_tertiles")
  mod_rev_immune_conc_df <- reshape2::melt(mod_rev_immune_conc,id.vars=c("constru_tertiles","cyt_tertiles"))
  colnames(mod_rev_immune_conc_df) <- c("Constru","Cyt","Immune_Celltype","Value")
  mod_rev_immune_conc_df$Cyt <- as.factor(as.vector(mod_rev_immune_conc_df$Cyt))
  mod_rev_immune_conc_df$Constru <- as.factor(as.vector(mod_rev_immune_conc_df$Constru))
  mod_rev_immune_conc_df$Value <- as.numeric(as.vector(mod_rev_immune_conc_df$Value))
  mod_rev_immune_conc_df$Immune_Celltype <- as.character(as.vector(mod_rev_immune_conc_df$Immune_Celltype))
  
  
  average_immune_conc_matrix <- matrix(0, nrow=nrow(rev_immune_cell_type_concentrations), ncol=9)
  rownames(average_immune_conc_matrix) <- rownames(rev_immune_cell_type_concentrations)
  immune_constru_high_vs_low_ratio_matrix <- matrix(1, nrow=nrow(rev_immune_cell_type_concentrations), ncol=3)
  immune_constru_high_vs_low_pval_matrix <- matrix(1, nrow=nrow(rev_immune_cell_type_concentrations), ncol=3)
  immune_cyt3_constru_high_vs_low_ratio_matrix <- matrix(1, nrow=nrow(rev_immune_cell_type_concentrations), ncol=3)
  immune_cyt3_constru_high_vs_low_pval_matrix <- matrix(1, nrow=nrow(rev_immune_cell_type_concentrations), ncol=3)
  rownames(immune_constru_high_vs_low_pval_matrix) <- rownames(rev_immune_cell_type_concentrations)
  rownames(immune_constru_high_vs_low_ratio_matrix) <- rownames(rev_immune_cell_type_concentrations)
  rownames(immune_cyt3_constru_high_vs_low_ratio_matrix) <- rownames(rev_immune_cell_type_concentrations)
  rownames(immune_cyt3_constru_high_vs_low_pval_matrix) <- rownames(rev_immune_cell_type_concentrations)
  
  unique_immune_conc <- rownames(average_immune_conc_matrix)
  colnames(average_immune_conc_matrix) <- NULL
  
  unique_constru_tertiles <- unique(mod_rev_immune_conc_df$Constru)
  unique_cyt_tertiles <- unique(mod_rev_immune_conc_df$Cyt)
  
  for (k in 1:length(unique_immune_conc))
  {
    immune_type <- unique_immune_conc[k]
    fixed_immune_conc_df <- mod_rev_immune_conc_df[mod_rev_immune_conc_df$Immune_Celltype==immune_type,]
    constru_high_conc <- fixed_immune_conc_df[fixed_immune_conc_df$Constru==3,]$Value
    constru_low_conc <- fixed_immune_conc_df[fixed_immune_conc_df$Constru==1,]$Value
    wx_test <- wilcox.test(constru_high_conc,constru_low_conc,paired=F)
    immune_constru_high_vs_low_pval_matrix[immune_type,2] <- wx_test$p.value
    immune_constru_high_vs_low_ratio_matrix[immune_type,2] <- mean(constru_high_conc)/mean(constru_low_conc)
    
    for (i in 1:length(unique_constru_tertiles))
    {
      cyt_constru_high_conc <- fixed_immune_conc_df[fixed_immune_conc_df$Constru==3 & fixed_immune_conc_df$Cyt==3,]$Value
      cyt_constru_low_conc <- fixed_immune_conc_df[fixed_immune_conc_df$Constru==1 & fixed_immune_conc_df$Cyt==3,]$Value
      rev_wx_test <- wilcox.test(cyt_constru_high_conc,cyt_constru_low_conc,paired=F,exact=F)
      
      immune_cyt3_constru_high_vs_low_pval_matrix[immune_type,2] <- rev_wx_test$p.value
      immune_cyt3_constru_high_vs_low_ratio_matrix[immune_type,2] <- mean(cyt_constru_high_conc)/mean(cyt_constru_low_conc)
      constru_tertile <- unique_constru_tertiles[i]
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
  immune_cyt3_constru_high_vs_low_pval_matrix[,2] <- p.adjust(immune_cyt3_constru_high_vs_low_pval_matrix[,2], method="fdr")
  
  return(list(average_immune_conc_matrix, immune_constru_high_vs_low_ratio_matrix, immune_constru_high_vs_low_pval_matrix, 
              immune_cyt3_constru_high_vs_low_ratio_matrix, immune_cyt3_constru_high_vs_low_pval_matrix))
}

#Make the combined heatmap of avg pathway activity of training vs test set
get_all_heatmaps <- function(train_pathway_out)
{
  average_pathway_activity_matrix <- train_pathway_out[[1]]
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
  ht_train_rev_pathway = Heatmap(matrix= average_pathway_activity_matrix, col_fun3, 
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
  
  ht_train_constru_mean_pathway <- make_bubble_plot_heatmap(activity_matrix = train_pathway_out[[2]], train_pathway_out[[3]], tertile=4, col_fun3)
  ht_train_constru_cyt3_pathway <- make_bubble_plot_heatmap(activity_matrix = train_pathway_out[[4]], train_pathway_out[[5]], tertile=5, col_fun3)
  return(list(ht_train_rev_pathway, ht_train_constru_mean_pathway, ht_train_constru_cyt3_pathway))
}

#Make the bubble plot of heatmap for immune concentrations
make_bubble_plot_heatmap_immune <- function(immune_constru_high_vs_low_ratio_matrix, immune_constru_high_vs_low_pval_matrix, tertile)
{
  if (tertile==4){
    ha2 = HeatmapAnnotation(
      empty = anno_empty(border = FALSE, height = unit(8, "mm")),
      Constru = anno_block(gp=gpar(fill=c("brown","white","#FEDD00"),col="white"),
                           labels=c("low","vs","high")),
      Cyt = c(rep("white",3)),
      col=list(Cyt=c("white"="white")),
      simple_anno_size = unit(8, "mm"),
      annotation_name_gp = gpar(fontsize=11,col="white"),
      annotation_legend_param = list(directon="horizontal"),
      show_legend=F)
  }else{
    ha2 = HeatmapAnnotation(
      empty = anno_empty(border = FALSE, height = unit(8, "mm")),
      Constru = anno_block(gp=gpar(fill=c("brown","white","#FEDD00"),col="white"),
                           labels=c("low","vs","high")),
      Cyt = c("black","red","darkgreen"),
      col=list(Cyt=c("black"="darkgreen","red"="white","darkgreen"="darkgreen")),
      simple_anno_size = unit(8, "mm"),
      annotation_name_gp = gpar(fontsize=11,col="white"),
      annotation_legend_param = list(directon="horizontal"),
      show_legend=F)
  }
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
  return(ht_constru_mean_immune_conc)
}

#Make the combined heatmp for immune cell concentration for training + test
get_immune_all_heatmaps <- function(train_immune_conc_out)
{
  average_immune_conc_matrix <- train_immune_conc_out[[1]]
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
  ht_train_constru_mean_immune_conc <- make_bubble_plot_heatmap_immune(train_immune_conc_out[[2]], train_immune_conc_out[[3]], tertile=4)
  ht_train_constru_cyt3_immune_conc <- make_bubble_plot_heatmap_immune(train_immune_conc_out[[4]], train_immune_conc_out[[5]], tertile=5)
  return(list(ht_rev_immune_conc, ht_train_constru_mean_immune_conc, ht_train_constru_cyt3_immune_conc))
}

get_cox_info <- function(x){
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
  return(res)
}