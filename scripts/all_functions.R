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
  }else{
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
                           r = min(unit.c(w, h))*0.5
                           grid.circle(x, y, r, gp = gpar(fill = "red"))
                         } 
                         else if(pval_matrix[i,j]<0.01 & pval_matrix[i,j]>0.001 & pval_matrix[i,j]>0) {
                           r = min(unit.c(w, h))*1
                           grid.circle(x, y, r, gp = gpar(fill = "red"))
                         }
                         else if(pval_matrix[i,j]<0.001 & pval_matrix[i,j]>0) {
                           r = min(unit.c(w, h))*2
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
