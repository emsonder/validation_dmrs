#' Methods used in several scripts

library(data.table)
library(pROC)

# Simes method, code borrowed from Stephany Orjuela: https://github.com/markrobinsonuzh/DAMEfinder/blob/master/R/simes_pval.R
simes <- function(pval) min((length(pval)*pval[order(pval)])/seq(from = 1, to = length(pval), by = 1))

# Function to calculate ROC curves
getRocs <- function(dt,
                    scores="mean_meth", 
                    labels="cond",
                    models="method",
                    posClass="cancer",
                    negClass="NM", 
                    subSample=FALSE){
  
  dt <- copy(dt)
  setnames(dt, scores, "scores")
  
  setorder(dt, -scores)
  if(subSample)
  {
    dt <- dt[sample(.N, 5*1e6)]
  }
  
  dt[,tpr:=cumsum(.SD==posClass)/sum(.SD==posClass), by=c(models), .SD=labels]
  dt[,fpr:=cumsum(.SD==negClass)/sum(.SD==negClass), by=c(models), .SD=labels]
  dt[,fdr:=cumsum(.SD==negClass)/seq_len(.N), by=c(models), .SD=labels]
  dt[,ppv:=cumsum(.SD==posClass)/seq_len(.N), by=c(models), .SD=labels]
  dt[,p:=seq_len(.N), by=c(models)]
  
  # Found no better solution for this
  setnames(dt, c(labels), c("labels"))
  dt[,auc_mod:=pROC::auc(pROC::roc(labels, scores)), by=c(models)]
  setnames(dt, c("labels"), c(labels))
  
  return(dt)
}

#Function to draw heatmap
drawHm <- function(obj, regions, colMeanCategory){
  
  gr_dmr <- regions[, !colnames(regions) %in% c("probe","category"), with=FALSE]
  
  # order probes by Category
  o <- order(regions$category)
  score <- as.matrix(gr_dmr)[o,]
  
  # Get colors I like
  col <- RColorBrewer::brewer.pal(n = 9, name = "YlGnBu")
  col_anot <- RColorBrewer::brewer.pal(n = 9, name = "Set1")
  purples <- RColorBrewer::brewer.pal(n = 3, name = "Purples")
  greens <- RColorBrewer::brewer.pal(n = 3, name = "Greens")
  set2 <- RColorBrewer::brewer.pal(n = 5, name = "Set2")
  
  # Set colors for methylation values in heatmap, yello high - blue low
  col_fun <- circlize::colorRamp2(c(0,0.1,1), c(col[9], col[7], col_anot[6]))
  
  # Set annotation colors
  #col_age <- purples[1:nlevels(obj$age_group)]
  #names(col_age) <- levels(obj$age_group)
  
  col_tis <- purples[1:nlevels(obj$condition)]
  names(col_tis) <- levels(obj$condition)
  
  col_cat <- set2
  #names(col_cat) <- regions$category
  names(col_cat) <- levels(factor(regions$category, 
                                  levels=unique(regions$category)))
  
  # Define column annotation
  column_ha <- HeatmapAnnotation(Tissue=obj$condition,
                                 col=list(Tissue=col_tis),
                                 gp=gpar(col="black"),
                                 CM=anno_barplot(colMeans(score[regions$category %in% colMeanCategory,], na.rm=TRUE), # this has a different order
                                                 gp=gpar(fill="#808080", 
                                                         col="#808080")))
  
  # Define row annotation 
  row_ha <- rowAnnotation(Category = as.character(regions$category)[o],
                          col = list(Category = col_cat))
  
  # Remove names
  #rownames(score) <- colnames(score) <- NULL
  rownames(score) <- NULL
  
  # Plot heatmap
  colnames(score) <- unlist(tstrsplit(colnames(score), "_", keep=1))
  hm <- Heatmap(score, 
                use_raster = TRUE,
                na_col = "white",
                column_split=obj$condition,
                row_split = factor(as.character(regions$category)[o],
                                   levels=c("Control regions", "Age dependent markers", 
                                            "cADN markers", "SSL markers", 
                                            "cADN and SSL specific markers")),
                row_order=o,
                top_annotation = column_ha,
                left_annotation = row_ha,
                col = viridis(100),#col_fun,
                clustering_distance_columns = "spearman",
                cluster_columns = TRUE,
                cluster_rows = FALSE,
                show_row_dend = FALSE,
                show_column_dend = FALSE,
                row_title = "Regions (n = 1096)", 
                column_title = "Samples",
                column_title_side = "bottom",
                column_names_gp = gpar(fontsize = 8),
                heatmap_legend_param = list(title = "methylation level",
                                            title_position = "lefttop-rot",
                                            grid_height = unit(1, "cm"),
                                            grid_width = unit(0.5, "cm")))
  return(hm)
}