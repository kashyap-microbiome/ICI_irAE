zicoseq_heatmap <- function(zico_list, grp.names, grp.labels, sign_cutoffs, single_cutoff = FALSE, p_or_q, column_name_rot = 0, 
                            rotate=FALSE, left_margin=4, new_names=FALSE, new_names_vec=NA) {
  
  require(ComplexHeatmap)
  require(circlize)
  

  q.cut1 <- min(sign_cutoffs); q.cut2 <- median(sign_cutoffs); q.cut3 <- max(sign_cutoffs)
  max.cutoff <- max(sign_cutoffs)
  
  names(zico_list) <- grp.names
  
  FDR_list <- zico_list
  p_list <- zico_list
  R2_list <- zico_list
  
  for (i in 1:length(zico_list)) {
    zico.obj <- zico_list[[i]]
    FDR_list[[i]] <- zico.obj$p.adj.fdr
    p_list[[i]] <- zico.obj$p.raw
    r2 <- zico.obj$R2[,'Func1']
    co <- zico.obj$coef.list[[1]][grep(grp.names[i],rownames(zico.obj$coef.list[[1]])),]
    R2_list[[i]] <-  r2 * sign(co)
    
  }
  
  unique_taxa <- unique(unlist(lapply(FDR_list, names)))
  
  df = data.frame(matrix(nrow = length(unique_taxa), ncol = length(grp.names))) 
  rownames(df) <- unique_taxa
  colnames(df) <- grp.names
  
  df[,1:length(grp.names)] <- 1
  FDR_df <- df
  p_df <- df
  df[,1:length(grp.names)] <- 0
  R2_df <- df
  
  for (i in 1:ncol(df)) {
    FDR_df[names(FDR_list[[i]]),i] <- FDR_list[[i]]
    p_df[names(p_list[[i]]),i] <- p_list[[i]]
    R2_df[names(R2_list[[i]]),i] <- R2_list[[i]]
  }
  

  sig.idx <- apply(FDR_df, 1, function(x) sum(x <= max.cutoff) > 0) # only include significant OTU with FDR < fdr.cutoff, i.e., 0.2
  FDR.df <- as.matrix(FDR_df[sig.idx,, drop =F])
  R2.df.FDR <- as.matrix(R2_df[sig.idx,, drop =F])
  
  sig.idx <- apply(p_df, 1, function(x) sum(x < max.cutoff) > 0) 
  p.df <- as.matrix(p_df[sig.idx,])
  R2.df.p <- as.matrix(R2_df[sig.idx,, drop =F])
  
  
  if (p_or_q == "p") {
    stats_df <- p.df
    R2.df <- R2.df.p
  }
  if (p_or_q == "q") {
    stats_df <- FDR.df
    R2.df <- R2.df.FDR
  }
  
  colnames(stats_df) <- grp.labels
  colnames(R2.df) <- grp.labels
  
  #remove rows that have NAs
  stats_df <- stats_df[!rowSums(is.na(stats_df)) > 0,]
  R2.df <- R2.df[rownames(stats_df),]
  
  #replace dots by spaces
  new_rownames <- gsub("^s__", "", rownames(stats_df))
  new_rownames <- gsub("_", " ", new_rownames)
  new_rownames <- gsub("\\.", " ", new_rownames)
  
  
  if (new_names==TRUE) {
    new_rownames <- as.character(new_names_vec[new_rownames])
  }
  
  rownames(stats_df) <- new_rownames
  rownames(R2.df) <- new_rownames
  
  legend_name <- paste0(p_or_q, "-cutoff")

  grid.text_1 <- '***'
  grid.text_2 <- '**'
  grid.text_3 <- '*'
  
  if (rotate == TRUE) {
    stats_df <- t(stats_df)
    R2.df <- t(R2.df)
  }
  
  lgd_sig_p = Legend(pch = c("*","**","***"), type = "points", labels = paste0("<", sign_cutoffs), background="white", title = p_or_q)

  if (single_cutoff == TRUE) {
    grid.text_1 <- '*'
    grid.text_2 <- '*'
  }
  
  #make the scale symmetric
  scale_max <- max(c(abs(min(R2.df)), max(R2.df)))
  
  heatmap_plot <- Heatmap(R2.df, name = "R2", 
          col = colorRamp2(c(-scale_max, 0, scale_max), 
                            c("#2171B5",'white', "#FEB24C")),
          column_gap = unit(1, "mm"), 
          border = TRUE,
          column_names_rot = column_name_rot,
          
          cell_fun = function(j, i, x, y, width, height, fill) {
            if(stats_df[i,j] < q.cut1) {
              grid.text(grid.text_1, x = x, y = y, hjust = 0.5, vjust=0.7)
            }else if(stats_df[i,j] < q.cut2){
              grid.text(grid.text_2, x = x, y = y, hjust = 0.5, vjust=0.7)
            }else if(stats_df[i,j] <= q.cut3){
              grid.text(grid.text_3, x = x, y = y, hjust = 0.5, vjust=0.7)
            }
          },
          rect_gp = gpar(col= "white"),
          show_column_names = T, show_row_names = T,
          show_column_dend = F, show_row_dend = F,
          show_heatmap_legend = T,
          row_names_gp = gpar(fontface = 'italic'),
          row_names_max_width = max_text_width(rownames(R2.df), gp = gpar(fontsize = 18))) 

  packed_legends <- packLegend(lgd_sig_p, direction = "vertical")
  heatmap_plot <- draw(heatmap_plot, annotation_legend_list = list(packed_legends), merge_legend=T, 
                       padding = unit(c(0.1, left_margin, 0.1, 0.1), "cm")) 

  return(list(heatmap_plot = heatmap_plot, sign_names = colnames(R2.df)))

}
