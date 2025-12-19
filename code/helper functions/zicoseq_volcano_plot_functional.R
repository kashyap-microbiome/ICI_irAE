zicoseq_volcano_plot_func <- function(zico.obj, grp.name, p.cutoff, top.n, p_value_class, log2FC=NA, update_names="No",
                                      update_df="", save_plot = FALSE) {
  
  
  # arrange the information for ZicoSeq plot
  df <- cbind.data.frame(p.raw=zico.obj$p.raw, # raw p value
                         p.adj.fdr=zico.obj$p.adj.fdr, # fdr adjusted p value
                         R2=zico.obj$R2[,'Func1'], # effect size 
                         coef = zico.obj$coef.list[[1]][grep(grp.name,rownames(zico.obj$coef.list[[1]])),] # coefficient signify the change direction
  ) %>% 
    rownames_to_column('otu') %>%
    mutate(R2.2 = sign(coef) * R2, # prepare x-axis
           `-log10(p-value)`=-log10(p.raw), # prepare y-axis
           `-log10(q-value)`=-log10(p.adj.fdr), # prepare y-axis
    )
  
  
  #-----------------------------------------------------------------------------
  #updating names / otu column
  if (class(update_df) == "character" & !"EC|KO" %in% update_names & !update_names == "pathways") 
    { print("please provide a dataframe to update names")}
  
  if (update_names == "KO") {
    rownames(update_df) <- make.unique(update_df$KO)
    rownames(df) <- df$otu
    name_vec <- update_df[rownames(df), "name_split_KO"]
    name_vec[is.na(name_vec)] <- rownames(df)[is.na(name_vec)]
    name_vec <- gsub(" NA", "", name_vec)
    df$otu <- name_vec
  }
  
  if (update_names == "EC") {
    rownames(update_df) <- make.unique(update_df$EC_number)
    rownames(df) <- df$otu
    name_vec <- update_df[rownames(df), "name"]
    name_vec[is.na(name_vec)] <- rownames(df)[is.na(name_vec)]
    name_vec <- gsub(" NA", "", name_vec)
    df$otu <- name_vec
  }
  
  #-----------------------------------------------------------------------------
  
  
  y_labs <- c('-log10(p-value)', '-log10(q-value)')
  names(y_labs) <- c("raw", "adj")
  y_lab <- as.character(y_labs[names(y_labs) == p_value_class])
  y_vals <- df[,y_lab]
  
  df$sig = ifelse(df[,y_lab] > -log10(p.cutoff),'sig','nosig') # prepare color group for dots

  sign_otu_df <- df[df$sig == "sig",]
  selected <- sign_otu_df$otu[order(sign_otu_df$R2, decreasing = T)[1:top.n]] # prepare top.n OTU names to be highlighted in the plot
  df$label <- ''
  df[df$otu %in% selected,'label'] <- df[df$otu %in% selected,'otu']
  

  size=0.8
  if (!unique(is.na(log2FC))) {
    size <- sqrt(abs(3*log2FC))
  }

  ## volcano plot
  output_plot <- ggplot(df, aes(x = R2.2, y = y_vals)) + 
    geom_point(aes(color = sig), size = size) + 
    geom_vline(aes(xintercept = 0), color = "gray") + 
    geom_hline(aes(yintercept = -log10(p.cutoff)), color = "gray") + 
    theme_bw() +
    scale_colour_manual(values = c('sig'='#c91f37','nosig'='darkgrey')) + 
    scale_y_continuous(limits = c(0,max(y_vals) * 1.1)) + 
    scale_x_continuous(limits = c(min(df$R2.2)*1.1,max(df$R2.2) * 1.1)) + 
    ggrepel::geom_text_repel(aes(label = label),max.overlaps = Inf, color = "black") + 
    theme(axis.text = element_text(color = "black", size = 12), 
          axis.title = element_text(color = "black", size = 12), 
          panel.grid = element_blank(),
          legend.position = 'none',
          legend.text = element_text(color = "black", size = 12), 
          legend.title = element_text(color = "black", size = 12),
          plot.margin = unit(c(0.5,1,0.5,0.5), "lines"))+
    labs(x = bquote(R^2~ 'x direction'), y = y_lab) +
    labs(title = paste0("functional ", update_names, " ", grp.name)) +
    labs(caption = paste0("raw or adj. p-values = ", p_value_class, "  ;  cutoff = ", p.cutoff))

  
  if (save_plot == TRUE) {
    filename <- paste0("./figures/func_volcano_", update_names, "_", grp.name ,".png")
    png(filename, width = 6*330, height = 5*330, res = 330)
    print(output_plot)
    dev.off()
  }
  
  
  return(list(output_plot = output_plot, res_df = df))
  
}