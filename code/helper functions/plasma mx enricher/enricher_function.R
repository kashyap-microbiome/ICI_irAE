enricher_function <- function(names_sig, level="pathway_name", matched_df, pvalueCutoff=0.25, qvalueCutoff = 1) {
  
  #Ruben Mars, 202501
  
  require(clusterProfiler)
  require(pheatmap)
  require(DOSE)
  require(enrichplot)
  require(ggupset)
  
  
  term2gene_df <- as.data.frame(cbind(term=matched_df[,level], gene=rownames(matched_df)))
  
  #packageVersion("enrichplot") #1.26.2
  

  enrichres <- enricher(
    names_sig,
    pvalueCutoff = pvalueCutoff,
    pAdjustMethod = "fdr",
    universe = NULL,
    minGSSize = 2,
    maxGSSize = 500,
    qvalueCutoff = qvalueCutoff,
    gson = NULL,
    TERM2GENE = term2gene_df
  )
  
  res_df <- as.data.frame(enrichres)
  if (nrow(res_df) > 0) { 
    res_barplot_n <- barplot(enrichres, showCategory = 20) 
    res_barplot_q <- mutate(enrichres, qscore = -log(p.adjust, base = 10)) %>% 
      barplot(x = "qscore")
    res_dotplot <- dotplot(enrichres, showCategory = 15) + ggtitle("")
    
    
    res_networkplot <- cnetplot(enrichres, layout = "kk")
    
    res_heatmap8 <- heatplot(enrichres, showCategory = 8)
  } else {
    res_barplot_n <- NULL
    res_barplot_q <- NULL
    res_dotplot <- NULL
    res_networkplot <- NULL
    res_heatmap8 <- NULL
  }
    
  
  return(list(enrichres=enrichres, res_df=res_df, res_barplot_n=res_barplot_n, res_barplot_q=res_barplot_q, res_dotplot=res_dotplot,
         res_networkplot=res_networkplot, res_heatmap8=res_heatmap8))
  
}