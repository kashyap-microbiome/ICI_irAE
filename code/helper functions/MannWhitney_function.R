MannWhitney_function <- function(tax_df, groups, meta) {

  col_ind <- which(colnames(meta) == groups)
  category <- as.factor(meta[,col_ind])
  
  if(length(levels(category)) == 2){
    ix1 <- category == levels(category)[1]
    ix2 <- category == levels(category)[2]
    pvals <- apply(tax_df,1,function(x) wilcox.test((x[ix1]**(1/2)), (x[ix2]**(1/2)),exact=FALSE)$p.value)
    stats <- apply(tax_df,1,function(x) wilcox.test((x[ix1]**(1/2)), (x[ix2]**(1/2)),exact=FALSE)$statistic)

    # classwise means
    classwise.means <- t(apply(tax_df,1,function(xx) sapply(split(xx,category),mean)))
    colnames(classwise.means) <- paste0(groups, colnames(classwise.means))
    
    log2FC <- log2(classwise.means[,1] / classwise.means[,2])
    
  }

  adj.pvals <- rep(NA,length(pvals))
  names(adj.pvals) <- names(pvals)
  na.ix <- is.na(pvals)
  adj.pvals[!na.ix] <- p.adjust(pvals[!na.ix],'fdr')

  res_list <- list(pvals = pvals, stats =  stats, adj.pvals = adj.pvals, classwise.means = classwise.means, log2FC = log2FC)
  res_df <- as.data.frame(do.call(cbind, res_list))
  
  return(res_df)
  
}