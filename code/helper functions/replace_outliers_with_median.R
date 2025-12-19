replace_outliers_with_median <- function(df, groups, n_sd=4, verbose=FALSE) {
  #metabolites assumed to be in rows, groups in columns
  
  #get the medians per group per metabolite
  median_values <- list()
  df_output <- df
  for (i in 1:nrow(df)) {
  x <- as.numeric(df[i,])
  median_values[[i]] <- aggregate(x[!is.na(x)], by=list(groups[!is.na(x)]), FUN=median)
  outlier_index <- which((x - mean(x)) / sd(x) > n_sd)
  if (length(outlier_index) > 0) {
    if  (verbose==FALSE) {
    cat("number of removed outliers out of:\n")
    print(length(outlier_index))
    print(length(x))
    }
    for (j in 1:length(outlier_index)) {
      na_group <- groups[outlier_index[j]]
      group_index <- which(groups == na_group)
      replace_val <- median_values[[i]]$x[median_values[[i]]$Group.1 == na_group]
      df_output[i,outlier_index[j]] <- replace_val
    }
  }
  }
  return(df_output)
}