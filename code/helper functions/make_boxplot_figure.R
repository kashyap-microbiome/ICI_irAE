make_boxplot_figure <- function(x, groups, project="exploring", plottitle="test", xlable="", ylable="value", 
                                plot_width=4, plot_height=4, full_groups=TRUE, horizontal=F, 
                                left_margin=3.5, bottom_margin=2.5, file_type="png", add_significance=FALSE, 
                                q.vals="", q_or_p = "q", add_significance_for_5=FALSE, symbol=F) {

  sel_colors <- adjustcolor(sel_colors_full_names, alpha.f = 0.8)
  groups <- gsub("_", "-", groups)
  
  sign_symbols <- c("***", "**", "*", "^", "^^")
  sign_cutoffs <- c(0.001, 0.01, 0.05, 0.102, 0.2) #.102 done for one specific plot
  
  n_groups <- length(unique(groups))
  if (n_groups == 5) {
    sel_colors <- c(sel_colors[1],sel_colors[2],sel_colors[2],sel_colors[3],sel_colors[3])
  }
  
  x <- as.numeric(x)
  #might be necessary later
  #data$names=factor(data$names , levels=levels(data$names)[c(1,4,3,2)])
  
  if (file_type == "png") {
    png(paste(project, "_", plottitle, ".png", sep=""), width=plot_width*250, height=plot_height*250, res=300)    
  }
  if (file_type == "pdf") {
    pdf(paste(project, "_", plottitle, ".pdf", sep=""), width=plot_width*0.8, height=plot_height*0.8)
  }
  par(mar=(c(bottom_margin,left_margin,2,1)))
  boxplot(x ~ groups, varwidth=T, las=1, xlab="", ylab="", outline=F, ylim=c(min(x), max(x)*1.3), border=sel_colors, 
          horizontal=horizontal, col="white")
  mtext(xlable, 1, line=2.7) 
  mtext(ylable, 2, line=left_margin-1)
  title(plottitle, line=0.5, cex.main=0.8)
  stripchart(x ~ groups, vertical = !horizontal, method = "jitter", jitter= 0.25, add = TRUE, pch=16, col=sel_colors, cex=0.8)

  sign_type <- "q="
  if (q_or_p == "p") {
    sign_type <- "P="
  }
  
  #add the significance indications, will be plotted if below 0.25
  if (add_significance == TRUE) {
    parusr <- par('usr')
    if (length(q.vals) == length(unique(groups))-1) {
    #order of q.vals vector is H-C, H-D
      q.vals <- round(q.vals, 3)
      draw_lines <- which(q.vals <= 0.25)
      if (1 %in% draw_lines) {
        scale_temp <- 1.17
        segments(1, parusr[4]/scale_temp, 2, parusr[4]/scale_temp) 
        scale_temp <- 1.15
        if (symbol==F) {
          text(1.5, parusr[4]/scale_temp, paste(sign_type, q.vals[1], sep=""), cex=0.8, font=3) }
        if (symbol==T) {
          temp_symb <- sign_symbols[min(which(q.vals[1] < sign_cutoffs))]
          text(1.5, parusr[4]/scale_temp, paste(temp_symb, sep=""), cex=1.5, font=3) }
      }
      if (2 %in% draw_lines) {
        scale_temp <- 1.13
        segments(1, parusr[4]/scale_temp, 3, parusr[4]/scale_temp)
        scale_temp <- 1.1
        if (symbol==F) {
          text(2, parusr[4]/scale_temp, paste(sign_type, q.vals[2], sep=""), cex=0.8, font=3) }
        if (symbol==T) {
          temp_symb <- sign_symbols[min(which(q.vals[2] < sign_cutoffs))]
          text(2, parusr[4]/scale_temp, paste(temp_symb, sep=""), cex=1.5, font=3) }
      }
    }
    #or when adding 3 q-values, also for C vs D
    if (length(q.vals) == length(unique(groups))) { 
      #order of q.vals vector is H-C, H-D, C-D
      q.vals_mod <- round(q.vals, 3)
      ind_temp <- which(q.vals < 0.001)
      if (length(ind_temp) > 0) {
        q.vals_mod[ind_temp] <- "<0.001"
      }
      draw_lines <- which(q.vals <= 0.2)
      if (1 %in% draw_lines) {
        scale_temp <- 1.17
        segments(1, parusr[4]/scale_temp, 2, parusr[4]/scale_temp) 
        scale_temp <- 1.15
        if (symbol==F) {
          text(1.5, parusr[4]/scale_temp, paste(sign_type, q.vals_mod[1], sep=""), cex=0.8, font=3) }
        if (symbol==T) {
          temp_symb <- sign_symbols[min(which(q.vals_mod[1] < sign_cutoffs))]
          text(1.5, parusr[4]/scale_temp, paste(temp_symb, sep=""), cex=1.5, font=3) }
      }
      if (2 %in% draw_lines) {
        scale_temp <- 1.12
        segments(1, parusr[4]/scale_temp, 3, parusr[4]/scale_temp)
        scale_temp <- 1.1
        if (symbol==F) {
          text(2, parusr[4]/scale_temp, paste(sign_type, q.vals_mod[2], sep=""), cex=0.8, font=3) }
        if (symbol==T) {
          temp_symb <- sign_symbols[min(which(q.vals_mod[2] < sign_cutoffs))]
          text(2, parusr[4]/scale_temp, paste(temp_symb, sep=""), cex=1.5, font=3) }
      }
      if (3 %in% draw_lines) {
        scale_temp <- 1.07
        segments(2, parusr[4]/scale_temp, 3, parusr[4]/scale_temp)
        scale_temp <- 1.05
        if (symbol==F) {
          text(2.5, parusr[4]/scale_temp, paste(sign_type, q.vals_mod[3], sep=""), cex=0.8, font=3) }
        if (symbol==T) {
          temp_symb <- sign_symbols[min(which(q.vals_mod[3] < sign_cutoffs))]
          text(2.5, parusr[4]/scale_temp, paste(temp_symb, sep=""), cex=1.5, font=3) }
      }
    }
  }

  
  if (add_significance_for_5 == TRUE) {
      parusr <- par('usr')
      #for now only works for 2 q-values which are D normal vs D flare and C normal vs C flare
      q.vals <- round(q.vals, 3)
      ind_temp <- which(q.vals < 0.001)
      if (length(ind_temp) > 0) {
        q.vals[ind_temp] <- "<0.001"
      }
      draw_lines <- which(q.vals <= 0.1)
      if (1 %in% draw_lines) {
        scale_temp <- 1.18
        segments(parusr[4]/scale_temp, 4, parusr[4]/scale_temp, 5) 
        scale_temp <- 1.12
        if (symbol==F) {
          text(parusr[4]/scale_temp, 4.5, paste(sign_type, q.vals[1], sep=""), cex=0.8, font=3, srt=270) }
        if (symbol==T) {
          temp_symb <- sign_symbols[min(which(q.vals[1] < sign_cutoffs))]
          text(parusr[4]/scale_temp, 4.5, paste(temp_symb, sep=""), cex=1.5, font=3, srt=270) }
      }
      if (2 %in% draw_lines) {
       scale_temp <- 1.18
       segments(parusr[4]/scale_temp, 2, parusr[4]/scale_temp, 3)
       scale_temp <- 1.12
       if (symbol==F) {
        text(parusr[4]/scale_temp, 2.5, paste(sign_type, q.vals[2], sep=""), cex=0.8, font=3 , srt=270) }
       if (symbol==T) {
         temp_symb <- sign_symbols[min(which(q.vals[2] < sign_cutoffs))]
         text(parusr[4]/scale_temp, 2.5, paste(temp_symb, sep=""), cex=1.5, font=3, srt=270) }
      }
  }

  dev.off()  
}

