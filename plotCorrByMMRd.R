plotCorr.MMRd <- function (x, y, z) {
  my_comparisons <- list()
  df_corr <- compare_means(as.formula(paste0(y," ~ ", z)),  data = x)
  idx <- 1
  for(sig in which(df_corr$p.signif != "ns")){
    #print(as.character(df_corr[sig,c("group1")]))
    
    my_comparisons[[idx]] <- c(as.character(df_corr[sig,c("group1")]), as.character(df_corr[sig,c("group2")]))
    idx <- idx + 1
  }
  n <- length(my_comparisons)
  write.table(df_corr, 
              file = file.path(getwd(),paste0(y,".",z,".corr.result.txt")), 
              quote = FALSE, 
              sep = "\t")
  
  ylabel <- y
  
  if (y == "Fraction_Genome_Altered") {
    max_y <- max(x[,y])
    min_y <- min(x[,y])
    assign(paste0("p_",y), ggboxplot(x, 
                                     x = z, y = y,
                                     #color = "Cancer_Type_NEW", 
                                     palette = "category20", 
                                     size = 0.5,
                                     add = c("mean","mean_sd"),
                                     add.params = list(color = "red")) +
             rotate_x_text(angle = 25) + 
             ylim(min_y, 2*max_y) + 
             ylab(ylabel) + 
             xlab("MMRd Type") +
             stat_compare_means(aes(label = ..p.signif..), comparisons = my_comparisons, label.y = 1.6*max_y*(1-0.02*(1:n)), vjust = 0.8, bracket.size = 0.1) + 
             stat_compare_means(label.y = 1.8*max_y) + 
             theme(axis.text = element_text(size = 10)) + 
             theme(axis.title = element_text(face = "bold")) + 
             theme(legend.text = element_text(size = 7)) + 
             theme(plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"))
    )
    return(get(paste0("p_",y)))    
  }else {
    max_y <- log10(max(x[,y]))
    min_y <- log10(min(x[,y]))
    assign(paste0("p_",y), ggboxplot(x, 
                                     x = z, y = y,
                                     #color = "Cancer_Type_NEW", 
                                     palette = "category20", 
                                     size = 0.5,
                                     add = c("mean","mean_sd"),
                                     add.params = list(color = "red")) + 
             #add = "jitter") + 
             rotate_x_text(angle = 25) + 
             ylim(min_y, 2*max_y) + 
             ylab(ylabel) + 
             xlab("MMRd Type") +
             stat_compare_means(aes(label = ..p.signif..), comparisons = my_comparisons, label.y = 1.6*max_y*(1-0.03*(1:n)), vjust = 0.8, bracket.size = 0.1) + 
             stat_compare_means(label.y = 2*max_y) + 
             scale_y_log10() +
             theme(axis.text = element_text(size = 10)) + 
             theme(axis.title = element_text(face = "bold")) + 
             theme(legend.text = element_text(size = 7)) + 
             theme(plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"))
    )
    return(get(paste0("p_",y)))
  }
  
}
