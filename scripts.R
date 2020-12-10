###########################################################################
# plot case frequency by MMRd and Cancer types
###########################################################################
plotCaseFreq <- function(x, y, z) {
  # By Cancer types
  p_freq_abs <- ggplot(x) +
    geom_bar(aes_string(x=y,fill=z)) +
    scale_fill_brewer(palette = "Dark2") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    xlab("Cancer Type") + labs(fill="MMRd Type")
  print(p_freq_abs)
  
  p_freq_pct <- ggplot(x) +
    geom_bar(aes_string(x=y,fill=z), position = "fill") +
    ylab("Percentage") +
    xlab("Cancer Type") + labs(fill="MMRd Type") +
    scale_fill_brewer(palette = "Dark2") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
  print(p_freq_pct)
  
  # By MMRd types
  library("ggsci")
  p_freq_abs.1 <- ggplot(x) +
    geom_bar(aes_string(x=z,fill=y)) +
    scale_fill_d3(palette = "category20") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    xlab("MMRd Type") + labs(fill="Cancer Type")
  print(p_freq_abs.1)
  
  p_freq_pct.1 <- ggplot(x) +
    geom_bar(aes_string(x=z,fill=y), position = "fill") +
    ylab("Percentage") +
    xlab("MMRd Type") + labs(fill="Cancer Type") +
    scale_fill_d3(palette = "category20") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
  print(p_freq_pct.1)
}

###########################################################################
# plot survival plots by MMRd and Cancer types and key oncogenic gene mutaions
###########################################################################


drawPairwiseSurvivalPlot <- function(data, group, ...) {
  #subset by cancer type if available
  args <- list(...)
  cancerType = ""

  if(length(args) == 1){
    cancerType = args[[length(args)]]
    #print(cancerType)
    data <- data[data$Cancer_Type == cancerType,]
  }
  #print(dim(data))
  #print(cancerType)
  # factorize the status
  data[,"status"] <- as.numeric(as.factor(data[,"status"]))
  
  # factorize the group
  data[,group] <- as.factor(data[,group])
  assign(paste0(group,"_lvls"), levels(data[,group]))
  data[,group] <- as.numeric(data[,group])
    
  # create pairs from all the levels
  d_cmp <- expand.grid(x = unique(data[,group]), y = unique(data[,group]), stringsAsFactors = F) 
  d_cmp <- d_cmp[d_cmp$x < d_cmp$y,]
  
  # generate survial plots in pairwise mode
  for (i in 1:nrow(d_cmp)) {

    subset <- which(data[,group] %in% as.numeric(d_cmp[i,]))

    df_pair <- data[subset,] 
    df_pair[,group] <- factor(df_pair[,group], levels = unique(df_pair[,group]))

    modform <- as.formula(paste0("Surv(time = time, event =  status) ~ ", group))

    fit <- eval(substitute(survfit(modform, data = df_pair), list(modform = modform)))
    sdf <- eval(substitute(survdiff(modform, data = df_pair), list(modform = modform)))
    p.Val <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
    #print(p.Val)
    
    if(p.Val < 0.1) {
      # Estimate hazard ratio pairwisely
      res.cox <- eval(substitute(coxph(formula = modform, data = df_pair), list(modform = modform)))
      
      #summary(res.cox) # summarize the result
      # collect the following values: HR, lower 95%, upper 95%
      estimate = exp(res.cox$coefficients)
      conf.low = exp(confint(res.cox)[1,1])
      conf.high = exp(confint(res.cox)[1,2])
      
      # Drawing curves
      
      #p <- ggsurvplot(fit, pval = T, pval.method = T, legend.labs = MMRd_lvls_1[as.numeric(d_cmp[i,])], 
      #legend.labs = paste0(group,"_lvls")[as.numeric(d_cmp[i,])],
      p <- ggsurvplot(fit, pval = T, pval.method = T, combine = T, legend.title = group, legend.labs = get(paste0(group,"_lvls"))[as.numeric(d_cmp[i,])],
                      conf.int = T,  risk.table = T, risk.table.height = 0.3, risk.table.col = "strata", 
                      title = paste("All",cancerType,"Cancer Samples", sep = " "),
                      risk.table.y.text.col = T, data = df_pair)
      p$plot <- p$plot + 
        annotate("text", x = Inf, y = Inf, vjust = 1, hjust = 1, 
                 label = paste0("HR = ", round(estimate,2), " ( ",round(conf.low,2), " - ", round(conf.high,2), " ) "),
                 size = 4)
      print(p)
      #assign(paste0("p",i),  p$plot)
      
    }
    
  }
  
}

###########################################################################
###########################################################################
plotCorr.Cancer <- function(x, y) {
  my_comparisons <- list()
  df_corr <- compare_means(as.formula(paste0(y," ~ ", "Cancer_Type_NEW")),  data = x)
  idx <- 1
  for(sig in which(df_corr$p.signif != "ns")){
    #print(as.character(df_corr[sig,c("group1")]))
    
    my_comparisons[[idx]] <- c(as.character(df_corr[sig,c("group1")]), as.character(df_corr[sig,c("group2")]))
    idx <- idx + 1
  }
  n <- length(my_comparisons)
  write.table(df_corr, 
              file = file.path(getwd(),paste0(y,".corr.result.by.Cancer.txt")), 
              quote = FALSE, 
              sep = "\t")
  
  ylabel <- y
  
  if (y == "Fraction_Genome_Altered") {
    max_y <- max(x[,y])
    min_y <- min(x[,y])
    assign(paste0("p_",y), ggboxplot(x, 
                                     x = "Cancer_Type_NEW", y = y,
                                     palette = "category20", 
                                     size = 0.5,
                                     add = c("mean","mean_sd"),
                                     add.params = list(color = "red")
                                      ) + 
             rotate_x_text(angle = 25) + 
             ylim(min_y, 2*max_y) + 
             ylab(ylabel) + 
             xlab("Cancer Type") +
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
                                     x = "Cancer_Type_NEW", y = y,
                                     palette = "category20", 
                                     size = 0.5,
                                     add = c("mean","mean_sd"),
                                     add.params = list(color = "red")) + 
             rotate_x_text(angle = 25) + 
             ylim(min_y, 2*max_y) + 
             ylab(ylabel) + 
             xlab("Cancer Type") +
             stat_compare_means(aes(label = ..p.signif..), comparisons = my_comparisons, label.y = 1.8*max_y*(1-0.02*(1:n)), vjust = 0.8, bracket.size = 0.1) + 
             stat_compare_means(label.y = 2*max_y) + 
             scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                           labels = trans_format("log10", math_format(10^.x))) +
             theme(axis.text = element_text(size = 10)) + 
             theme(axis.title = element_text(face = "bold")) + 
             theme(legend.text = element_text(size = 7)) + 
             theme(plot.margin = unit(c(0.5,0.5,0.5,0.5),"cm"))
    )
    return(get(paste0("p_",y)))
  }
  
}


########################################################
# Generate Signature mutations (SBS&ID) bubble plot by MMRd types
########################################################
drawSigMutBubblePlot <- function (files, sigs, sigType) {
  df.samples <- data.frame()
  
  for (file in files) {
    df <- read.delim2(file = file, quote = "", stringsAsFactors = FALSE)
    row.names(df) <- df$Samples
    # remove the "Samples" column
    df <- df[,! names(df) %in% c("Samples","X96A","IDA")]
    
    #df_pad <- data.frame()
    for (sig in sigs[! sigs %in% colnames(df)]) {
      df <- cbind(df, new = 0)
      names(df)[dim(df)[2]] <- sig
    }
    #tmp <- cbind(df, df_pad)
    df <- df[,order(colnames(df))]
    # concatenate to the final data frame
    df.samples <- rbind(df.samples, df)
  }
  
  
  # MSK-IMPACT target size : ~ 1.5 Mb in human genome
  df_MMRd <- read.table(file = "MMRd_type.txt", header = T, sep = "\t")
  row.names(df_MMRd) <- df_MMRd$caseID
  df_MMRd <- df_MMRd[, ! names(df_MMRd) %in% c("X")]
  
  tumor.number <- as.data.frame(as.matrix(table(df_MMRd$MMRd_type)))
  colnames(tumor.number) <- c("Tumor.number")
  
  df.samples.MMRd <- merge(df.samples, df_MMRd, by = 0)
  row.names(df.samples.MMRd) <- df.samples.MMRd$Row.names
  df.samples.MMRd <- df.samples.MMRd[,! names(df.samples.MMRd) %in% c("Row.names","caseID")]
  
  # write to file
  write.table(df.samples.MMRd, file = file.path(getwd(),paste0(sigType,"_mut_No_in_sample_by_MMRd.txt")), quote = F, sep = "\t", col.names = NA, row.names = T)
  
  # estimate the number of samples relevant to any Cosmic signature
  df.samples.MMRd$SUM <- apply(df.samples.MMRd[,!names(df.samples.MMRd) %in% c("MMRd_type")],1,sum)
  
  sig.tumor.count <- df.samples.MMRd %>% group_by(MMRd_type) %>% summarise(count = sum(SUM>0, na.rm = TRUE))
  # aggregate median # of mutation due to the SBS signatures by MMRd types
  median_mut <- aggregate(.~MMRd_type, 
                          df.samples.MMRd[,!names(df.samples.MMRd) %in% c("SUM")],
                          FUN=(function(x){
                            ifelse(sum(x==0)>0 & sum(x !=0) >0, median(x[x>0]), median(x))
                          })
  )
  row.names(median_mut) <- median_mut$MMRd_type
  #median_mut <- median_mut[,!names(median_mut) %in% c("MMRd_type")]
  #colnames(median_mut) <- paste0(names(median_mut),".median")
  # normalize by the size of MSK-IMPACT-468genes size 1.5Mb
  median_mut <- log2(median_mut[,!names(median_mut) %in% c("MMRd_type")]/1.5)
  median_mut$MMRd_type <- row.names(median_mut)
  
  median_mut_melt <- melt(median_mut, 
                          variable.name = sigType,
                          value.name = "Median_mut",
                          id.vars = c("MMRd_type"))
  
  
  # dplyr way
  median_mut_dplyr <- df.samples.MMRd %>% 
    group_by(MMRd_type) %>% 
    summarise_each(funs(median))
  
  # aggregate mean # of mutation due to the SBS signatures by MMRd types
  mean_mut <- aggregate(.~MMRd_type, 
                        df.samples.MMRd[,!names(df.samples.MMRd) %in% c("SUM")],
                        FUN=(function(x){
                          ifelse(sum(x==0)>0 & sum(x !=0) >0, mean(x[x>0]), mean(x))
                        })
  )
  
  row.names(mean_mut) <- mean_mut$MMRd_type
  #mean_mut <- mean_mut[,!names(mean_mut) %in% c("MMRd_type")]
  #colnames(mean_mut) <- paste0(names(mean_mut),".mean")
  mean_mut <- log2(mean_mut[,!names(mean_mut) %in% c("MMRd_type")]/1.5)
  mean_mut$MMRd_type <- row.names(mean_mut)
  
  
  mean_mut_melt <- melt(mean_mut, 
                        variable.name = sigType,
                        value.name = "Mean_mut",
                        id.vars = c("MMRd_type"))
  # function count portion of tumors with a signature
  portion <- function(x) {
    result <- sum(x>0)/length(x)
    return(result)
  }
  
  portion_tumor <- aggregate(.~MMRd_type, 
                             df.samples.MMRd[,!names(df.samples.MMRd) %in% c("SUM")],
                             portion)
  row.names(portion_tumor) <- portion_tumor$MMRd_type
  
  portion_tumor_melt <- melt(portion_tumor, 
                             variable.name = sigType,
                             value.name = "Tumor_portion",
                             id.vars = c("MMRd_type"))
  
  
  ##############################################################################
  # Combine all the summary data, mean/meadian of Mut #/Mb, portion.tumor, tumor.number
  ##############################################################################
  #df_plot <- expand.grid(as.factor(row.names(portion_tumor)), as.factor(sigs))
  
  df_plot <- merge(mean_mut_melt, median_mut_melt, by = c("MMRd_type",sigType))
  df_plot <- merge(df_plot, portion_tumor_melt, by = c("MMRd_type",sigType))
  #levels(df_plot$MMRd_type) <- sigs # reorder the levels
  #df_plot <- df_plot[order(df_plot$SBS),]
  
  options(repr.plot.width = 14, repr.plot.height = 8)
  
  df_plot %>% mutate(SBS = factor(SBS, levels=sigs)) %>%
    ggplot(aes(x=MMRd_type, y=SBS, size = Tumor_portion, color = Median_mut)) + 
    geom_point(alpha=0.7) +
    scale_size(range = c(.1, 10), name="Proportion of tumours \nwith the signature") + 
    scale_fill_viridis(discrete=TRUE, guide=FALSE, option="B") + 
    scale_colour_viridis_c("log2 Median mutations \n per Mb due to signature \n(among tumours with \nthe signature)") +
    theme(plot.margin = margin(16,16,16,16)) + 
    theme(legend.position="bottom", legend.justification = "left", legend.margin = margin(6,6,6,6)) +
    scale_x_discrete(labels=c("Complex" = "Complex\n8/11", "MLH1_PMS2" = "MLH1_PMS2\n164/246","MSH2_MSH6" = "MSH2_MSH6\n68/69",
                              "MSH6_only" = "MSH6_only\n12/12", "NORMAL_only" = "NORMAL_only\n17/19", "PMS2_only" = "PMS2_only\n3/13", "UNK" = "UNK\n9/11")) + 
    theme(axis.text.x = element_text(size = 8)) +
    theme(legend.title = element_text(size = 8, face = "bold"),
          legend.text = element_text(size = 8)) + 
    ylab("SBS") +
    xlab("MMRd Types") +
    theme(axis.title = element_text(size = 15, face = "bold"))
}