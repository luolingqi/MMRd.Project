
drawSigMutDotplotByMMRd <- function(file_pattern, sigs, sigType){
  # read in all the Decomposed_Solution_Activities_SBS96 files
  # padding with all the cosmic SBS signatures (add 0 value)
  # concatenate all the dataframe into one
  
  
  files <- list.files(path = getwd(), 
                      pattern = file_pattern,
                      #pattern = 'Decomposed_Solution_Activities_SBS96.txt', 
                      full.names = TRUE, recursive = TRUE)
  
  
  #sigs <- c("SBS1","SBS2","SBS4","SBS5","SBS6","SBS10a","SBS10b","SBS13","SBS14","SBS15","SBS17a","SBS17b","SBS20","SBS21","SBS24","SBS26","SBS33","SBS42","SBS44","SBS45","SBS46","SBS53","SBS54","SBS59")
  
  #sigs <- c("SBS1","SBS5","SBS6","SBS10b","SBS14","SBS15","SBS20","SBS21","SBS26","SBS44")
  
  # Prepare the a grand data frame of signatures for all the samples
  df.samples <- data.frame()
  for (file in files) {
    df <- read.delim2(file = file, quote = "", stringsAsFactors = FALSE)
    row.names(df) <- df$Samples
    
    #Subset for key signatures only
    df <- df[, names(df) %in% sigs]
    
    # padding with 0 for those signatures not in the data
    for (sig in sigs[! sigs %in% colnames(df)]) {
      df <- cbind(df, new = 0)
      names(df)[dim(df)[2]] <- sig
    }
    
    df <- df[,order(colnames(df))]
    # concatenate to the final data frame
    df.samples <- rbind(df.samples, df)
  }
  
  # merge in the MMRd type information
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
  write.table(df.samples.MMRd, 
              file = file.path(getwd(),paste0(sigType,"_mut_No_in_sample_by_MMRd.txt")), 
              quote = F, sep = "\t", col.names = NA, row.names = T)
  
  
  
  # estimate the number of samples relevant to any Cosmic signature
  df.samples.MMRd$SUM <- apply(df.samples.MMRd[,!names(df.samples.MMRd) %in% c("MMRd_type")],1,sum)
  
  sig.tumor.count <- df.samples.MMRd %>% group_by(MMRd_type) %>% summarise(count = sum(SUM>0, na.rm = TRUE))
  # aggregate median # of mutation due to the SBS signatures by MMRd types
  median_mut <- aggregate(.~MMRd_type, 
                          df.samples.MMRd[,!names(df.samples.MMRd) %in% c("SUM")],
                          FUN=(function(x){
                            # if value set contains both 0 and not 0, get median from non_0 values only
                            ifelse(sum(x==0)>0 & sum(x !=0) >0, median(x[x>0]), median(x))
                          })
  )
  row.names(median_mut) <- median_mut$MMRd_type
  
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
  
  df_plot <- merge(mean_mut_melt, median_mut_melt, by = c("MMRd_type",sigType))
  df_plot <- merge(df_plot, portion_tumor_melt, by = c("MMRd_type",sigType))
  # Remove "UNK" MMRd Type
  df_plot <- df_plot[df_plot$MMRd_type != "UNK",]
  df_plot[[sigType]] <- factor(df_plot[,sigType], levels = sigs)
  #levels(df_plot$MMRd_type) <- sigs # reorder the levels
  #df_plot <- df_plot[order(df_plot$SBS),]
  
  options(repr.plot.width = 14, repr.plot.height = 8)
  
  #Prepare labels
  for(type in unique(df_plot$MMRd_type)){
    assign(paste("n",type,sep = "."), as.character(sig.tumor.count[sig.tumor.count$MMRd_type==type,"count"]))
  }
  
  
  
  # Plot
  #df_plot %>% mutate(SBS = factor(SBS, levels=sigs)) %>%
  #print(head(df_plot))
  
  #p <- df_plot %>% mutate("{sigType}" := factor("{sigType}", levels = sigs)) %>%
  p <- df_plot %>%
    ggplot(aes_string(x="MMRd_type", y=sigType, size = "Tumor_portion", color = "Median_mut")) + 
    geom_point(alpha=0.7) +
    scale_size(range = c(.1, 10), name="Proportion of tumours \nwith the signature") + 
    scale_fill_viridis(discrete=TRUE, guide=FALSE, option="B") + 
    scale_colour_viridis_c("log2 Median mutations \n per Mb due to signature \n(among tumours with \nthe signature)") +
    theme(plot.margin = margin(16,16,16,16)) + 
    theme(legend.position="bottom", legend.justification = "left", legend.margin = margin(6,6,6,6)) +
    scale_x_discrete(labels=c("Complex" = paste0("Complex\n",n.Complex,"/11"), 
                              "MLH1_PMS2" = paste0("MLH1_PMS2\n",n.MLH1_PMS2,"/246"),
                              "MSH2_MSH6" = paste0("MSH2_MSH6\n",n.MSH2_MSH6,"/69"),
                              "MSH6_only" = paste0("MSH6_only\n",n.MSH6_only,"/12"), 
                              "NORMAL_only" = paste0("NORMAL_only\n",n.NORMAL_only,"/19"), 
                              "PMS2_only" = paste0("PMS2_only\n",n.PMS2_only,"/13"))) + 
    theme(axis.text.x = element_text(size = 8)) +
    theme(legend.title = element_text(size = 8, face = "bold"),
          legend.text = element_text(size = 8)) + 
    ylab(sigType) +
    xlab("MMRd Types") +
    theme(axis.title = element_text(size = 15, face = "bold"))
  
  print(p)
}
