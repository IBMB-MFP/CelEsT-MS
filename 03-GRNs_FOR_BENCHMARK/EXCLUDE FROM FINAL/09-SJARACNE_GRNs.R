# read in SJARACNe output file
SJARACNeage_out <- read.table("output/SJARACNe/age_corrected_non_zero/consensus_network_ncol_.txt",
                              header = TRUE)

# filter by p < 1e-5

SJARACNEage_out_filt <- SJARACNeage_out[SJARACNeage_out$p.value < 0.00001, ]

SJARACNEage_out_filt_separated <- lapply(unique(SJARACNEage_out_filt$source), function(x){
  
  thisone_SJ <- SJARACNEage_out_filt[SJARACNEage_out_filt$source == x, ]
  
  thisone_SJ[order(thisone_SJ$MI, decreasing = TRUE), ]
  
})

names(SJARACNEage_out_filt_separated) <- unique(SJARACNEage_out_filt$source)

View(SJARACNEage_out_filt_separated[[1]])
sum(SJARACNEage_out_filt_separated[[1]]$slope > 0)

hist(SJARACNEage_out_filt_separated[[1]]$pearson)

SJARACNEage_GRN_list <- lapply(SJARACNEage_out_filt_separated, function(x){
  
  temp_df <- x[, c("source", "target")]

  temp_df[, "weight"] <- sign(x$pearson)
  
  temp_df
  
})

SJARACNEage_GRN <- do.call(rbind, SJARACNEage_GRN_list)

write.table(SJARACNEage_GRN,
            "output/GRNs/SJARACNEage_GRN.txt",
            col.names = TRUE,
            row.names = FALSE,
            sep = "\t")

