library(stringr)
library(ggplot2)

setwd("~/Cel_GRN_manuscript/")
STREMEtopmotif_orthology_probs <- readRDS("output/STREMEtopmotif_orthology_probs10k.rds")

allTFS_manualtoptargets_HOTexcl <- readRDS(file = "output/MODernENCODE_manualtoptargets_operon&HOTexcluded.rds")
allTFS_manualtoptargets <- readRDS(file = "output/MODernENCODE_manualtoptargets_operonexcluded.rds")


allTFs_labelused <- readRDS("output/allTFs_labelused.rds")

allmodERN_TFsONLY_BM <- readRDS("output/allmodERN_TFsONLY_BM.rds")

TFs <- names(allTFs_labelused)[match(allTFs_labelused[allTFs_labelused != "cfi-1"], names(allTFS_manualtoptargets_HOTexcl))]

neworder_targets_Hotexcl <- lapply(TFs, function(x){
print(x)
  temp <- allTFS_manualtoptargets_HOTexcl[[allTFs_labelused[x]]]
  
  if(!is.null(STREMEtopmotif_orthology_probs[[x]][match(temp, names(STREMEtopmotif_orthology_probs[[x]]))])){
  
  tempout <- data.frame("target" = temp, "probs" = 1) 
  tempout[, "probs"] <- STREMEtopmotif_orthology_probs[[x]][match(temp, names(STREMEtopmotif_orthology_probs[[x]]))]
  
  tempout[is.na(tempout$probs), "probs"] <- 1
  
  tempout <- tempout[order(tempout$probs), ]
  
  # order scores of 1 according to ChIP score. although it seems this order is retained, anyway
  
  tempout[tempout$probs == 1, "target"] <- temp[temp %in% tempout[tempout$probs == 1, "target"]]
  
  return(tempout$target)
  
} else {
    
  return(temp)

  }
  
})

names(neworder_targets_Hotexcl) <- TFs

neworder_targets <- lapply(TFs, function(x){
  print(x)
  temp <- allTFS_manualtoptargets[[allTFs_labelused[x]]]
  
  if(!is.null(STREMEtopmotif_orthology_probs[[x]][match(temp, names(STREMEtopmotif_orthology_probs[[x]]))])){
    
    tempout <- data.frame("target" = temp, "probs" = 1) 
    tempout[, "probs"] <- STREMEtopmotif_orthology_probs[[x]][match(temp, names(STREMEtopmotif_orthology_probs[[x]]))]
    
    tempout[is.na(tempout$probs), "probs"] <- 1
    
    tempout <- tempout[order(tempout$probs), ]
    
    # order scores of 1 according to ChIP score. although it seems this order is retained, anyway
    
    tempout[tempout$probs == 1, "target"] <- temp[temp %in% tempout[tempout$probs == 1, "target"]]
    
    return(tempout$target)
    
  } else {
    
    return(temp)
    
  }
  
})

names(neworder_targets) <- TFs

make.ChIP.GRN <- function(cutoff = 1000,
                          targets,
                          suffix = NULL,
                          prefix = "myGRN",
                          min_targets = NULL){
  
  tempGRN <- do.call(rbind, lapply(names(targets), function(thisTF){
    
    x <- targets[[thisTF]]
    thisTF_gseq <- allmodERN_TFsONLY_BM[apply(allmodERN_TFsONLY_BM, 1, function(y){thisTF %in% unlist(y)}), "wormbase_gseq"]
    
    thisonecutoff <- min(length(x), cutoff)
    
    cut_x <- x[1:thisonecutoff]
    
    data.frame("source" = rep(thisTF_gseq, times = length(cut_x)),
               "target" = cut_x,
               "weight" = rep(1, times = length(cut_x)))
    
  }))
  
  # let's exclude TFs with fewer than 15 targets
  if(!is.null(min_targets)){
    tempGRN <- tempGRN[!tempGRN$source %in% names(table(tempGRN$source))[table(tempGRN$source) < min_targets], ]
    
    write.table(tempGRN,
                file = paste0("output/GRNs/", prefix, "_", cutoff, "_", suffix, ".txt"),
                sep = "\t",
                row.names = FALSE,
                col.names = TRUE)
    
  } else {
    
    write.table(tempGRN,
                file = paste0("output/GRNs/", prefix, "_", cutoff, "_", suffix, "_unfiltered.txt"),
                sep = "\t",
                row.names = FALSE,
                col.names = TRUE)
    
  }
  
}

cutoffs_vec <- c(100, 500, 1000, 1500, 2000, 2500, 3000, 5000, 10000)

for(i in 1:length(cutoffs_vec)){
  
  message(paste0("Now doing cutoff ", cutoffs_vec[i]))
  make.ChIP.GRN(cutoff = cutoffs_vec[i],
                targets = neworder_targets,
                prefix = "ChIPorth",
                suffix = "HOTincl",
                min_targets = 15)
  
}

for(i in 1:length(cutoffs_vec)){

  message(paste0("Now doing cutoff ", cutoffs_vec[i]))
  make.ChIP.GRN(cutoff = cutoffs_vec[i],
                targets = neworder_targets_Hotexcl,
                prefix = "ChIPorth",
                suffix = "HOTexcl",
                min_targets = 15)
  
}

make.ChIP.GRN(cutoff = 1000,
              targets = neworder_targets_Hotexcl,
              prefix = "ChIPorth",
              suffix = "HOTexcl",
              min_targets = NULL)

bench_out_filelist <- list.files("output/benchmark_out/")

ChIPorth_HOTincl_benchraptor_files <- bench_out_filelist[str_detect(bench_out_filelist, "benchRAPToR") & str_detect(bench_out_filelist, "HOTincl") & str_detect(bench_out_filelist, "ChIPorth")]

ChIPorth_HOTincl_benchraptor_files <- ChIPorth_HOTincl_benchraptor_files[!str_detect(ChIPorth_HOTincl_benchraptor_files, "dense")]
ChIPorth_HOTincl_benchraptor_files <- ChIPorth_HOTincl_benchraptor_files[!str_detect(ChIPorth_HOTincl_benchraptor_files, "shufflestat")]

ChIPorth_HOTincl_benchraptor_list <- lapply(ChIPorth_HOTincl_benchraptor_files, function(x){
  
  read.table(paste0("output/benchmark_out/", x),
             sep = "\t",
             header = TRUE)
  
})

ChIPorth_HOTincl_benchraptor_all <- do.call(rbind, ChIPorth_HOTincl_benchraptor_list)

ChIPorth_HOTincl_benchraptor_allmethods_plot <- ChIPorth_HOTincl_benchraptor_all[ChIPorth_HOTincl_benchraptor_all$metric %in% c("auroc", "auprc"), c("net", "method", "metric", "score")]
ChIPorth_HOTincl_benchraptor_allmethods_plot <- tidyr::pivot_wider(ChIPorth_HOTincl_benchraptor_allmethods_plot, names_from = c("metric"), values_from = "score")

ChIPorth_HOTincl_benchraptor_allmethods_plot[, "label"] <- as.character(ChIPorth_HOTincl_benchraptor_allmethods_plot$net)
# ChIPorth_HOTincl_benchraptor_allmethods_plot[str_detect(unlist(ChIPorth_HOTincl_benchraptor_allmethods_plot$net), "shuffle"), "label"] <- ""

ChIPorth_HOTincl_benchraptor_allmethods_plot <- ChIPorth_HOTincl_benchraptor_allmethods_plot[ChIPorth_HOTincl_benchraptor_allmethods_plot$method == "mlm_estimate", ]

allChIP_HOTincl_benchraptor_allmethods_plot <- readRDS("plotdata/allChIP_HOTincl_benchraptor_allmethods_plot.rds")
allChIP_HOTincl_benchraptor_allmethods_plot <- allChIP_HOTincl_benchraptor_allmethods_plot[allChIP_HOTincl_benchraptor_allmethods_plot$method == "mlm_estimate", ]

allChIP_HOTincl_benchraptor_allmethods_plot <- allChIP_HOTincl_benchraptor_allmethods_plot[, colnames(allChIP_HOTincl_benchraptor_allmethods_plot) %in% colnames(ChIPorth_HOTincl_benchraptor_allmethods_plot)]

allChIP_HOTincl_benchraptor_allmethods_plot[, "order"] <- "ChIPscore"
ChIPorth_HOTincl_benchraptor_allmethods_plot[, "order"] <- "orth"

both_HOTincl_benchraptor_allmethods_plot <- rbind(allChIP_HOTincl_benchraptor_allmethods_plot,
                                                  ChIPorth_HOTincl_benchraptor_allmethods_plot)

ggplot(data = both_HOTincl_benchraptor_allmethods_plot, aes(x = as.numeric(net), y = auprc, colour = order)) + 
  geom_point() + 
  geom_line() + 
  scale_x_log10(breaks = c(100, 500, 1000, 1500, 2000, 2500, 3000, 5000)) +
  ylab("Precision (AUPRC)") +
  xlab("Cutoff (log"[10]~"scale)") +
  coord_cartesian(ylim = c(0.5, 0.75)) +
  theme(strip.text = element_text(colour = 'white'),
        strip.background = element_rect(fill = "black"),
        legend.position = "none",
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_line(colour = "grey", linewidth = 0.2),
        # panel.grid.minor = element_line(colour = "grey", linewidth = 0.1),
        axis.text = element_text(colour = "black"),
        axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.5,
                                   size = 6) )

ggplot(data = both_HOTincl_benchraptor_allmethods_plot, aes(x = as.numeric(net), y = auroc, colour = order)) + 
  geom_point() + 
  geom_line() +
  scale_x_log10(breaks = c(100, 500, 1000, 1500, 2000, 2500, 3000, 5000)) +
  ylab("Precision (AUPRC)") +
  xlab("Cutoff (log"[10]~"scale)") +
  coord_cartesian(ylim = c(0.5, 0.75)) +
  theme(strip.text = element_text(colour = 'white'),
        strip.background = element_rect(fill = "black"),
        legend.position = "none",
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_line(colour = "grey", linewidth = 0.2),
        # panel.grid.minor = element_line(colour = "grey", linewidth = 0.1),
        axis.text = element_text(colour = "black"),
        axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.5,
                                   size = 6) )

ChIPorth_HOTexcl_benchraptor_files <- bench_out_filelist[str_detect(bench_out_filelist, "benchRAPToR") & str_detect(bench_out_filelist, "HOTexcl") & str_detect(bench_out_filelist, "ChIPorth")]

ChIPorth_HOTexcl_benchraptor_files <- ChIPorth_HOTexcl_benchraptor_files[!str_detect(ChIPorth_HOTexcl_benchraptor_files, "dense")]
ChIPorth_HOTexcl_benchraptor_files <- ChIPorth_HOTexcl_benchraptor_files[!str_detect(ChIPorth_HOTexcl_benchraptor_files, "shufflestat")]

ChIPorth_HOTexcl_benchraptor_list <- lapply(ChIPorth_HOTexcl_benchraptor_files, function(x){
  
  read.table(paste0("output/benchmark_out/", x),
             sep = "\t",
             header = TRUE)
  
})

ChIPorth_HOTexcl_benchraptor_all <- do.call(rbind, ChIPorth_HOTexcl_benchraptor_list)

ChIPorth_HOTexcl_benchraptor_allmethods_plot <- ChIPorth_HOTexcl_benchraptor_all[ChIPorth_HOTexcl_benchraptor_all$metric %in% c("auroc", "auprc"), c("net", "method", "metric", "score")]
ChIPorth_HOTexcl_benchraptor_allmethods_plot <- tidyr::pivot_wider(ChIPorth_HOTexcl_benchraptor_allmethods_plot, names_from = c("metric"), values_from = "score")

ChIPorth_HOTexcl_benchraptor_allmethods_plot[, "label"] <- as.character(ChIPorth_HOTexcl_benchraptor_allmethods_plot$net)
# ChIPorth_HOTexcl_benchraptor_allmethods_plot[str_detect(unlist(ChIPorth_HOTexcl_benchraptor_allmethods_plot$net), "shuffle"), "label"] <- ""

ChIPorth_HOTexcl_benchraptor_allmethods_plot <- ChIPorth_HOTexcl_benchraptor_allmethods_plot[ChIPorth_HOTexcl_benchraptor_allmethods_plot$method == "mlm_estimate", ]

allChIP_HOTexcl_benchraptor_allmethods_plot <- readRDS("plotdata/allChIP_HOTexcl_benchraptor_allmethods_plot.rds")
allChIP_HOTexcl_benchraptor_allmethods_plot <- allChIP_HOTexcl_benchraptor_allmethods_plot[allChIP_HOTexcl_benchraptor_allmethods_plot$method == "mlm_estimate", ]

allChIP_HOTexcl_benchraptor_allmethods_plot <- allChIP_HOTexcl_benchraptor_allmethods_plot[, colnames(allChIP_HOTexcl_benchraptor_allmethods_plot) %in% colnames(ChIPorth_HOTexcl_benchraptor_allmethods_plot)]

allChIP_HOTexcl_benchraptor_allmethods_plot[, "order"] <- "ChIPscore"
ChIPorth_HOTexcl_benchraptor_allmethods_plot[, "order"] <- "orth"

both_HOTexcl_benchraptor_allmethods_plot <- rbind(allChIP_HOTexcl_benchraptor_allmethods_plot,
                                                  ChIPorth_HOTexcl_benchraptor_allmethods_plot)

ggplot(data = both_HOTexcl_benchraptor_allmethods_plot, aes(x = as.numeric(net), y = auprc, colour = order)) + 
  geom_point() + 
  geom_line() + 
  scale_x_log10(breaks = c(100, 500, 1000, 1500, 2000, 2500, 3000, 5000)) +
  ylab("Precision (AUPRC)") +
  xlab("Cutoff (log"[10]~"scale)") +
  coord_cartesian(ylim = c(0.5, 0.75)) +
  theme(strip.text = element_text(colour = 'white'),
        strip.background = element_rect(fill = "black"),
        legend.position = "none",
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_line(colour = "grey", linewidth = 0.2),
        # panel.grid.minor = element_line(colour = "grey", linewidth = 0.1),
        axis.text = element_text(colour = "black"),
        axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.5,
                                   size = 6) )

ggplot(data = both_HOTexcl_benchraptor_allmethods_plot, aes(x = as.numeric(net), y = auroc, colour = order)) + 
  geom_point() +
  geom_line() +
  scale_x_log10(breaks = c(100, 500, 1000, 1500, 2000, 2500, 3000, 5000)) +
  ylab("Precision (AUPRC)") +
  xlab("Cutoff (log"[10]~"scale)") +
  coord_cartesian(ylim = c(0.5, 0.75)) +
  theme(strip.text = element_text(colour = 'white'),
        strip.background = element_rect(fill = "black"),
        legend.position = "none",
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_line(colour = "grey", linewidth = 0.2),
        # panel.grid.minor = element_line(colour = "grey", linewidth = 0.1),
        axis.text = element_text(colour = "black"),
        axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.5,
                                   size = 6) )


#### TRY WITH CISBP PROBS FOR THOSE ONES ####

TF_orthology_probs <- readRDS("output/TF_orthology_probs10k.rds")

STREMECisBP_orthology_probs <- STREMEtopmotif_orthology_probs

shared_TFs <- intersect(names(TF_orthology_probs),
                        names(STREMECisBP_orthology_probs))

STREMECisBP_orthology_probs[shared_TFs] <- TF_orthology_probs[shared_TFs]

STREMECisBPorder_targets_Hotexcl <- lapply(TFs, function(x){
  print(x)
  temp <- allTFS_manualtoptargets_HOTexcl[[allTFs_labelused[x]]]
  
  if(!is.null(STREMECisBP_orthology_probs[[x]][match(temp, names(STREMECisBP_orthology_probs[[x]]))])){
    
    tempout <- data.frame("target" = temp, "probs" = 1) 
    tempout[, "probs"] <- STREMECisBP_orthology_probs[[x]][match(temp, names(STREMECisBP_orthology_probs[[x]]))]
    
    tempout[is.na(tempout$probs), "probs"] <- 1
    
    tempout <- tempout[order(tempout$probs), ]
    
    # order scores of 1 according to ChIP score. although it seems this order is retained, anyway
    
    tempout[tempout$probs == 1, "target"] <- temp[temp %in% tempout[tempout$probs == 1, "target"]]
    
    return(tempout$target)
    
  } else {
    
    return(temp)
    
  }
  
})

names(STREMECisBPorder_targets_Hotexcl) <- TFs

for(i in 1:length(cutoffs_vec)){
  
  message(paste0("Now doing cutoff ", cutoffs_vec[i]))
  make.ChIP.GRN(cutoff = cutoffs_vec[i],
                targets = STREMECisBPorder_targets_Hotexcl,
                prefix = "STREMECisBP",
                suffix = "HOTexcl",
                min_targets = 15)
  
}

### NOW RUN PLOT IT ####

bench_out_filelist <- list.files("output/benchmark_out/")

STREMECisBP_HOTexcl_benchraptor_files <- bench_out_filelist[str_detect(bench_out_filelist, "benchRAPToR") &
                                                              # str_detect(bench_out_filelist, "HOTexcl") &
                                                              str_detect(bench_out_filelist, "STRCisBP")]

STREMECisBP_HOTexcl_benchraptor_files <- STREMECisBP_HOTexcl_benchraptor_files[!str_detect(STREMECisBP_HOTexcl_benchraptor_files, "dense")]
STREMECisBP_HOTexcl_benchraptor_files <- STREMECisBP_HOTexcl_benchraptor_files[!str_detect(STREMECisBP_HOTexcl_benchraptor_files, "shufflestat")]

STREMECisBP_HOTexcl_benchraptor_list <- lapply(STREMECisBP_HOTexcl_benchraptor_files, function(x){
  
  read.table(paste0("output/benchmark_out/", x),
             sep = "\t",
             header = TRUE)
  
})

STREMECisBP_HOTexcl_benchraptor_all <- do.call(rbind, STREMECisBP_HOTexcl_benchraptor_list)

STREMECisBP_HOTexcl_benchraptor_allmethods_plot <- STREMECisBP_HOTexcl_benchraptor_all[STREMECisBP_HOTexcl_benchraptor_all$metric %in% c("auroc", "auprc"), c("net", "method", "metric", "score")]
STREMECisBP_HOTexcl_benchraptor_allmethods_plot <- tidyr::pivot_wider(STREMECisBP_HOTexcl_benchraptor_allmethods_plot, names_from = c("metric"), values_from = "score")

STREMECisBP_HOTexcl_benchraptor_allmethods_plot[, "label"] <- as.character(STREMECisBP_HOTexcl_benchraptor_allmethods_plot$net)
# STREMECisBP_HOTexcl_benchraptor_allmethods_plot[str_detect(unlist(STREMECisBP_HOTexcl_benchraptor_allmethods_plot$net), "shuffle"), "label"] <- ""

STREMECisBP_HOTexcl_benchraptor_allmethods_plot <- STREMECisBP_HOTexcl_benchraptor_allmethods_plot[STREMECisBP_HOTexcl_benchraptor_allmethods_plot$method == "mlm_estimate", ]

allChIP_HOTexcl_benchraptor_allmethods_plot <- readRDS("plotdata/allChIP_HOTexcl_benchraptor_allmethods_plot.rds")
allChIP_HOTexcl_benchraptor_allmethods_plot <- allChIP_HOTexcl_benchraptor_allmethods_plot[allChIP_HOTexcl_benchraptor_allmethods_plot$method == "mlm_estimate", ]

allChIP_HOTexcl_benchraptor_allmethods_plot <- allChIP_HOTexcl_benchraptor_allmethods_plot[, colnames(allChIP_HOTexcl_benchraptor_allmethods_plot) %in% colnames(STREMECisBP_HOTexcl_benchraptor_allmethods_plot)]

allChIP_HOTexcl_benchraptor_allmethods_plot[, "order"] <- "ChIPscore"
STREMECisBP_HOTexcl_benchraptor_allmethods_plot[, "order"] <- "STRCisBP"

three_HOTexcl_benchraptor_allmethods_plot <- rbind(allChIP_HOTexcl_benchraptor_allmethods_plot,
                                                  STREMECisBP_HOTexcl_benchraptor_allmethods_plot,
                                                  ChIPorth_HOTexcl_benchraptor_allmethods_plot)

ggplot(data = three_HOTexcl_benchraptor_allmethods_plot, aes(x = as.numeric(net), y = auprc, colour = order)) + 
  geom_point() + 
  geom_line() + 
  scale_x_log10(breaks = c(100, 500, 1000, 1500, 2000, 2500, 3000, 5000)) +
  ylab("Precision (AUPRC)") +
  xlab("Cutoff (log"[10]~"scale)") +
  coord_cartesian(ylim = c(0.5, 0.75)) +
  theme(strip.text = element_text(colour = 'white'),
        strip.background = element_rect(fill = "black"),
        # legend.position = "none",
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_line(colour = "grey", linewidth = 0.2),
        # panel.grid.minor = element_line(colour = "grey", linewidth = 0.1),
        axis.text = element_text(colour = "black"),
        axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.5,
                                   size = 6) )

ggplot(data = three_HOTexcl_benchraptor_allmethods_plot, aes(x = as.numeric(net), y = auroc, colour = order)) + 
  geom_point() + 
  geom_line() +
  scale_x_log10(breaks = c(100, 500, 1000, 1500, 2000, 2500, 3000, 5000)) +
  ylab("Precision (AUPRC)") +
  xlab("Cutoff (log"[10]~"scale)") +
  coord_cartesian(ylim = c(0.5, 0.75)) +
  theme(strip.text = element_text(colour = 'white'),
        strip.background = element_rect(fill = "black"),
        # legend.position = "none",
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_line(colour = "grey", linewidth = 0.2),
        # panel.grid.minor = element_line(colour = "grey", linewidth = 0.1),
        axis.text = element_text(colour = "black"),
        axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.5,
                                   size = 6) )
