#start with empty workspace

rm(list = ls(all = TRUE))


# clear all loaded packages
# invisible(lapply(paste0("package:", names(sessionInfo()$otherPkgs)),
#                  detach,
#                  character.only = TRUE, unload = TRUE))

# turn off scientific notation for plots

options(scipen=10000)

#### set working directory ####

# here create new folder and set working directory within it

setwd("~/Cel_GRN_manuscript")

#### DEFINE FUNCTIONS ####

ecdf_fun <- function(x,perc) ecdf(x)(perc)

make.TF.orth.GRN <- function(cutoff,
                             fdrcut,
                             TF_orthology_prob_list = TF_orthology_probs,
                             elegans_targets = elegans_target_table,
                             suffix = NULL){
  # cutoff = 2000
  # fdrcut = 0
  tempGRN <- do.call(rbind, lapply(names(TF_orthology_probs), function(thisTF){
    # thisTF <- names(TF_orthology_probs)[5]
    x <- TF_orthology_probs[[thisTF]]
    
    fdrcutx <- names(x[!x > fdrcut])
    
    elegans_x <- elegans_target_table[, thisTF]
    names(elegans_x) <- row.names(elegans_target_table)
    
    elegans_x_order <- elegans_x[order(elegans_x, decreasing = TRUE)]
    
    thiscutoff <- min(cutoff, sum(elegans_x_order != 0))
    
    fdrcut_names <- names(elegans_x_order[1:thiscutoff])[names(elegans_x_order[1:thiscutoff]) %in% fdrcutx]
    
    data.frame("source" = rep(thisTF, times = length(fdrcut_names)),
               "target" = fdrcut_names,
               "weight" = rep(1, times = length(fdrcut_names)))
    
  }))
  
  # let's exclude fewer than 15
  tempGRN <- tempGRN[!tempGRN$source %in% names(table(tempGRN$source))[table(tempGRN$source) < 15], ]
  
  write.table(tempGRN,
              file = paste0("output/GRNs/TF_orthprobs_cut", cutoff, "_fdr", fdrcut, suffix,".txt"),
              sep = "\t",
              row.names = FALSE,
              col.names = TRUE)
  
}

#### LOAD PACKAGES & FUNCTIONS ####

## First specify the packages of interest

packages <- c("stringr",
              "biomaRt",
              "readxl", # for reading xls files (NOTE not xlsx)
              "ggplot2",
              "ggrepel") 

## Now load or install&load all
package.check <- lapply(
  
  packages,
  
  FUN = function(x) {
    
    if (!require(x, character.only = TRUE)) {
      
      install.packages(x, dependencies = TRUE)
      
      library(x, character.only = TRUE)
      
    }
    
  }
  
)

# Here define packages which need to be loaded through biocmanager

biocmanager_packages <- c() 

bioc_package.check <- lapply(
  
  biocmanager_packages,
  
  FUN = function(x) {
    
    if (!require(x, character.only = TRUE)) {
      
      if (!requireNamespace("BiocManager", quietly = TRUE)){
        
        install.packages("BiocManager")
        
      }
      
      BiocManager::install(x, dependencies = TRUE)
      
      library(x, character.only = TRUE)
      
    }
  }
)

# Here define packages to be loaded through Github

github_packages <- c("r-lib/conflicted",
                     "etam4260/kneedle") # kneedle to automatially find elbow, if that works

if (!requireNamespace("devtools", quietly = TRUE)){
  install.packages("devtools")}

gh_package.check <- lapply(
  github_packages,
  FUN = function(x) {
    if (!require(str_remove(x, ".*\\/"), character.only = TRUE)) {
      
      if (!requireNamespace("devtools", quietly = TRUE)){
        install.packages("devtools")}
      
      devtools::install_github(x, build_vignettes = TRUE)
      
      library(str_remove(x, ".*\\/"), character.only = TRUE)
      
    }
  }
)

#### INPUT DATA ####

bench_out_filelist <- list.files("output/benchmark_out/")
GRN_filelist <- list.files("output/GRNs") 

TForthprobs_files <- bench_out_filelist[str_detect(bench_out_filelist, "TF_orth")]
TForthprobs_GRNs <- GRN_filelist[str_detect(GRN_filelist, "TF_orth")]

cutoffs_vec <- c(500, 1000, 1500, 2000, 2500, 3000, 5000, 10000)
fdr_vec <- c(0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8)

method_labels <- c("Consensus", "MLM", "ULM", "WSum")
names(method_labels) <- c("consensus_estimate", "mlm_estimate", "ulm_estimate", "wsum_norm")

TFnumbers <- sapply(cutoffs_vec, function(thiscutoff){

  sapply(fdr_vec, function(thisfdr){

    tempGRN <- read.table(paste0("output/GRNs/TF_orthprobs_cut", thiscutoff, "_fdr", thisfdr, ".txt"),
                          header = TRUE)
    
    length(unique(tempGRN$source))
    
  })
  
})

colnames(TFnumbers) <- cutoffs_vec
row.names(TFnumbers) <- fdr_vec

TF_meantarget_numbers <- sapply(cutoffs_vec, function(thiscutoff){
  
  sapply(fdr_vec, function(thisfdr){
    
    tempGRN <- read.table(paste0("output/GRNs/TF_orthprobs_cut", thiscutoff, "_fdr", thisfdr, ".txt"),
                          header = TRUE)
    
    mean(table(tempGRN$source))
    
  })
  
})

colnames(TF_meantarget_numbers) <- cutoffs_vec
row.names(TF_meantarget_numbers) <- fdr_vec

TF_mediantarget_numbers <- sapply(cutoffs_vec, function(thiscutoff){
  
  sapply(fdr_vec, function(thisfdr){
    
    tempGRN <- read.table(paste0("output/GRNs/TF_orthprobs_cut", thiscutoff, "_fdr", thisfdr, ".txt"),
                          header = TRUE)
    
    median(table(tempGRN$source))
    
  })
  
})

colnames(TF_mediantarget_numbers) <- cutoffs_vec
row.names(TF_mediantarget_numbers) <- fdr_vec

for(i in 1:length(fdr_vec)){

TForthprobs_tempfdr_benchraptor_list <- lapply(TForthprobs_files[str_detect(TForthprobs_files, paste0("fdr", fdr_vec[i], "_bench"))], function(x){
         
read.table(paste0("output/benchmark_out/", x),
                    sep = "\t",
                    header = TRUE)
         
})
       
TForthprobs_tempfdr_benchraptor_all <- do.call(rbind, TForthprobs_tempfdr_benchraptor_list)
       
TForthprobs_tempfdr_benchraptor_allmethods_plot <- TForthprobs_tempfdr_benchraptor_all[TForthprobs_tempfdr_benchraptor_all$metric %in% c("auroc", "auprc"), c("net", "method", "metric", "score")]
TForthprobs_tempfdr_benchraptor_allmethods_plot <- tidyr::pivot_wider(TForthprobs_tempfdr_benchraptor_allmethods_plot, names_from = c("metric"), values_from = "score")
       
TForthprobs_tempfdr_benchraptor_allmethods_plot[, "label"] <- as.character(TForthprobs_tempfdr_benchraptor_allmethods_plot$net)

pdf(paste0("graphics/TForth_methodsplot_fdr", fdr_vec[i], ".pdf"),
    width = 5,
    height = 4)

print(
  ggplot(data = TForthprobs_tempfdr_benchraptor_allmethods_plot[!TForthprobs_tempfdr_benchraptor_allmethods_plot$method %in% c("wsum_estimate", "wsum_corr"), ], aes(x = auroc, y = auprc, colour = net)) + 
  geom_point(size = 2) +
  geom_label_repel(aes(label = net),
                   label.size = 0.1,
                   max.overlaps = 15,
                   box.padding = 0.3,
                   label.padding = 0.1,
                   point.padding = 0,
                   size = 2) +
  geom_vline(xintercept = 0.5,
             linetype = "dashed",
             col = "black") +
  geom_hline(yintercept = 0.5,
             linetype = "dashed",
             col = "black") + 
  facet_wrap(~ method, labeller = labeller(method = method_labels)) + 
  xlab("Sensitivity (AUROC)") + 
  ylab("Precision (AUPRC)") +
  coord_cartesian(xlim = c(0.45, 0.75),
                  ylim = c(0.45, 0.8)) +
  theme(strip.text = element_text(colour = 'white'),
        strip.background = element_rect(fill = "black"),
        legend.position = "none",
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_line(colour = "grey", linewidth = 0.2),
        panel.grid.minor = element_line(colour = "grey", linewidth = 0.1),
        axis.text = element_text(colour = "black"))
)
  
dev.off()

}

for(i in 1:length(cutoffs_vec)){

  TForthprobs_tempcut_benchraptor_list <- lapply(TForthprobs_files[str_detect(TForthprobs_files, paste0("cut", cutoffs_vec[i], "_"))], function(x){
    
    read.table(paste0("output/benchmark_out/", x),
               sep = "\t",
               header = TRUE)
    
  })
  
  TForthprobs_tempcut_benchraptor_all <- do.call(rbind, TForthprobs_tempcut_benchraptor_list)
  
  TForthprobs_tempcut_benchraptor_allmethods_plot <- TForthprobs_tempcut_benchraptor_all[TForthprobs_tempcut_benchraptor_all$metric %in% c("auroc", "auprc"), c("net", "method", "metric", "score")]
  TForthprobs_tempcut_benchraptor_allmethods_plot <- tidyr::pivot_wider(TForthprobs_tempcut_benchraptor_allmethods_plot, names_from = c("metric"), values_from = "score")
  
  TForthprobs_tempcut_benchraptor_allmethods_plot[, "label"] <- as.character(TForthprobs_tempcut_benchraptor_allmethods_plot$net)
  
  pdf(paste0("graphics/TForth_methodsplot_fdr", cutoffs_vec[i], ".pdf"),
      width = 5,
      height = 4)
  
  print(
    ggplot(data = TForthprobs_tempcut_benchraptor_allmethods_plot[!TForthprobs_tempcut_benchraptor_allmethods_plot$method %in% c("wsum_estimate", "wsum_corr"), ], aes(x = auroc, y = auprc, colour = net)) + 
      geom_point(size = 2) +
      geom_label_repel(aes(label = net),
                       label.size = 0.1,
                       max.overlaps = 15,
                       box.padding = 0.3,
                       label.padding = 0.1,
                       point.padding = 0,
                       size = 2) +
      geom_vline(xintercept = 0.5,
                 linetype = "dashed",
                 col = "black") +
      geom_hline(yintercept = 0.5,
                 linetype = "dashed",
                 col = "black") + 
      facet_wrap(~ method, labeller = labeller(method = method_labels)) + 
      xlab("Sensitivity (AUROC)") + 
      ylab("Precision (AUPRC)") +
      coord_cartesian(xlim = c(0.45, 0.75),
                      ylim = c(0.45, 0.8)) +
      theme(strip.text = element_text(colour = 'white'),
            strip.background = element_rect(fill = "black"),
            legend.position = "none",
            panel.background = element_rect(fill = "white", colour = "black"),
            panel.grid.major = element_line(colour = "grey", linewidth = 0.2),
            panel.grid.minor = element_line(colour = "grey", linewidth = 0.1),
            axis.text = element_text(colour = "black"))
  )
  
  dev.off()
  
}


TF_auroc_mlm <- sapply(cutoffs_vec, function(thiscutoff){
  
  sapply(fdr_vec, function(thisfdr){

    tempbench <- read.table(paste0("output/benchmark_out/TF_orthprobs_cut", thiscutoff, "_fdr", thisfdr, "_benchRAPToR.tsv"),
                          sep = "\t",
                          header = TRUE)
    
    tempbench[tempbench$method == "mlm_estimate" & tempbench$metric == "auroc", "score"]
    
  })
  
})

colnames(TF_auroc_mlm) <- cutoffs_vec
row.names(TF_auroc_mlm) <- fdr_vec

TF_auroc_ulm <- sapply(cutoffs_vec, function(thiscutoff){
  
  sapply(fdr_vec, function(thisfdr){
    
    tempbench <- read.table(paste0("output/benchmark_out/TF_orthprobs_cut", thiscutoff, "_fdr", thisfdr, "_benchRAPToR.tsv"),
                            sep = "\t",
                            header = TRUE)
    
    tempbench[tempbench$method == "ulm_estimate" & tempbench$metric == "auroc", "score"]
    
  })
  
})

TF_auroc_wsum <- sapply(cutoffs_vec, function(thiscutoff){
  
  sapply(fdr_vec, function(thisfdr){
    
    tempbench <- read.table(paste0("output/benchmark_out/TF_orthprobs_cut", thiscutoff, "_fdr", thisfdr, "_benchRAPToR.tsv"),
                            sep = "\t",
                            header = TRUE)
    
    tempbench[tempbench$method == "wsum_estimate" & tempbench$metric == "auroc", "score"]
    
  })
  
})

colnames(TF_auroc_wsum) <- cutoffs_vec
row.names(TF_auroc_wsum) <- fdr_vec

TF_auroc_consensus <- sapply(cutoffs_vec, function(thiscutoff){
  
  sapply(fdr_vec, function(thisfdr){
    
    tempbench <- read.table(paste0("output/benchmark_out/TF_orthprobs_cut", thiscutoff, "_fdr", thisfdr, "_benchRAPToR.tsv"),
                            sep = "\t",
                            header = TRUE)
    
    tempbench[tempbench$method == "consensus_estimate" & tempbench$metric == "auroc", "score"]
    
  })
  
})

colnames(TF_auroc_consensus) <- cutoffs_vec
row.names(TF_auroc_consensus) <- fdr_vec

TF_auprc_mlm <- sapply(cutoffs_vec, function(thiscutoff){
  
  sapply(fdr_vec, function(thisfdr){
    
    tempbench <- read.table(paste0("output/benchmark_out/TF_orthprobs_cut", thiscutoff, "_fdr", thisfdr, "_benchRAPToR.tsv"),
                            sep = "\t",
                            header = TRUE)
    
    tempbench[tempbench$method == "mlm_estimate" & tempbench$metric == "auprc", "score"]
    
  })
  
})

colnames(TF_auprc_mlm) <- cutoffs_vec
row.names(TF_auprc_mlm) <- fdr_vec

TF_auprc_ulm <- sapply(cutoffs_vec, function(thiscutoff){
  
  sapply(fdr_vec, function(thisfdr){
    
    tempbench <- read.table(paste0("output/benchmark_out/TF_orthprobs_cut", thiscutoff, "_fdr", thisfdr, "_benchRAPToR.tsv"),
                            sep = "\t",
                            header = TRUE)
    
    tempbench[tempbench$method == "ulm_estimate" & tempbench$metric == "auprc", "score"]
    
  })
  
})

colnames(TF_auprc_ulm) <- cutoffs_vec
row.names(TF_auprc_ulm) <- fdr_vec

TF_auprc_consensus <- sapply(cutoffs_vec, function(thiscutoff){
  
  sapply(fdr_vec, function(thisfdr){
    
    tempbench <- read.table(paste0("output/benchmark_out/TF_orthprobs_cut", thiscutoff, "_fdr", thisfdr, "_benchRAPToR.tsv"),
                            sep = "\t",
                            header = TRUE)
    
    tempbench[tempbench$method == "consensus_estimate" & tempbench$metric == "auprc", "score"]
    
  })
  
})

colnames(TF_auprc_consensus) <- cutoffs_vec
row.names(TF_auprc_consensus) <- fdr_vec

library(gplots)
mainpal <- (colorRampPalette(c("red", "orange", "yellow", "white"))(100))
?redblue
colors = seq(0.6, 0.7, length.out=101)
redblue(100)

pdf("graphics/TForth_cutoffheatmap_auprc_mlm.pdf",
    width = 5,
    height = 5)

gplots::heatmap.2(x = t(TF_auprc_mlm),
                  dendrogram = "none",
                  trace = "none",
                  Rowv = FALSE,
                  Colv = FALSE,
                  density.info = "none",
                  col = mainpal,
                  breaks = colors,
                  key = FALSE,
                  xlab = "TF-target conservation FDR",
                  ylab = "Max. targets cutoff")

dev.off()

pdf("graphics/TForth_cutoffheatmap_auroc_mlm.pdf",
    width = 5,
    height = 5)

gplots::heatmap.2(x = t(TF_auroc_mlm),
                  dendrogram = "none",
                  trace = "none",
                  Rowv = FALSE,
                  Colv = FALSE,
                  density.info = "none",
                  col = mainpal,
                  breaks = colors,
                  key = FALSE)

dev.off()

pdf("graphics/TForth_cutoffheatmap_auprc_mlm_KEY.pdf",
    width = 5,
    height = 5)

gplots::heatmap.2(x = t(TF_auprc_mlm),
                  dendrogram = "none",
                  trace = "none",
                  Rowv = FALSE,
                  Colv = FALSE,
                  density.info = "none",
                  col = mainpal,
                  breaks = colors,
                  key = TRUE,
                  xlab = "TF-target conservation FDR",
                  ylab = "Max. targets cutoff")

dev.off()

# other methods


TF_auroc_ulm <- sapply(cutoffs_vec, function(thiscutoff){
  
  sapply(fdr_vec, function(thisfdr){
    
    tempbench <- read.table(paste0("output/benchmark_out/TF_orthprobs_cut", thiscutoff, "_fdr", thisfdr, "_benchRAPToR.tsv"),
                            sep = "\t",
                            header = TRUE)
    
    tempbench[tempbench$method == "ulm_estimate" & tempbench$metric == "auroc", "score"]
    
  })
  
})

colnames(TF_auroc_ulm) <- cutoffs_vec
row.names(TF_auroc_ulm) <- fdr_vec

TF_auprc_ulm <- sapply(cutoffs_vec, function(thiscutoff){
  
  sapply(fdr_vec, function(thisfdr){
    thiscutoff <- cutoffs_vec[1]
    thisfdr <- fdr_vec[1]
    tempbench <- read.table(paste0("output/benchmark_out/TF_orthprobs_cut", thiscutoff, "_fdr", thisfdr, "_benchRAPToR.tsv"),
                            sep = "\t",
                            header = TRUE)
    
    tempbench[tempbench$method == "ulm_estimate" & tempbench$metric == "auprc", "score"]
    
  })
  
})

colnames(TF_auprc_ulm) <- cutoffs_vec
row.names(TF_auprc_ulm) <- fdr_vec




#### MAIN PLOT ####

TForth_benchraptor_shufflestat_files <- paste0(c("TForth", "CelEsT", "control", "FIMOfull"), "_benchRAPToR_shufflestats.tsv")

TForth_benchraptor <- read.table("output/benchmark_out/FIMO_TForth_df.tsv",
                                 sep= '\t',
                                 header = TRUE)

TForth_benchraptormethods_plot <- TForth_benchraptor[TForth_benchraptor$metric %in% c("auroc", "auprc"), c("net", "method", "metric", "score")]
TForth_benchraptormethods_plot <- tidyr::pivot_wider(TForth_benchraptormethods_plot, names_from = c("metric"), values_from = "score")

TForth_benchraptormethods_plot[, "label"] <- TForth_benchraptormethods_plot$net
TForth_benchraptormethods_plot[str_detect(unlist(TForth_benchraptormethods_plot$net), "shuffle"), "label"] <- ""

TForth_labelvec <- c("OrthNet", "ctrlNet", "full1500", "CelEsT")
names(TForth_labelvec) <- unique(TForth_benchraptormethods_plot$net)

TForth_benchraptormethods_plot[, "label"] <- TForth_labelvec[TForth_benchraptormethods_plot$net]

TForth_benchraptor_shufflestats_lists <- lapply(TForth_benchraptor_shufflestat_files, function(file){
  
  thisdata <- read.table(paste0("output/benchmark_out/", file),
                         sep = "\t",
                         header = TRUE)
  
  auprc_stats <- thisdata[thisdata$metric == "auprc", ] %>% group_by(method) %>% summarise(mean_auprc = mean(score),
                                                                                           sd_auprc = sd(score))
  
  auprc_stats[, "auprc_95CI_margin"] <- qt(0.975, df = 99) * auprc_stats$sd_auprc / sqrt(99)
  
  auroc_stats <- thisdata[thisdata$metric == "auroc", ] %>% group_by(method) %>% summarise(mean_auroc = mean(score),
                                                                                           sd_auroc = sd(score))
  
  auroc_stats[, "auroc_95CI_margin"] <- qt(0.975, df = 99) * auroc_stats$sd_auroc / sqrt(99)
  
  tempout <- merge(auprc_stats,
                   auroc_stats)
  
  tempout["net"] <- str_remove(file, "_benchRAPToR_shufflestats.tsv")
  
  tempout
  
})

TForth_benchraptor_shufflestats_df <- do.call(rbind, TForth_benchraptor_shufflestats_lists)

TForth_benchraptor_shufflestats_df[, "label"] <- TForth_labelvec[TForth_benchraptor_shufflestats_df$net]

saveRDS(TForth_benchraptormethods_plot,
        "plotdata/TForth_benchraptormethods_plot.rds")

saveRDS(TForth_benchraptor_shufflestats_df,
        "plotdata/TForth_benchraptor_shufflestats_df.rds")

pdf("graphics/TForth_benchRAPToR_methodsplot.pdf",
    height = 4,
    width = 5)

ggplot(data = TForth_benchraptormethods_plot[!TForth_benchraptormethods_plot$method %in% c("wsum_estimate", "wsum_corr") &
                                               TForth_benchraptormethods_plot$net != "TForth_unequal", ], aes(x = auroc, y = auprc, colour = net)) + 
  geom_point(size = 2) +
  geom_label_repel(aes(label = label),
                   label.size = 0.1,
                   max.overlaps = 15,
                   box.padding = 0.3,
                   label.padding = 0.1,
                   point.padding = 0,
                   size = 2) +
  geom_vline(xintercept = 0.5,
             linetype = "dashed",
             col = "black") +
  geom_hline(yintercept = 0.5,
             linetype = "dashed",
             col = "black") + 
  facet_wrap(~ method, labeller = labeller(method = method_labels)) + 
  xlab("Sensitivity (AUROC)") + 
  ylab("Precision (AUPRC)") +
  coord_cartesian(xlim = c(0.45, 0.75),
                  ylim = c(0.45, 0.8)) +
  theme(strip.text = element_text(colour = 'white'),
        strip.background = element_rect(fill = "black"),
        legend.position = "none",
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_line(colour = "grey", linewidth = 0.2),
        panel.grid.minor = element_line(colour = "grey", linewidth = 0.1),
        axis.text = element_text(colour = "black")) +
  geom_point(data = TForth_benchraptor_shufflestats_df[!TForth_benchraptor_shufflestats_df$method %in% c("wsum_estimate", "wsum_corr"), ], aes(x = mean_auroc, y = mean_auprc), alpha = 0.5) + 
  geom_errorbar(data = TForth_benchraptor_shufflestats_df[!TForth_benchraptor_shufflestats_df$method %in% c("wsum_estimate", "wsum_corr"), ], aes(x = mean_auroc, y = mean_auprc,
                                                                                                                                                  ymin = (mean_auprc - sd_auprc),
                                                                                                                                                  ymax = (mean_auprc + sd_auprc)), alpha = 0.5) +
  geom_errorbar(data = TForth_benchraptor_shufflestats_df[!TForth_benchraptor_shufflestats_df$method %in% c("wsum_estimate", "wsum_corr"), ], aes(x = mean_auroc, y = mean_auprc,
                                                                                                                                                  xmin = (mean_auroc - sd_auroc),
                                                                                                                                                  xmax = (mean_auroc + sd_auroc)), alpha = 0.5)

dev.off()

pdf("graphics/TForth_benchraptor_mlmplot.pdf",
    height = 3.5,
    width = 3.5)

ggplot(data = TForth_benchraptormethods_plot[!str_detect(TForth_benchraptormethods_plot$net, "shuffle") & 
                                               TForth_benchraptormethods_plot$method == "mlm_estimate" &
                                               TForth_benchraptormethods_plot$net != "TForth_unequal", ], aes(x = auroc, y = auprc, colour = net)) + 
  geom_point(size = 2) +
  geom_label_repel(aes(label = label),
                   label.size = 0.1,
                   max.overlaps = 10,
                   box.padding = 0.2,
                   label.padding = 0.1,
                   point.padding = 0,
                   size = 3) +
  geom_vline(xintercept = 0.5,
             linetype = "dashed",
             col = "black") +
  geom_hline(yintercept = 0.5,
             linetype = "dashed",
             col = "black") + 
  xlab("Sensitivity (AUROC)") + 
  ylab("Precision (AUPRC)") +
  coord_cartesian(xlim = c(0.45, 0.75),
                  ylim = c(0.45, 0.8)) +
  theme(strip.text = element_text(colour = 'white'),
        strip.background = element_rect(fill = "black"),
        legend.position = "none",
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_line(colour = "grey", linewidth = 0.2),
        panel.grid.minor = element_line(colour = "grey", linewidth = 0.1),
        axis.text = element_text(colour = "black")) +
  geom_point(data = TForth_benchraptor_shufflestats_df[TForth_benchraptor_shufflestats_df$method == "mlm_estimate", ], aes(x = mean_auroc, y = mean_auprc), alpha = 0.5) + 
  geom_errorbar(data = TForth_benchraptor_shufflestats_df[TForth_benchraptor_shufflestats_df$method == "mlm_estimate", ], aes(x = mean_auroc, y = mean_auprc,
                                                                                                                              ymin = (mean_auprc - sd_auprc),
                                                                                                                              ymax = (mean_auprc + sd_auprc)), alpha = 0.5) +
  geom_errorbar(data = TForth_benchraptor_shufflestats_df[TForth_benchraptor_shufflestats_df$method == "mlm_estimate", ], aes(x = mean_auroc, y = mean_auprc,
                                                                                                                              xmin = (mean_auroc - sd_auroc),
                                                                                                                              xmax = (mean_auroc + sd_auroc)), alpha = 0.5)
dev.off()
