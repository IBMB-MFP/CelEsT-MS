#start with empty workspace

rm(list = ls(all = TRUE))

# turn off scientific notation for plots

options(scipen=10000)

#### set working directory ####

# here create new folder and set working directory within it

dir.create("~/Cel_GRN_revisions/")
setwd("~/Cel_GRN_revisions/")

# create subfolders for input, output and graphics

dir.create("input")

# into input folder, add input files 

dir.create("output")

dir.create("output/GRNs")

dir.create("graphics")

#### DEFINE FUNCTIONS ####


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

#### LOAD PACKAGES & FUNCTIONS ####

## First specify the packages of interest

packages <- c("ggplot2",
              "stringr",
              "ggrepel",
              "dplyr")

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

biocmanager_packages <- c("reshape2") 

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

github_packages <- c("r-lib/conflicted")

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

# allmodERN_TFsONLY_BM <- readRDS("~/Cel_GRN_manuscript/output/allmodERN_TFsONLY_BM.rds")
allmodERN_TFsONLY_BM <- readRDS("output/allmodERN_TFsONLY_BM.rds")


# allTFS_manualtoptargets <- readRDS("~/Cel_GRN_manuscript/output/MODernENCODE_manualtoptargets_operonexcluded.rds")
allTFS_manualtoptargets <- readRDS("output/MODernENCODE_manualtoptargets_operonexcluded.rds")

# allTFS_manualtoptargets_HOTexcl <- readRDS("~/Cel_GRN_manuscript/output/MODernENCODE_manualtoptargets_operon&HOTexcluded.rds")
allTFS_manualtoptargets_HOTexcl <- readRDS("output/MODernENCODE_manualtoptargets_operon&HOTexcluded.rds")

method_labels <- c("Consensus", "MLM", "ULM", "WSum")
names(method_labels) <- c("consensus_estimate", "mlm_estimate", "ulm_estimate", "wsum_norm")

bench_out_filelist <- list.files("output/benchmark_out/")
shufflestats_files <- bench_out_filelist[str_detect(bench_out_filelist, "_shufflestats")]

#### CHIP PLOTS ####

## and with RAPToR-corrected benchmark stats

allChIP_HOTincl_benchraptor_files <- bench_out_filelist[str_detect(bench_out_filelist, "benchRAPToR") &
                                                          str_detect(bench_out_filelist, "HOTincl") &
                                                                       str_detect(bench_out_filelist, "allChIP")
                                                          ]

allChIP_HOTincl_benchraptor_files <- allChIP_HOTincl_benchraptor_files[!str_detect(allChIP_HOTincl_benchraptor_files, "dense")]
allChIP_HOTincl_benchraptor_files <- allChIP_HOTincl_benchraptor_files[!str_detect(allChIP_HOTincl_benchraptor_files, "shufflestat")]

allChIP_HOTincl_benchraptor_list <- lapply(allChIP_HOTincl_benchraptor_files, function(x){
  
  read.table(paste0("output/benchmark_out/", x),
             sep = "\t",
             header = TRUE)
  
})

allChIP_HOTincl_benchraptor_all <- do.call(rbind, allChIP_HOTincl_benchraptor_list)

allChIP_HOTincl_benchraptor_allmethods_plot <- allChIP_HOTincl_benchraptor_all[allChIP_HOTincl_benchraptor_all$metric %in% c("auroc", "auprc"), c("net", "method", "metric", "score")]
allChIP_HOTincl_benchraptor_allmethods_plot <- tidyr::pivot_wider(allChIP_HOTincl_benchraptor_allmethods_plot, names_from = c("metric"), values_from = "score")

allChIP_HOTincl_benchraptor_allmethods_plot[, "label"] <- as.character(allChIP_HOTincl_benchraptor_allmethods_plot$net)
# allChIP_HOTincl_benchraptor_allmethods_plot[str_detect(unlist(allChIP_HOTincl_benchraptor_allmethods_plot$net), "shuffle"), "label"] <- ""

allChIP_HOTincl_benchraptor_allmethods_plot[, "cutoff"] <- str_remove(allChIP_HOTincl_benchraptor_allmethods_plot$net, "_shuffle")
allChIP_HOTincl_benchraptor_allmethods_plot[, "shuffled"] <- str_detect(allChIP_HOTincl_benchraptor_allmethods_plot$net, "_shuffle")

benchraptor_HOTincl_shufflestats_files <- shufflestats_files[str_detect(shufflestats_files, "RAPToR") &
                                                               str_detect(shufflestats_files, "HOTincl") &
                                                               str_detect(shufflestats_files, "allChIP")]

benchraptor_HOTincl_shufflestats_lists <- lapply(benchraptor_HOTincl_shufflestats_files, function(file){
  
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
  
  tempout["net"] <- paste0(str_extract(file, "[0-9]{3,5}"), "_shuffle")
  
  tempout
  
})

benchraptor_HOTincl_shufflestats_df <- do.call(rbind, benchraptor_HOTincl_shufflestats_lists)
benchraptor_HOTincl_shufflestats_df[, "cutoff"] <- str_remove(benchraptor_HOTincl_shufflestats_df$net, "_shuffle")

allChIP_HOTincl_benchraptor_allmethods_plot[allChIP_HOTincl_benchraptor_allmethods_plot$label == "5000", "label"] <- "none"

saveRDS(allChIP_HOTincl_benchraptor_allmethods_plot,
        "plotdata/allChIP_HOTincl_benchraptor_allmethods_plot.rds")
# allChIP_HOTincl_benchraptor_allmethods_plot <- readRDS("plotdata/allChIP_HOTincl_benchraptor_allmethods_plot.rds")

saveRDS(benchraptor_HOTincl_shufflestats_df,
        "plotdata/benchRAPToR_HOTincl_shufflestats_df.rds")
# benchraptor_HOTincl_shufflestats_df <- readRDS("plotdata/benchRAPToR_HOTincl_shufflestats_df.rds")

pdf("graphics/allChIP_HOTincl_benchRAPToR_methodsplot.pdf",
    height = 4,
    width = 5)

ggplot(data = allChIP_HOTincl_benchraptor_allmethods_plot[!str_detect(allChIP_HOTincl_benchraptor_allmethods_plot$net, "shuffle") & 
                                                            !allChIP_HOTincl_benchraptor_allmethods_plot$cutoff %in% c("1500", "3000", "10000") &
                                                            !allChIP_HOTincl_benchraptor_allmethods_plot$method %in% c("wsum_estimate", "wsum_corr"), ], aes(x = auroc, y = auprc, colour = cutoff)) + 
  geom_point(size = 2) +
  geom_label_repel(aes(label = label,
                       fill = cutoff),
                   colour = "black",
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
  geom_point(data = benchraptor_HOTincl_shufflestats_df[!benchraptor_HOTincl_shufflestats_df$cutoff %in% c("1500", "3000", "10000") &
                                                          !benchraptor_HOTincl_shufflestats_df$method %in% c("wsum_estimate", "wsum_corr"), ], aes(x = mean_auroc, y = mean_auprc), alpha = 0.5) + 
  geom_errorbar(data = benchraptor_HOTincl_shufflestats_df[!benchraptor_HOTincl_shufflestats_df$cutoff %in% c("1500", "3000", "10000") &
                                                             !benchraptor_HOTincl_shufflestats_df$method %in% c("wsum_estimate", "wsum_corr"), ], aes(x = mean_auroc, y = mean_auprc,
                                                                                                                                                      ymin = (mean_auprc - sd_auprc),
                                                                                                                                                      ymax = (mean_auprc + sd_auprc)), alpha = 0.5, width = 0.01) +
  geom_errorbar(data = benchraptor_HOTincl_shufflestats_df[!benchraptor_HOTincl_shufflestats_df$cutoff %in% c("1500", "3000", "10000") &
                                                             !benchraptor_HOTincl_shufflestats_df$method %in% c("wsum_estimate", "wsum_corr"), ], aes(x = mean_auroc, y = mean_auprc,
                                                                                                                                                      xmin = (mean_auroc - sd_auroc),
                                                                                                                                                      xmax = (mean_auroc + sd_auroc)), alpha = 0.5, width = 0.01)

dev.off()

pdf("graphics/allChIP_HOTincl_benchraptor_cutoffs_methodlines_AUPRC.pdf",
    height = 2.5,
    width = 2.5)

ggplot(data = allChIP_HOTincl_benchraptor_allmethods_plot[allChIP_HOTincl_benchraptor_allmethods_plot$shuffled == FALSE &
                                                            !allChIP_HOTincl_benchraptor_allmethods_plot$method %in% c("wsum_estimate", "wsum_corr") &
                                                            allChIP_HOTincl_benchraptor_allmethods_plot$cutoff != "10000", ], aes(x = as.numeric(cutoff), y = auprc, colour = method)) + 
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

dev.off()

pdf("graphics/allChIP_HOTincl_benchraptor_cutoffs_methodlines_AUROC.pdf",
    height = 2.5,
    width = 2.5)

ggplot(data = allChIP_HOTincl_benchraptor_allmethods_plot[allChIP_HOTincl_benchraptor_allmethods_plot$shuffled == FALSE &
                                                            !allChIP_HOTincl_benchraptor_allmethods_plot$method %in% c("wsum_estimate", "wsum_corr") &
                                                            allChIP_HOTincl_benchraptor_allmethods_plot$cutoff != "10000", ], aes(x = as.numeric(cutoff), y = auroc, colour = method)) + 
  geom_point() + 
  geom_line() +
  scale_x_log10(breaks = c(100, 500, 1000, 1500, 2000, 2500, 3000, 5000, 10000)) +
  ylab("Sensitivity (AUROC)") +
  xlab("Cutoff (log"[10]~"scale)") +
  coord_cartesian(ylim = c(0.5, 0.7)) +
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

dev.off()

#### benchraptor hot exclusion 

allChIP_HOTexcl_benchraptor_files <- bench_out_filelist[str_detect(bench_out_filelist, "benchRAPToR") &
                                                          str_detect(bench_out_filelist, "HOTexcl") &
                                                          str_detect(bench_out_filelist, "allChIP")]

allChIP_HOTexcl_benchraptor_files <- allChIP_HOTexcl_benchraptor_files[!str_detect(allChIP_HOTexcl_benchraptor_files, "dense")]
allChIP_HOTexcl_benchraptor_files <- allChIP_HOTexcl_benchraptor_files[!str_detect(allChIP_HOTexcl_benchraptor_files, "shufflestat")]

allChIP_HOTexcl_benchraptor_list <- lapply(allChIP_HOTexcl_benchraptor_files, function(x){
  
  read.table(paste0("output/benchmark_out/", x),
             sep = "\t",
             header = TRUE)
  
})

allChIP_HOTexcl_benchraptor_all <- do.call(rbind, allChIP_HOTexcl_benchraptor_list)

allChIP_HOTexcl_benchraptor_allmethods_plot <- allChIP_HOTexcl_benchraptor_all[allChIP_HOTexcl_benchraptor_all$metric %in% c("auroc", "auprc"), c("net", "method", "metric", "score")]
allChIP_HOTexcl_benchraptor_allmethods_plot <- tidyr::pivot_wider(allChIP_HOTexcl_benchraptor_allmethods_plot, names_from = c("metric"), values_from = "score")

allChIP_HOTexcl_benchraptor_allmethods_plot[, "label"] <- as.character(allChIP_HOTexcl_benchraptor_allmethods_plot$net)
# allChIP_HOTexcl_benchraptor_allmethods_plot[str_detect(unlist(allChIP_HOTexcl_benchraptor_allmethods_plot$net), "shuffle"), "label"] <- ""

allChIP_HOTexcl_benchraptor_allmethods_plot[, "cutoff"] <- str_remove(allChIP_HOTexcl_benchraptor_allmethods_plot$net, "_shuffle")
allChIP_HOTexcl_benchraptor_allmethods_plot[, "shuffled"] <- str_detect(allChIP_HOTexcl_benchraptor_allmethods_plot$net, "_shuffle")

allChIP_HOTexcl_benchraptor_allmethods_plot[allChIP_HOTexcl_benchraptor_allmethods_plot$label == "5000", "label"] <- "none"

benchraptor_HOTexcl_shufflestats_files <- shufflestats_files[str_detect(shufflestats_files, "RAPToR") &
                                                               str_detect(shufflestats_files, "HOTexcl") &
                                                               str_detect(shufflestats_files, "allChIP")]

benchraptor_HOTexcl_shufflestats_lists <- lapply(benchraptor_HOTexcl_shufflestats_files, function(file){
  
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
  
  tempout["net"] <- paste0(str_extract(file, "[0-9]{3,5}"), "_shuffle")
  
  tempout
  
})

benchraptor_HOTexcl_shufflestats_df <- do.call(rbind, benchraptor_HOTexcl_shufflestats_lists)
benchraptor_HOTexcl_shufflestats_df[, "cutoff"] <- str_remove(benchraptor_HOTexcl_shufflestats_df$net, "_shuffle")

saveRDS(allChIP_HOTexcl_benchraptor_allmethods_plot,
        "plotdata/allChIP_HOTexcl_benchraptor_allmethods_plot.rds")
# allChIP_HOTexcl_benchraptor_allmethods_plot <- readRDS("plotdata/allChIP_HOTexcl_benchraptor_allmethods_plot.rds")

saveRDS(benchraptor_HOTexcl_shufflestats_df,
        "plotdata/benchRAPToR_HOTexcl_shufflestats_df.rds")
# benchraptor_HOTexcl_shufflestats_df <- readRDS("plotdata/benchRAPToR_HOTexcl_shufflestats_df.rds")

pdf("graphics/allChIP_HOTexcl_benchRAPToR_mlmplot.pdf",
    height = 2.5,
    width = 2.5)

ggplot(data = allChIP_HOTexcl_benchraptor_allmethods_plot[!str_detect(allChIP_HOTexcl_benchraptor_allmethods_plot$net, "shuffle") & 
                                                            !allChIP_HOTexcl_benchraptor_allmethods_plot$cutoff %in% c("1500", "3000", "10000") &
                                                            allChIP_HOTexcl_benchraptor_allmethods_plot$method == "mlm_estimate", ], aes(x = auroc, y = auprc, colour = cutoff)) + 
  geom_point(size = 2) +
  geom_label_repel(aes(label = label,
                       fill = cutoff),
                   colour = "black",
                   label.size = 0.1,
                   max.overlaps = 10,
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
  geom_point(data = benchraptor_HOTexcl_shufflestats_df[!benchraptor_HOTexcl_shufflestats_df$cutoff %in% c("1500", "3000", "10000") &
                                                          benchraptor_HOTexcl_shufflestats_df$method == "mlm_estimate", ], aes(x = mean_auroc, y = mean_auprc), alpha = 0.5) + 
  geom_errorbar(data = benchraptor_HOTexcl_shufflestats_df[!benchraptor_HOTexcl_shufflestats_df$cutoff %in% c("1500", "3000", "10000") &
                                                             benchraptor_HOTexcl_shufflestats_df$method == "mlm_estimate", ], aes(x = mean_auroc, y = mean_auprc,
                                                                                                                                  ymin = (mean_auprc - sd_auprc),
                                                                                                                                  ymax = (mean_auprc + sd_auprc)), alpha = 0.5, width = 0.01) +
  geom_errorbar(data = benchraptor_HOTexcl_shufflestats_df[!benchraptor_HOTexcl_shufflestats_df$cutoff %in% c("1500", "3000", "10000") &
                                                             benchraptor_HOTexcl_shufflestats_df$method == "mlm_estimate", ], aes(x = mean_auroc, y = mean_auprc,
                                                                                                                                  xmin = (mean_auroc - sd_auroc),
                                                                                                                                  xmax = (mean_auroc + sd_auroc)), alpha = 0.5, width = 0.01)
dev.off()

pdf("graphics/allChIP_HOTexcl_benchraptor_cutoffs_methodlines_AUPRC.pdf",
    height = 2,
    width = 2)

ggplot(data = allChIP_HOTexcl_benchraptor_allmethods_plot[allChIP_HOTexcl_benchraptor_allmethods_plot$shuffled == FALSE &
                                                            !allChIP_HOTexcl_benchraptor_allmethods_plot$method %in% c("wsum_estimate", "wsum_corr") &
                                                            allChIP_HOTexcl_benchraptor_allmethods_plot$cutoff != "10000", ], aes(x = as.numeric(cutoff), y = auprc, colour = method)) + 
  geom_point() + 
  geom_line() +
  scale_x_log10(breaks = c(100, 500, 1000, 1500, 2000, 2500, 3000, 5000, 10000)) +
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
                                   size = 6),
        axis.title = element_text(size = 10))

dev.off()

pdf("graphics/allChIP_HOTexcl_benchraptor_cutoffs_methodlines_AUROC.pdf",
    height = 2,
    width = 2)

ggplot(data = allChIP_HOTexcl_benchraptor_allmethods_plot[allChIP_HOTexcl_benchraptor_allmethods_plot$shuffled == FALSE &
                                                            !allChIP_HOTexcl_benchraptor_allmethods_plot$method %in% c("wsum_estimate", "wsum_corr") &
                                                            allChIP_HOTexcl_benchraptor_allmethods_plot$cutoff != "10000", ], aes(x = as.numeric(cutoff), y = auroc, colour = method)) + 
  geom_point() + 
  geom_line() +
  scale_x_log10(breaks = c(100, 500, 1000, 1500, 2000, 2500, 3000, 5000, 10000)) +
  ylab("Sensitivity (AUROC)") +
  xlab("Cutoff (log"[10]~"scale)") +
  coord_cartesian(ylim = c(0.5, 0.7)) +
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
                                   size = 6),
        axis.title = element_text(size = 10))

dev.off()

pdf("graphics/allChIP_HOTexcl_benchRAPToR_methodsplot.pdf",
    height = 4,
    width = 4)

ggplot(data = allChIP_HOTexcl_benchraptor_allmethods_plot[!str_detect(allChIP_HOTexcl_benchraptor_allmethods_plot$net, "shuffle") & 
                                                            !allChIP_HOTexcl_benchraptor_allmethods_plot$cutoff %in% c("1500", "3000", "10000") &
                                                            !allChIP_HOTexcl_benchraptor_allmethods_plot$method %in% c("wsum_estimate", "wsum_corr"), ], aes(x = auroc, y = auprc, colour = cutoff)) + 
  geom_point(size = 2) +
  geom_label_repel(aes(label = label,
                       fill  = cutoff),
                   colour = "black",
                   label.size = 0.1,
                   max.overlaps = 15,
                   box.padding = 0.3,
                   label.padding = 0.1,
                   point.padding = 0,
                   size = 2
  ) +
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
        axis.text = element_text(colour = "black"),
        axis.title = element_text(size = 10)) +
  geom_point(data = benchraptor_HOTexcl_shufflestats_df[!benchraptor_HOTexcl_shufflestats_df$cutoff %in% c("1500", "3000", "10000") &
                                                          !benchraptor_HOTexcl_shufflestats_df$method %in% c("wsum_estimate", "wsum_corr"), ], aes(x = mean_auroc, y = mean_auprc), alpha = 0.5) + 
  geom_errorbar(data = benchraptor_HOTexcl_shufflestats_df[!benchraptor_HOTexcl_shufflestats_df$cutoff %in% c("1500", "3000", "10000") &
                                                             !benchraptor_HOTexcl_shufflestats_df$method %in% c("wsum_estimate", "wsum_corr"), ], aes(x = mean_auroc, y = mean_auprc,
                                                                                                                                                      ymin = (mean_auprc - sd_auprc),
                                                                                                                                                      ymax = (mean_auprc + sd_auprc)), alpha = 0.5, width = 0.01) +
  geom_errorbar(data = benchraptor_HOTexcl_shufflestats_df[!benchraptor_HOTexcl_shufflestats_df$cutoff %in% c("1500", "3000", "10000") &
                                                             !benchraptor_HOTexcl_shufflestats_df$method %in% c("wsum_estimate", "wsum_corr"), ], aes(x = mean_auroc, y = mean_auprc,
                                                                                                                                                      xmin = (mean_auroc - sd_auroc),
                                                                                                                                                      xmax = (mean_auroc + sd_auroc)), alpha = 0.5, width = 0.01)

dev.off()

#### CHIP HOT OR NOT MLM LINE PLOT ####
hotornot_plot <- rbind(allChIP_HOTexcl_benchraptor_allmethods_plot[allChIP_HOTexcl_benchraptor_allmethods_plot$method == "mlm_estimate", ],
                       allChIP_HOTincl_benchraptor_allmethods_plot[allChIP_HOTincl_benchraptor_allmethods_plot$method == "mlm_estimate", ])


hotornot_plot[1: nrow(allChIP_HOTexcl_benchraptor_allmethods_plot[allChIP_HOTexcl_benchraptor_allmethods_plot$method == "mlm_estimate", ]), "group"] <- "HOTexcl"
hotornot_plot[(nrow(allChIP_HOTexcl_benchraptor_allmethods_plot[allChIP_HOTexcl_benchraptor_allmethods_plot$method == "mlm_estimate", ]) + 1):nrow(hotornot_plot), "group"] <- "HOTincl"

hotornot_shuffleplot <- rbind(benchraptor_HOTexcl_shufflestats_df[benchraptor_HOTexcl_shufflestats_df$method == "mlm_estimate", ],
                              benchraptor_HOTincl_shufflestats_df[benchraptor_HOTincl_shufflestats_df$method == "mlm_estimate", ])

hotornot_shuffleplot[1: nrow(benchraptor_HOTexcl_shufflestats_df[benchraptor_HOTexcl_shufflestats_df$method == "mlm_estimate", ]), "group"] <- "HOTexcl"
hotornot_shuffleplot[(nrow(benchraptor_HOTexcl_shufflestats_df[benchraptor_HOTexcl_shufflestats_df$method == "mlm_estimate", ]) + 1):nrow(hotornot_shuffleplot), "group"] <- "HOTincl"

saveRDS(hotornot_shuffleplot, "plotdata/hotornot_plot")

pdf("graphics/Hotornot_mlm_plot.pdf",
    width = 2.5,
    height = 2.5)

ggplot(data = hotornot_plot, aes(x = as.numeric(cutoff), y = auprc, colour = group)) + 
  geom_point() + 
  geom_line() +
  scale_x_log10(breaks = c(100, 500, 1000, 1500, 2000, 2500, 3000, 5000, 10000)) +
  ylab("Precision (AUPRC)") +
  xlab("Max. targets cutoff (log"[10]~"scale)") +
  coord_cartesian(ylim = c(0.45, 0.8)) +
  theme(strip.text = element_text(colour = 'white'),
        strip.background = element_rect(fill = "black"),
        legend.position = "none",
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_line(colour = "grey", linewidth = 0.2),
        # panel.grid.minor = element_line(colour = "grey", linewidth = 0.1),
        axis.text = element_text(colour = "black"),
        axis.title = element_text(size = 10),
        axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.5,
                                   size = 6) 
  ) +
  
  geom_point(data = hotornot_shuffleplot, aes(x = as.numeric(cutoff), y = mean_auprc, colour = group),
             alpha = 0.3) +
  geom_errorbar(data = hotornot_shuffleplot, aes(x = as.numeric(cutoff), y = mean_auprc,
                                                                                                                                        ymin = (mean_auprc - sd_auprc),
                                                                                                                                        ymax = (mean_auprc + sd_auprc),
                                                 colour = group), alpha = 0.3                ) +
  geom_line(data = hotornot_shuffleplot, aes(x = as.numeric(cutoff), y = mean_auprc, colour = group), alpha = 0.3)

dev.off()

#### CHIP HOT CUTOFF PLOTS ####

hotcutfiles <- list.files("output/GRNs")[str_detect(list.files("output/GRNs"), "HOTcut")]
# hotcutfiles <- list.files("~/Cel_GRN_manuscript/output/GRNs")[str_detect(list.files("~/Cel_GRN_manuscript/output/GRNs"), "HOTcut")]

HOTcut_benchraptor_files <- bench_out_filelist[str_detect(bench_out_filelist, "HotCut")]

HOTcut_TF_numbers <- sapply(hotcutfiles, function(thisfile){
  
  Hotcut <- str_extract(thisfile, "HOTcut[0-9]{2}")
  
  thisGRN <- read.table(paste0("~/Cel_GRN_manuscript/output/GRNs/allChIP_10000_", Hotcut, ".txt"),
                        header = TRUE)
  
  length(unique(thisGRN$source))
  
})

names(HOTcut_TF_numbers) <- str_remove(str_remove(names(HOTcut_TF_numbers), "allChIP_10000_HOTcut"), "\\.txt")

HOTcut_benchraptor_list <- lapply(HOTcut_benchraptor_files, function(x){
  
  read.table(paste0("output/benchmark_out/", x),
             sep = "\t",
             header = TRUE)
  
})

HOTcut_benchraptor_all <- do.call(rbind, HOTcut_benchraptor_list)

HOTcut_benchraptor_allmethods_plot <- HOTcut_benchraptor_all[HOTcut_benchraptor_all$metric %in% c("auroc", "auprc"), c("net", "method", "metric", "score")]
HOTcut_benchraptor_allmethods_plot <- tidyr::pivot_wider(HOTcut_benchraptor_allmethods_plot, names_from = c("metric"), values_from = "score")

HOTcut_benchraptor_allmethods_plot[, "label"] <- as.character(HOTcut_benchraptor_allmethods_plot$net)

pdf("graphics/allChIP_HOTcut__benchraptor_cutoffs_methodlines_AUPRC.pdf",
    height = 2,
    width = 2)

ggplot(data = HOTcut_benchraptor_allmethods_plot[!HOTcut_benchraptor_allmethods_plot$method %in% c("wsum_estimate", "wsum_corr"), ], aes(x = net, y = auprc, colour = method)) + 
  geom_point() + 
  geom_line() +
  scale_x_continuous(breaks = c(20, 30, 40, 50, 70)) +
  ylab("Precision (AUPRC)") +
  xlab("HOT regions max.\nTFs cutoff") +
  coord_cartesian(ylim = c(0.5, 0.75),
                  xlim = c(18, 72)) +
  theme(strip.text = element_text(colour = 'white'),
        strip.background = element_rect(fill = "black"),
        legend.position = "none",
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_line(colour = "grey", linewidth = 0.2),
        # panel.grid.minor = element_line(colour = "grey", linewidth = 0.1),
        axis.text = element_text(colour = "black"),
        axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.5),
        axis.title = element_text(size = 10) ) 
# + 
  # annotate(geom = "text",
  #          x = 45,
  #          y = 0.575,
  #          label = "TFs in GRN:",
  #          fontface = "bold",
  #          cex = 2.5) +
  # annotate(geom = "text",
  #          x = as.numeric(names(HOTcut_TF_numbers)),
  #          y = 0.525,
  #          label = HOTcut_TF_numbers,
  #          cex = 2)

dev.off()

pdf("graphics/allChIP_HOTcut_methodlines_AUROC.pdf",
    height = 2,
    width = 2)

ggplot(data = HOTcut_benchraptor_allmethods_plot[!HOTcut_benchraptor_allmethods_plot$method %in% c("wsum_estimate", "wsum_corr"), ], aes(x = net, y = auroc, colour = method)) + 
  geom_point() + 
  geom_line() +
  scale_x_continuous(breaks = c(20, 30, 40, 50, 70)) +
  ylab("Sensitivity (AUROC)") +
  xlab("HOT regions max.\nTFs cutoff") +
  coord_cartesian(ylim = c(0.5, 0.7),
                  xlim = c(18, 72)) +
  theme(strip.text = element_text(colour = 'white'),
        strip.background = element_rect(fill = "black"),
        legend.position = "none",
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_line(colour = "grey", linewidth = 0.2),
        # panel.grid.minor = element_line(colour = "grey", linewidth = 0.1),
        axis.text = element_text(colour = "black"),
        axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.5),
        axis.title = element_text(size = 10)) +
  annotate(geom = "text",
           x = 45,
           y = 0.575,
           label = "TFs in GRN:",
           fontface = "bold",
           cex = 2.5) +
  annotate(geom = "text",
           x = as.numeric(names(HOTcut_TF_numbers)),
           y = 0.525,
           label = HOTcut_TF_numbers,
           cex = 2)

dev.off()

#### FIMO plots ####

FIMO_nohomo_benchraptor_files <- bench_out_filelist[str_detect(bench_out_filelist, "benchRAPToR") &
                                                      str_detect(bench_out_filelist, "nohomo")]

FIMO_nohomo_benchraptor_files <- FIMO_nohomo_benchraptor_files[!str_detect(FIMO_nohomo_benchraptor_files, "shufflestat")]

FIMO_nohomo_benchraptor_list <- lapply(FIMO_nohomo_benchraptor_files, function(x){
  
  read.table(paste0("output/benchmark_out/", x),
             sep = "\t",
             header = TRUE)
  
})

FIMO_nohomo_benchraptor_all <- do.call(rbind, FIMO_nohomo_benchraptor_list)

FIMO_nohomo_benchraptor_allmethods_plot <- FIMO_nohomo_benchraptor_all[FIMO_nohomo_benchraptor_all$metric %in% c("auroc", "auprc"), c("net", "method", "metric", "score")]
FIMO_nohomo_benchraptor_allmethods_plot <- tidyr::pivot_wider(FIMO_nohomo_benchraptor_allmethods_plot, names_from = c("metric"), values_from = "score")

FIMO_nohomo_benchraptor_allmethods_plot[, "label"] <- as.character(FIMO_nohomo_benchraptor_allmethods_plot$net)
FIMO_nohomo_benchraptor_allmethods_plot[str_detect(unlist(FIMO_nohomo_benchraptor_allmethods_plot$net), "shuffle"), "label"] <- ""

FIMO_nohomo_benchraptor_allmethods_plot[, "cutoff"] <- str_remove(FIMO_nohomo_benchraptor_allmethods_plot$net, "_shuffle")
FIMO_nohomo_benchraptor_allmethods_plot[, "shuffled"] <- str_detect(FIMO_nohomo_benchraptor_allmethods_plot$net, "_shuffle")

FIMO_nohomo_benchraptor_shufflestats_files <- shufflestats_files[str_detect(shufflestats_files, "RAPToR") & str_detect(shufflestats_files, "HOTincl")]

FIMO_nohomo_benchraptor_shufflestats_lists <- lapply(FIMO_nohomo_benchraptor_shufflestats_files, function(file){

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
  
  tempout["net"] <- paste0(str_extract(file, "[0-9]{3,5}"), "_shuffle")
  
  tempout
  
})

FIMO_nohomo_benchraptor_shufflestats_df <- do.call(rbind, FIMO_nohomo_benchraptor_shufflestats_lists)
FIMO_nohomo_benchraptor_shufflestats_df[, "cutoff"] <- str_remove(FIMO_nohomo_benchraptor_shufflestats_df$net, "_shuffle")

method_labels <- c("Consensus", "MLM", "ULM", "WSum")
names(method_labels) <- c("consensus_estimate", "mlm_estimate", "ulm_estimate", "wsum_norm")

saveRDS(FIMO_nohomo_benchraptor_allmethods_plot,
        "plotdata/FIMO_nohomo_benchraptor_allmethods_plot.rds")
# FIMO_nohomo_benchraptor_allmethods_plot <- readRDS("plotdata/FIMO_nohomo_benchraptor_allmethods_plot.rds")

saveRDS(FIMO_nohomo_benchraptor_shufflestats_df,
        "plotdata/FIMO_nohomo_benchraptor_shufflestats_df.rds")
# FIMO_nohomo_benchraptor_shufflestats_df <- readRDS("plotdata/FIMO_nohomo_benchraptor_shufflestats_df.rds")

pdf("graphics/FIMO_nohomo_benchraptor_methodsplot.pdf",
    height = 4,
    width = 5)

ggplot(data = FIMO_nohomo_benchraptor_allmethods_plot[!str_detect(FIMO_nohomo_benchraptor_allmethods_plot$net, "shuffle") & 
                                                        !FIMO_nohomo_benchraptor_allmethods_plot$cutoff %in% c("2500", "10000") &
                                                        !FIMO_nohomo_benchraptor_allmethods_plot$method %in% c("wsum_estimate", "wsum_corr"), ], aes(x = auroc, y = auprc, colour = cutoff)) + 
  geom_point(size = 2) +
  geom_label_repel(aes(label = label,
                       fill = cutoff),
                   colour = "black",
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
  geom_point(data = FIMO_nohomo_benchraptor_shufflestats_df[!FIMO_nohomo_benchraptor_shufflestats_df$cutoff %in% c("2500", "10000") &
                                                              !FIMO_nohomo_benchraptor_shufflestats_df$method %in% c("wsum_estimate", "wsum_corr"), ], aes(x = mean_auroc, y = mean_auprc), alpha = 0.5) + 
  geom_errorbar(data = FIMO_nohomo_benchraptor_shufflestats_df[!FIMO_nohomo_benchraptor_shufflestats_df$cutoff %in% c("2500", "10000") &
                                                                 !FIMO_nohomo_benchraptor_shufflestats_df$method %in% c("wsum_estimate", "wsum_corr"), ], aes(x = mean_auroc, y = mean_auprc,
                                                                                                                                                              ymin = (mean_auprc - sd_auprc),
                                                                                                                                                              ymax = (mean_auprc + sd_auprc)), alpha = 0.5, width = 0.01) +
  geom_errorbar(data = FIMO_nohomo_benchraptor_shufflestats_df[!FIMO_nohomo_benchraptor_shufflestats_df$cutoff %in% c("2500", "10000") &
                                                                 !FIMO_nohomo_benchraptor_shufflestats_df$method %in% c("wsum_estimate", "wsum_corr"), ], aes(x = mean_auroc, y = mean_auprc,
                                                                                                                                                              xmin = (mean_auroc - sd_auroc),
                                                                                                                                                              xmax = (mean_auroc + sd_auroc)), alpha = 0.5, width = 0.01)

dev.off()

View(FIMO_nohomo_benchraptor_allmethods_plot[FIMO_nohomo_benchraptor_allmethods_plot$shuffled == FALSE & FIMO_nohomo_benchraptor_allmethods_plot$method == "mlm_estimate", ])

pdf("graphics/FIMO_nohomo_benchraptor_mlmplot.pdf",
    height = 2.5,
    width = 2.5)

ggplot(data = FIMO_nohomo_benchraptor_allmethods_plot[!str_detect(FIMO_nohomo_benchraptor_allmethods_plot$net, "shuffle") & 
                                                        !FIMO_nohomo_benchraptor_allmethods_plot$cutoff %in% c("2500", "10000") &
                                                        FIMO_nohomo_benchraptor_allmethods_plot$method == "mlm_estimate", ], aes(x = auroc, y = auprc, colour = cutoff)) + 
  geom_point(size = 2) +
  geom_label_repel(aes(label = label),
                   label.size = 0.1,
                   max.overlaps = 10,
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
  geom_point(data = FIMO_nohomo_benchraptor_shufflestats_df[!FIMO_nohomo_benchraptor_shufflestats_df$cutoff %in% c("2500", "10000") &
                                                              FIMO_nohomo_benchraptor_shufflestats_df$method == "mlm_estimate", ], aes(x = mean_auroc, y = mean_auprc), alpha = 0.5) + 
  geom_errorbar(data = FIMO_nohomo_benchraptor_shufflestats_df[!FIMO_nohomo_benchraptor_shufflestats_df$cutoff %in% c("2500", "10000") &
                                                                 FIMO_nohomo_benchraptor_shufflestats_df$method == "mlm_estimate", ], aes(x = mean_auroc, y = mean_auprc,
                                                                                                                                          ymin = (mean_auprc - sd_auprc),
                                                                                                                                          ymax = (mean_auprc + sd_auprc)), alpha = 0.5) +
  geom_errorbar(data = FIMO_nohomo_benchraptor_shufflestats_df[!FIMO_nohomo_benchraptor_shufflestats_df$cutoff %in% c("2500", "10000") &
                                                                 FIMO_nohomo_benchraptor_shufflestats_df$method == "mlm_estimate", ], aes(x = mean_auroc, y = mean_auprc,
                                                                                                                                          xmin = (mean_auroc - sd_auroc),
                                                                                                                                          xmax = (mean_auroc + sd_auroc)), alpha = 0.5)
dev.off()

pdf("graphics/FIMO_benchraptor_cutoffs_methodlines_AUPRC_legend.pdf",
    height = 5,
    width = 5)

ggplot(data = FIMO_nohomo_benchraptor_allmethods_plot[FIMO_nohomo_benchraptor_allmethods_plot$shuffled == FALSE &
                                                        !FIMO_nohomo_benchraptor_allmethods_plot$method %in% c("wsum_estimate", 'wsum_corr'), ], aes(x = as.numeric(cutoff), y = auprc, colour = method)) + 
  geom_point() + 
  geom_line() +
  scale_x_log10(breaks = c(100, 500, 1000, 1500, 2000, 2500, 3000, 5000, 10000)) +
  ylab("Precision (AUPRC)") +
  xlab("Cutoff (log"[10]~"scale)") +
  coord_cartesian(ylim = c(0.55, 0.7)) +
  theme(strip.text = element_text(colour = 'white'),
        strip.background = element_rect(fill = "black"),
        # legend.position = "none",
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_line(colour = "grey", linewidth = 0.2),
        # panel.grid.minor = element_line(colour = "grey", linewidth = 0.1),
        axis.text = element_text(colour = "black"),
        axis.text.x = element_text(angle = 30,
                                   hjust = 0.5) 
  )

dev.off()

pdf("graphics/FIMO_benchraptor_cutoffs_methodlines_AUPRC.pdf",
    height = 2.5,
    width = 2.5)

ggplot(data = FIMO_nohomo_benchraptor_allmethods_plot[FIMO_nohomo_benchraptor_allmethods_plot$shuffled == FALSE &
                                                        !FIMO_nohomo_benchraptor_allmethods_plot$method %in% c("wsum_estimate", 'wsum_corr'), ], aes(x = as.numeric(cutoff), y = auprc, colour = method)) + 
  geom_point() + 
  geom_line() +
  scale_x_log10(breaks = c(100, 500, 1000, 1500, 2000, 2500, 3000, 5000, 10000)) +
  ylab("Precision (AUPRC)") +
  xlab("Cutoff (log"[10]~"scale)") +
  coord_cartesian(ylim = c(0.55, 0.7)) +
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
                                   size = 6) 
  )

dev.off()

pdf("graphics/FIMO_benchraptor_cutoffs_methodlines_AUROC.pdf",
    height = 2.5,
    width = 2.5)

ggplot(data = FIMO_nohomo_benchraptor_allmethods_plot[FIMO_nohomo_benchraptor_allmethods_plot$shuffled == FALSE &
                                                        !FIMO_nohomo_benchraptor_allmethods_plot$method %in% c("wsum_estimate", 'wsum_corr'), ], aes(x = as.numeric(cutoff), y = auroc, colour = method)) + 
  geom_point() + 
  geom_line() +
  scale_x_log10(breaks = c(100, 500, 1000, 1500, 2000, 2500, 3000, 5000, 10000)) +
  ylab("Sensitivity (AUROC)") +
  xlab("Cutoff (log"[10]~"scale)") +
  coord_cartesian(ylim = c(0.55, 0.7)) +
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
                                   size = 6) 
  )

dev.off()

## plots with MLM lines and random ##

pdf("graphics/FIMO_benchraptor_cutoffs_mlmAUPRC_line.pdf",
    height = 2.5,
    width = 2.5)

ggplot(data = FIMO_nohomo_benchraptor_allmethods_plot[FIMO_nohomo_benchraptor_allmethods_plot$shuffled == FALSE &
                                                        FIMO_nohomo_benchraptor_allmethods_plot$method %in% c('mlm_estimate'), ], aes(x = as.numeric(cutoff), y = auprc, colour = method)) + 
  geom_point() + 
  geom_line() +
  scale_x_log10(breaks = c(100, 500, 1000, 1500, 2000, 2500, 3000, 5000, 10000)) +
  ylab("Precision (AUPRC)") +
  xlab("Max. targets cutoff (log"[10]~"scale)") +
  coord_cartesian(ylim = c(0.45, 0.7)) +
  theme(strip.text = element_text(colour = 'white'),
        strip.background = element_rect(fill = "black"),
        legend.position = "none",
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_line(colour = "grey", linewidth = 0.2),
        # panel.grid.minor = element_line(colour = "grey", linewidth = 0.1),
        axis.text = element_text(colour = "black"),
        axis.title = element_text(size = 10),
        axis.text.x = element_text(angle = 90,
                                   hjust = 1,
                                   vjust = 0.5,
                                   size = 6) 
  ) +

geom_point(data = FIMO_nohomo_benchraptor_shufflestats_df[FIMO_nohomo_benchraptor_shufflestats_df$method == "mlm_estimate", ], aes(x = as.numeric(cutoff), y = mean_auprc, colour = method),
           alpha = 0.3) + 
  geom_errorbar(data = FIMO_nohomo_benchraptor_shufflestats_df[FIMO_nohomo_benchraptor_shufflestats_df$method == "mlm_estimate", ], aes(x = as.numeric(cutoff), y = mean_auprc,
                                                                                                                                          ymin = (mean_auprc - sd_auprc),
                                                                                                                                          ymax = (mean_auprc + sd_auprc),
                                                                                                                                        colour = method), alpha = 0.5) +
  geom_line(data = FIMO_nohomo_benchraptor_shufflestats_df[FIMO_nohomo_benchraptor_shufflestats_df$method == "mlm_estimate", ], aes(x = as.numeric(cutoff), y = mean_auprc, colour = method), alpha = 0.3, )
                                                                                                                                          

dev.off()

#### FIMO HOMOTYPIC BINDING PLOTS ####

#### homotypic 

param5 <- 5

FIMO_param5_benchraptor_files <- bench_out_filelist[str_detect(bench_out_filelist, "benchRAPToR") &
                                                      str_detect(bench_out_filelist, paste0("param", param5))]

FIMO_param5_benchraptor_files <- FIMO_param5_benchraptor_files[!str_detect(FIMO_param5_benchraptor_files, "shufflestat")]

FIMO_param5_benchraptor_list <- lapply(FIMO_param5_benchraptor_files, function(x){
  
  read.table(paste0("output/benchmark_out/", x),
             sep = "\t",
             header = TRUE)
  
})

FIMO_param5_benchraptor_all <- do.call(rbind, FIMO_param5_benchraptor_list)

FIMO_param5_benchraptor_allmethods_plot <- FIMO_param5_benchraptor_all[FIMO_param5_benchraptor_all$metric %in% c("auroc", "auprc"), c("net", "method", "metric", "score")]
FIMO_param5_benchraptor_allmethods_plot <- tidyr::pivot_wider(FIMO_param5_benchraptor_allmethods_plot, names_from = c("metric"), values_from = "score")

FIMO_param5_benchraptor_allmethods_plot[, "label"] <- FIMO_param5_benchraptor_allmethods_plot$net
FIMO_param5_benchraptor_allmethods_plot[str_detect(unlist(FIMO_param5_benchraptor_allmethods_plot$net), "shuffle"), "label"] <- ""

FIMO_param5_benchraptor_allmethods_plot[, "cutoff"] <- str_remove(FIMO_param5_benchraptor_allmethods_plot$net, "_shuffle")
FIMO_param5_benchraptor_allmethods_plot[, "shuffled"] <- str_detect(FIMO_param5_benchraptor_allmethods_plot$net, "_shuffle")

FIMO_param5_benchraptor_allmethods_plot[, "param"] <- as.character(param5)

# add no homotypic binding df
FIMO_noparam <- FIMO_nohomo_benchraptor_allmethods_plot
FIMO_noparam[, "param"] <- "none"

FIMO_param5_benchraptor_allmethods_plot <- rbind(FIMO_param5_benchraptor_allmethods_plot, FIMO_noparam)

FIMO_param5_benchraptor_allmethods_plot[FIMO_param5_benchraptor_allmethods_plot$param == "5", "param"] <- "+ homotypic binding"

FIMO_param5_benchraptor_shufflestats_files <- shufflestats_files[str_detect(shufflestats_files, "RAPToR") & str_detect(shufflestats_files, "homoparam5")]

FIMO_param5_benchraptor_shufflestats_lists <- lapply(FIMO_param5_benchraptor_shufflestats_files, function(file){
  
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
  
  tempout["net"] <- paste0(str_extract(file, "[0-9]{3,4}"), "_shuffle")
  
  tempout
  
})

FIMO_param5_benchraptor_shufflestats_df <- do.call(rbind, FIMO_param5_benchraptor_shufflestats_lists)
FIMO_param5_benchraptor_shufflestats_df[, "cutoff"] <- str_remove(FIMO_param5_benchraptor_shufflestats_df$net, "_shuffle")
FIMO_param5_benchraptor_shufflestats_df[, "cutoff"] <- paste0("cutoff: ", FIMO_param5_benchraptor_shufflestats_df[, "cutoff"])
FIMO_param5_benchraptor_shufflestats_df$cutoff <- factor(FIMO_param5_benchraptor_shufflestats_df$cutoff, levels = paste0("cutoff: ", c("100", "500", "1000", "2000")))

FIMO_param5_benchraptor_allmethods_plot_df <- FIMO_param5_benchraptor_allmethods_plot[FIMO_param5_benchraptor_allmethods_plot$method == "mlm_estimate" &
                                                                                        FIMO_param5_benchraptor_allmethods_plot$cutoff %in% c("100", "500", "1000", "2000") &
                                                                                        !FIMO_param5_benchraptor_allmethods_plot$shuffled, ]

FIMO_param5_benchraptor_allmethods_plot_df$cutoff <- paste0("cutoff: ", FIMO_param5_benchraptor_allmethods_plot_df$cutoff)
FIMO_param5_benchraptor_allmethods_plot_df$cutoff <- factor(FIMO_param5_benchraptor_allmethods_plot_df$cutoff, levels = paste0("cutoff: ", c("100", "500", "1000", "2000")))

saveRDS(FIMO_param5_benchraptor_allmethods_plot_df, 
        "plotdata/FIMO_param5_benchraptor_allmethods_plot_df.rds")

saveRDS(FIMO_param5_benchraptor_shufflestats_df,
        "plotdata/FIMO_param5_benchraptor_shufflestats_df.rds")

pdf("graphics/FIMO_homotypic_methodsplot.pdf",
    height = 4, 
    width = 5)

ggplot(data = FIMO_param5_benchraptor_allmethods_plot_df, aes(x = auroc, y = auprc, alpha = shuffled)) + 
  geom_point(aes(color = param)) + 
  geom_label_repel(aes(label = stringr::str_wrap(param, 11),
                       fill = param),
                   colour = "black",
                   label.size = 0.1,
                   max.overlaps = 50,
                   box.padding = 0.1,
                   label.padding = 0.1,
                   point.padding = 1,
                   size = 3) +
  scale_colour_manual(values = c("red", "grey")) + 
  scale_fill_manual(values = c("red", "grey")) +
  scale_alpha_discrete(range = c(1, 0.2)) +
  theme_classic() + 
  geom_vline(xintercept = 0.5,
             linetype = "dashed",
             col = "grey") +
  geom_hline(yintercept = 0.5,
             linetype = "dashed",
             col = "grey") + 
  facet_wrap(~ cutoff) +
  xlab("Sensitivity (AUROC)") + 
  ylab("Precision (AUROC)") +
  coord_cartesian(xlim = c(0.45, 0.7),
                  ylim = c(0.45, 0.7)) +
  theme(strip.text = element_text(colour = 'white'),
        strip.background = element_rect(fill = "black"),
        legend.position = "none",
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_line(colour = "grey", linewidth = 0.2),
        panel.grid.minor = element_line(colour = "grey", linewidth = 0.1),
        axis.text = element_text(colour = "black")) + 
  geom_point(data = FIMO_param5_benchraptor_shufflestats_df[FIMO_param5_benchraptor_shufflestats_df$method == "mlm_estimate", ], aes(x = mean_auroc, y = mean_auprc), alpha = 0.5) + 
  geom_errorbar(data = FIMO_param5_benchraptor_shufflestats_df[FIMO_param5_benchraptor_shufflestats_df$method == "mlm_estimate", ], aes(x = mean_auroc, y = mean_auprc,
                                                                                                                                        ymin = (mean_auprc - sd_auprc),
                                                                                                                                        ymax = (mean_auprc + sd_auprc)), alpha = 0.5,
                width = 0.01) +
  geom_errorbar(data = FIMO_param5_benchraptor_shufflestats_df[FIMO_param5_benchraptor_shufflestats_df$method == "mlm_estimate", ], aes(x = mean_auroc, y = mean_auprc,
                                                                                                                                        xmin = (mean_auroc - sd_auroc),
                                                                                                                                        xmax = (mean_auroc + sd_auroc)), alpha = 0.5,
                width = 0.01)

dev.off()

#### eY1H plots ####

eY1H_benchraptor <- read.table("output/benchmark_out/e1YH_benchRAPToR.tsv",
                               sep= '\t',
                               header = TRUE)

eY1H_shufflestats <- read.table("output/benchmark_out/e1YH_GRN_benchRAPToR_shufflestats.tsv",
                                sep= '\t',
                                header = TRUE)

eY1H_auprc_stats <- eY1H_shufflestats[eY1H_shufflestats$metric == "auprc", ] %>% group_by(method) %>% summarise(mean_auprc = mean(score),
                                                                                                                sd_auprc = sd(score))

eY1H_auprc_stats[, "auprc_95CI_margin"] <- qt(0.975, df = 99) * eY1H_auprc_stats$sd_auprc / sqrt(99)

eY1H_auroc_stats <- eY1H_shufflestats[eY1H_shufflestats$metric == "auroc", ] %>% group_by(method) %>% summarise(mean_auroc = mean(score),
                                                                                                                sd_auroc = sd(score))

eY1H_auroc_stats[, "auroc_95CI_margin"] <- qt(0.975, df = 99) * eY1H_auroc_stats$sd_auroc / sqrt(99)

eY1H_shufflestats_df <- merge(eY1H_auprc_stats,
                              eY1H_auroc_stats)

eY1H_benchraptormethods_plot <- eY1H_benchraptor[eY1H_benchraptor$metric %in% c("auroc", "auprc"), c("net", "method", "metric", "score")]
eY1H_benchraptormethods_plot <- tidyr::pivot_wider(eY1H_benchraptormethods_plot, names_from = c("metric"), values_from = "score")

eY1H_benchraptormethods_plot[, "label"] <- eY1H_benchraptormethods_plot$net
eY1H_benchraptormethods_plot[str_detect(unlist(eY1H_benchraptormethods_plot$net), "shuffle"), "label"] <- ""

method_labels <- c("Consensus", "MLM", "ULM", "WSum")
names(method_labels) <- c("consensus_estimate", "mlm_estimate", "ulm_estimate", "wsum_norm")

saveRDS(eY1H_benchraptormethods_plot,
        "plotdata/eY1H_benchraptormethods_plot.rds")
saveRDS(eY1H_shufflestats_df,
        "plotdata/eY1H_shufflestats_df.rds")

# eY1H_benchraptormethods_plot <- readRDS("plotdata/eY1H_benchraptormethods_plot.rds")
# eY1H_shufflestats_df <- readRDS("plotdata/eY1H_benchraptor_shufflestats_df.rds")

pdf("graphics/eY1H_benchraptor_methodsplot.pdf",
    height = 4,
    width = 5)

ggplot(data = eY1H_benchraptormethods_plot[!eY1H_benchraptormethods_plot$method %in% c("wsum_estimate", "wsum_corr"), ], aes(x = auroc, y = auprc)) + 
  geom_point(size = 2) +
  # geom_label_repel(aes(label = method),
  #                  label.size = 0.1,
  #                  max.overlaps = 15,
  #                  box.padding = 0.3,
  #                  label.padding = 0.1,
  #                  point.padding = 0,
  #                  size = 2) +
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
        panel.grid.minor = element_line(colour = "grey", linewidth = 0.1)) +
  geom_point(data = eY1H_shufflestats_df[!eY1H_shufflestats_df$method %in% c("wsum_estimate", "wsum_corr"), ], aes(x = mean_auroc, y = mean_auprc), alpha = 0.5) + 
  geom_errorbar(data = eY1H_shufflestats_df[!eY1H_shufflestats_df$method %in% c("wsum_estimate", "wsum_corr"), ], aes(x = mean_auroc, y = mean_auprc,
                                                                                                                      ymin = (mean_auprc - sd_auprc),
                                                                                                                      ymax = (mean_auprc + sd_auprc)), alpha = 0.5, width = 0.01) +
  geom_errorbar(data = eY1H_shufflestats_df[!eY1H_shufflestats_df$method %in% c("wsum_estimate", "wsum_corr"), ], aes(x = mean_auroc, y = mean_auprc,
                                                                                                                      xmin = (mean_auroc - sd_auroc),
                                                                                                                      xmax = (mean_auroc + sd_auroc)), alpha = 0.5, width = 0.01)


dev.off()

pdf("graphics/eY1H_benchraptor_mlmplot.pdf",
    height = 2.5,
    width = 2.5)

ggplot(data = eY1H_benchraptormethods_plot[eY1H_benchraptormethods_plot$method == "mlm_estimate", ], aes(x = auroc, y = auprc)) + 
  geom_point(size = 2) +
  # geom_label_repel(aes(label = label),
  #                  label.size = 0.1,
  #                  max.overlaps = 10,
  #                  box.padding = 0.3,
  #                  label.padding = 0.1,
  #                  point.padding = 0,
  #                  size = 2) +
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
  geom_point(data = eY1H_shufflestats_df[eY1H_shufflestats_df$method == "mlm_estimate", ], aes(x = mean_auroc, y = mean_auprc), alpha = 0.5) + 
  geom_errorbar(data = eY1H_shufflestats_df[eY1H_shufflestats_df$method == "mlm_estimate", ], aes(x = mean_auroc, y = mean_auprc,
                                                                                                  ymin = (mean_auprc - sd_auprc),
                                                                                                  ymax = (mean_auprc + sd_auprc)), alpha = 0.5, width = 0.01) +
  geom_errorbar(data = eY1H_shufflestats_df[eY1H_shufflestats_df$method == "mlm_estimate", ], aes(x = mean_auroc, y = mean_auprc,
                                                                                                  xmin = (mean_auroc - sd_auroc),
                                                                                                  xmax = (mean_auroc + sd_auroc)), alpha = 0.5, width = 0.01)
dev.off()

#### SJARACNE plots ####

SJARACNE_benchraptor <- read.table("output/benchmark_out/SJARACNE_benchRAPToR.tsv",
                                   sep= '\t',
                                   header = TRUE)

SJARACNE_shufflestats <- read.table("output/benchmark_out/SJARACNE_GRN_benchRAPToR_shufflestats.tsv",
                                    sep= '\t',
                                    header = TRUE)

SJARACNE_auprc_stats <- SJARACNE_shufflestats[SJARACNE_shufflestats$metric == "auprc", ] %>% group_by(method) %>% summarise(mean_auprc = mean(score),
                                                                                                                            sd_auprc = sd(score))

SJARACNE_auprc_stats[, "auprc_95CI_margin"] <- qt(0.975, df = 99) * SJARACNE_auprc_stats$sd_auprc / sqrt(99)

SJARACNE_auroc_stats <- SJARACNE_shufflestats[SJARACNE_shufflestats$metric == "auroc", ] %>% group_by(method) %>% summarise(mean_auroc = mean(score),
                                                                                                                            sd_auroc = sd(score))

SJARACNE_auroc_stats[, "auroc_95CI_margin"] <- qt(0.975, df = 99) * SJARACNE_auroc_stats$sd_auroc / sqrt(99)

SJARACNE_shufflestats_df <- merge(SJARACNE_auprc_stats,
                                  SJARACNE_auroc_stats)

SJARACNE_benchraptormethods_plot <- SJARACNE_benchraptor[SJARACNE_benchraptor$metric %in% c("auroc", "auprc"), c("net", "method", "metric", "score")]
SJARACNE_benchraptormethods_plot <- tidyr::pivot_wider(SJARACNE_benchraptormethods_plot, names_from = c("metric"), values_from = "score")

SJARACNE_benchraptormethods_plot[, "label"] <- SJARACNE_benchraptormethods_plot$net
SJARACNE_benchraptormethods_plot[str_detect(unlist(SJARACNE_benchraptormethods_plot$net), "shuffle"), "label"] <- ""

method_labels <- c("Consensus", "MLM", "ULM", "WSum")
names(method_labels) <- c("consensus_estimate", "mlm_estimate", "ulm_estimate", "wsum_norm")

saveRDS(SJARACNE_benchraptormethods_plot,
        "plotdata/SJARACNE_benchraptormethods_plot.rds")
# SJARACNE_benchraptormethods_plot <- readRDS("plotdata/SJARACNE_benchraptormethods_plot.rds")

saveRDS(SJARACNE_shufflestats_df,
        "plotdata/SJARACNE_benchraptor_shufflestats_df.rds")
SJARACNE_shufflestats_df <- readRDS("plotdata/SJARACNE_benchraptor_shufflestats_df.rds")

pdf("graphics/SJARACNE_benchraptor_methodsplot.pdf",
    height = 4,
    width = 5)

ggplot(data = SJARACNE_benchraptormethods_plot[!SJARACNE_benchraptormethods_plot$method %in% c("wsum_estimate", "wsum_corr"), ], aes(x = auroc, y = auprc)) + 
  geom_point(size = 2) +
  # geom_label_repel(aes(label = method),
  #                  label.size = 0.1,
  #                  max.overlaps = 15,
  #                  box.padding = 0.3,
  #                  label.padding = 0.1,
  #                  point.padding = 0,
  #                  size = 2) +
  geom_vline(xintercept = 0.5,
             linetype = "dashed",
             col = "black") +
  geom_hline(yintercept = 0.5,
             linetype = "dashed",
             col = "black") + 
  facet_wrap(~ method, labeller = labeller(method = method_labels)) + 
  xlab("Sensitivity (AUROC)") + 
  ylab("Precision (AUPRC)") +
  coord_cartesian(xlim = c(0.40, 0.75),
                  ylim = c(0.40, 0.8)) +
  theme(strip.text = element_text(colour = 'white'),
        strip.background = element_rect(fill = "black"),
        legend.position = "none",
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_line(colour = "grey", linewidth = 0.2),
        panel.grid.minor = element_line(colour = "grey", linewidth = 0.1),
        axis.text = element_text(colour = "black")) +
  geom_point(data = SJARACNE_shufflestats_df[!SJARACNE_shufflestats_df$method %in% c("wsum_estimate", "wsum_corr"), ], aes(x = mean_auroc, y = mean_auprc), alpha = 0.5) + 
  geom_errorbar(data = SJARACNE_shufflestats_df[!SJARACNE_shufflestats_df$method %in% c("wsum_estimate", "wsum_corr"), ], aes(x = mean_auroc, y = mean_auprc,
                                                                                                                              ymin = (mean_auprc - sd_auprc),
                                                                                                                              ymax = (mean_auprc + sd_auprc)), alpha = 0.5, width = 0.01) +
  geom_errorbar(data = SJARACNE_shufflestats_df[!SJARACNE_shufflestats_df$method %in% c("wsum_estimate", "wsum_corr"), ], aes(x = mean_auroc, y = mean_auprc,
                                                                                                                              xmin = (mean_auroc - sd_auroc),
                                                                                                                              xmax = (mean_auroc + sd_auroc)), alpha = 0.5, width = 0.01)


dev.off()

pdf("graphics/SJARACNE_benchraptor_mlmplot.pdf",
    height = 2.5,
    width = 2.5)

ggplot(data = SJARACNE_benchraptormethods_plot[SJARACNE_benchraptormethods_plot$method == "mlm_estimate", ], aes(x = auroc, y = auprc)) + 
  geom_point(size = 2) +
  # geom_label_repel(aes(label = label),
  #                  label.size = 0.1,
  #                  max.overlaps = 10,
  #                  box.padding = 0.3,
  #                  label.padding = 0.1,
  #                  point.padding = 0,
  #                  size = 2) +
  geom_vline(xintercept = 0.5,
             linetype = "dashed",
             col = "black") +
  geom_hline(yintercept = 0.5,
             linetype = "dashed",
             col = "black") + 
  xlab("Sensitivity (AUROC)") + 
  ylab("Precision (AUPRC)") +
  coord_cartesian(xlim = c(0.40, 0.75),
                  ylim = c(0.40, 0.8)) +
  theme(strip.text = element_text(colour = 'white'),
        strip.background = element_rect(fill = "black"),
        legend.position = "none",
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_line(colour = "grey", linewidth = 0.2),
        panel.grid.minor = element_line(colour = "grey", linewidth = 0.1),
        axis.text = element_text(colour = "black")) +
  geom_point(data = SJARACNE_shufflestats_df[SJARACNE_shufflestats_df$method == "mlm_estimate", ], aes(x = mean_auroc, y = mean_auprc), alpha = 0.5) + 
  geom_errorbar(data = SJARACNE_shufflestats_df[SJARACNE_shufflestats_df$method == "mlm_estimate", ], aes(x = mean_auroc, y = mean_auprc,
                                                                                                          ymin = (mean_auprc - sd_auprc),
                                                                                                          ymax = (mean_auprc + sd_auprc)), alpha = 0.5, width = 0.01) +
  geom_errorbar(data = SJARACNE_shufflestats_df[SJARACNE_shufflestats_df$method == "mlm_estimate", ], aes(x = mean_auroc, y = mean_auprc,
                                                                                                          xmin = (mean_auroc - sd_auroc),
                                                                                                          xmax = (mean_auroc + sd_auroc)), alpha = 0.5, width = 0.01)
dev.off()

#### NETWORK COMBINATION PLOTS ####

#### COMBO SHARED SET ####

combo_restrictedobs_benchraptor_files <- bench_out_filelist[str_detect(bench_out_filelist, "benchRAPToR") &
                                                              str_detect(bench_out_filelist, "restrictedobs") &
                                                              !str_detect(bench_out_filelist, "shufflestats")]

combo_restrictedobs_benchraptor_shufflestat_files <- bench_out_filelist[str_detect(bench_out_filelist, "benchRAPToR") &
                                                                          str_detect(bench_out_filelist, "restrictedobs") &
                                                                          str_detect(bench_out_filelist, "shufflestats")]

combo_benchraptor <- read.table(paste0("output/benchmark_out/", combo_restrictedobs_benchraptor_files),
                                sep= '\t',
                                header = TRUE)

combo_benchraptormethods_plot <- combo_benchraptor[combo_benchraptor$metric %in% c("auroc", "auprc"), c("net", "method", "metric", "score")]
combo_benchraptormethods_plot <- tidyr::pivot_wider(combo_benchraptormethods_plot, names_from = c("metric"), values_from = "score")

combo_benchraptormethods_plot[, "label"] <- combo_benchraptormethods_plot$net
combo_benchraptormethods_plot[str_detect(unlist(combo_benchraptormethods_plot$net), "shuffle"), "label"] <- ""

method_labels <- c("Consensus", "MLM", "ULM", "WSum")
names(method_labels) <- c("consensus_estimate", "mlm_estimate", "ulm_estimate", "wsum_norm")

combo_restrictedobs_labelvec <- c("motif", "ChIP", "not weighted", "weighted", "uneq. weights")
names(combo_restrictedobs_labelvec) <- unique(combo_benchraptormethods_plot$net)

combo_benchraptormethods_plot[, "label"] <- combo_restrictedobs_labelvec[combo_benchraptormethods_plot$net]

saveRDS(combo_benchraptormethods_plot,
        "plotdata/combo_benchraptormethods_plot.rds")
# combo_benchraptormethods_plot <- readRDS("plotdata/combo_benchraptormethods_plot.rds")

combo_restrictedobs_benchraptor_shufflestats_lists <- lapply(combo_restrictedobs_benchraptor_shufflestat_files, function(file){
  
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
  
  tempout["net"] <- str_remove(file, "restrictedobs_benchRAPToR_shufflestats.tsv")
  
  tempout
  
})

combo_restrictedobs_benchraptor_shufflestats_df <- do.call(rbind, combo_restrictedobs_benchraptor_shufflestats_lists)

combo_restrictedobs_benchraptor_shufflestats_df[, "label"] <- combo_restrictedobs_labelvec[combo_restrictedobs_benchraptor_shufflestats_df$net]

saveRDS(combo_restrictedobs_benchraptor_shufflestats_df,
        "plotdata/combo_restrictedobs_benchraptor_shufflestats_df.rds")
# combo_restrictedobs_benchraptor_shufflestats_df <- readRDS("plotdata/combo_restrictedobs_benchraptor_shufflestats_df.rds")

combo_benchraptormethods_plot$net <- factor(combo_benchraptormethods_plot$net, 
                                            levels = c("ChIPnocut", "FIMO1000", "FIMO1000ChIPnocut_noweights", "FIMO1000ChIPnocut_weighted"))

pdf("graphics/combo_restrictedobs_benchRAPToR_methodsplot.pdf",
    height = 4,
    width = 5)

ggplot(data = combo_benchraptormethods_plot[!combo_benchraptormethods_plot$method %in% c("wsum_estimate", "wsum_corr"), ], aes(x = auroc, y = auprc, colour = net)) + 
  geom_point(size = 2) +
  geom_label_repel(aes(label = label,
                       fill = net),
                   label.size = 0.1,
                   max.overlaps = 15,
                   box.padding = 0.5,
                   label.padding = 0.1,
                   point.padding = 0.5,
                   size = 2,
                   colour = "black") +
  scale_colour_manual(values = c("orange", "magenta", "dodger blue", "green")) +
  scale_fill_manual(values =c("orange", "magenta", "dodger blue", "green")) +
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
  geom_point(data = combo_restrictedobs_benchraptor_shufflestats_df[!combo_restrictedobs_benchraptor_shufflestats_df$method %in% c("wsum_estimate", "wsum_corr"), ], aes(x = mean_auroc, y = mean_auprc), alpha = 0.5) + 
  geom_errorbar(data = combo_restrictedobs_benchraptor_shufflestats_df[!combo_restrictedobs_benchraptor_shufflestats_df$method %in% c("wsum_estimate", "wsum_corr"), ], aes(x = mean_auroc, y = mean_auprc,
                                                                                                                                                                            ymin = (mean_auprc - sd_auprc),
                                                                                                                                                                            ymax = (mean_auprc + sd_auprc)), alpha = 0.5,
                width = 0.01) +
  geom_errorbar(data = combo_restrictedobs_benchraptor_shufflestats_df[!combo_restrictedobs_benchraptor_shufflestats_df$method %in% c("wsum_estimate", "wsum_corr"), ], aes(x = mean_auroc, y = mean_auprc,
                                                                                                                                                                            xmin = (mean_auroc - sd_auroc),
                                                                                                                                                                            xmax = (mean_auroc + sd_auroc)), alpha = 0.5, 
                width = 0.01)

dev.off()

pdf("graphics/combo_restrictedobs_benchraptor_mlmplot.pdf",
    height = 2.5,
    width = 2.5)

ggplot(data = combo_benchraptormethods_plot[!str_detect(combo_benchraptormethods_plot$net, "shuffle") & 
                                              combo_benchraptormethods_plot$method == "mlm_estimate", ], aes(x = auroc, y = auprc, colour = net)) + 
  geom_point(size = 2) +
  geom_label_repel(aes(label = label,
                       fill = net),
                   colour = "black", 
                   label.size = 0.1,
                   max.overlaps = 10,
                   box.padding = 0,
                   label.padding = 0.1,
                   point.padding = 2,
                   size = 1.5) +
  scale_colour_manual(values = c("orange", "magenta", "dodger blue", "green")) +
  scale_fill_manual(values =c("orange", "magenta", "dodger blue", "green")) +
  geom_vline(xintercept = 0.5,
             linetype = "dashed",
             col = "black") +
  geom_hline(yintercept = 0.5,
             linetype = "dashed",
             col = "black") + 
  xlab("Sensitivity (AUROC)") + 
  ylab("Precision (AUPRC)") +
  coord_cartesian(xlim = c(0.45, 0.75),
                  ylim = c(0.45, 0.85)) +
  theme(strip.text = element_text(colour = 'white'),
        strip.background = element_rect(fill = "black"),
        legend.position = "none",
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_line(colour = "grey", linewidth = 0.2),
        panel.grid.minor = element_line(colour = "grey", linewidth = 0.1),
        axis.text = element_text(colour = "black")) +
  geom_point(data = combo_restrictedobs_benchraptor_shufflestats_df[combo_restrictedobs_benchraptor_shufflestats_df$method == "mlm_estimate", ], aes(x = mean_auroc, y = mean_auprc), alpha = 0.5) + 
  geom_errorbar(data = combo_restrictedobs_benchraptor_shufflestats_df[combo_restrictedobs_benchraptor_shufflestats_df$method == "mlm_estimate", ], aes(x = mean_auroc, y = mean_auprc,
                                                                                                                                                        ymin = (mean_auprc - sd_auprc),
                                                                                                                                                        ymax = (mean_auprc + sd_auprc)), alpha = 0.5,
                width = 0.01) +
  geom_errorbar(data = combo_restrictedobs_benchraptor_shufflestats_df[combo_restrictedobs_benchraptor_shufflestats_df$method == "mlm_estimate", ], aes(x = mean_auroc, y = mean_auprc,
                                                                                                                                                        xmin = (mean_auroc - sd_auroc),
                                                                                                                                                        xmax = (mean_auroc + sd_auroc)), alpha = 0.5,
                width = 0.01)
dev.off()

#### NOW FOR COMBO FULL SET ####

combo3_benchraptor_files <- bench_out_filelist[str_detect(bench_out_filelist, "benchRAPToR") &
                                                 str_detect(bench_out_filelist, "combo3") &
                                                 !str_detect(bench_out_filelist, "shufflestats")]

combo3_benchraptor_shufflestat_files <- bench_out_filelist[str_detect(bench_out_filelist, "benchRAPToR") &
                                                             str_detect(bench_out_filelist, "combo3") &
                                                             str_detect(bench_out_filelist, "shufflestats")]

combo3_benchraptor <- read.table(paste0("output/benchmark_out/", combo3_benchraptor_files),
                                 sep= '\t',
                                 header = TRUE)

combo3_benchraptormethods_plot <- combo3_benchraptor[combo3_benchraptor$metric %in% c("auroc", "auprc"), c("net", "method", "metric", "score")]
combo3_benchraptormethods_plot <- tidyr::pivot_wider(combo3_benchraptormethods_plot, names_from = c("metric"), values_from = "score")

combo3_benchraptormethods_plot[, "label"] <- combo3_benchraptormethods_plot$net
combo3_benchraptormethods_plot[str_detect(unlist(combo3_benchraptormethods_plot$net), "shuffle"), "label"] <- ""

combo3_labelvec <- c("weighted", "not weighted")
names(combo3_labelvec) <- unique(combo3_benchraptormethods_plot$net)[1:2]

combo3_benchraptormethods_plot[, "label"] <- combo3_labelvec[combo3_benchraptormethods_plot$net]

combo3_benchraptor_shufflestats_lists <- lapply(combo3_benchraptor_shufflestat_files, function(file){
  
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

combo3_benchraptor_shufflestats_df <- do.call(rbind, combo3_benchraptor_shufflestats_lists)

combo3_benchraptor_shufflestats_df[, "label"] <- combo3_labelvec[combo3_benchraptor_shufflestats_df$net]

saveRDS(combo3_benchraptormethods_plot,
        "plotdata/combo3_benchraptormethods_plot.rds")
# combo3_benchraptormethods_plot <- readRDS("plotdata/combo3_benchraptormethods_plot.rds")
saveRDS(combo3_benchraptor_shufflestats_df,
        "plotdata/combo3_benchraptor_shufflestats_df.rds")
# combo3_benchraptor_shufflestats_df <- readRDS(
#         "plotdata/combo3_benchraptor_shufflestats_df.rds")

pdf("graphics/combo3_benchRAPToR_methodsplot.pdf",
    height = 4,
    width = 5)

ggplot(data = combo3_benchraptormethods_plot[!combo3_benchraptormethods_plot$method %in% c("wsum_estimate", "wsum_corr") &
                                               combo3_benchraptormethods_plot$net != "combo3_unequal", ], aes(x = auroc, y = auprc, colour = net)) + 
  geom_point(size = 2) +
  geom_label_repel(aes(label = label,
                       fill = net),
                   colour = "black",
                   label.size = 0.1,
                   max.overlaps = 15,
                   box.padding = 0.3,
                   label.padding = 0.1,
                   point.padding = 0,
                   size = 2) +
  scale_colour_manual(values = c("green", "dodger blue")) +
  scale_fill_manual(values = c("green", "dodger blue")) +
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
  geom_point(data = combo3_benchraptor_shufflestats_df[!combo3_benchraptor_shufflestats_df$method %in% c("wsum_estimate", "wsum_corr"), ], aes(x = mean_auroc, y = mean_auprc), alpha = 0.5) + 
  geom_errorbar(data = combo3_benchraptor_shufflestats_df[!combo3_benchraptor_shufflestats_df$method %in% c("wsum_estimate", "wsum_corr"), ], aes(x = mean_auroc, y = mean_auprc,
                                                                                                                                                  ymin = (mean_auprc - sd_auprc),
                                                                                                                                                  ymax = (mean_auprc + sd_auprc)), alpha = 0.5, width = 0.01) +
  geom_errorbar(data = combo3_benchraptor_shufflestats_df[!combo3_benchraptor_shufflestats_df$method %in% c("wsum_estimate", "wsum_corr"), ], aes(x = mean_auroc, y = mean_auprc,
                                                                                                                                                  xmin = (mean_auroc - sd_auroc),
                                                                                                                                                  xmax = (mean_auroc + sd_auroc)), alpha = 0.5, width = 0.01)

dev.off()

pdf("graphics/combo3_benchraptor_mlmplot.pdf",
    height = 2.5,
    width = 2.5)

ggplot(data = combo3_benchraptormethods_plot[!str_detect(combo3_benchraptormethods_plot$net, "shuffle") & 
                                               combo3_benchraptormethods_plot$method == "mlm_estimate" &
                                               combo3_benchraptormethods_plot$net != "combo3_unequal", ], aes(x = auroc, y = auprc, colour = net)) + 
  geom_point(size = 2) +
  geom_label_repel(aes(label = label,
                       fill = net),
                   colour = "black",
                   label.size = 0.1,
                   max.overlaps = 10,
                   box.padding = 0.2,
                   label.padding = 0.1,
                   point.padding = 0,
                   size = 3) +
  scale_colour_manual(values = c("green", "dodger blue")) +
  scale_fill_manual(values = c("green", "dodger blue")) +
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
  geom_point(data = combo3_benchraptor_shufflestats_df[combo3_benchraptor_shufflestats_df$method == "mlm_estimate", ], aes(x = mean_auroc, y = mean_auprc), alpha = 0.5) + 
  geom_errorbar(data = combo3_benchraptor_shufflestats_df[combo3_benchraptor_shufflestats_df$method == "mlm_estimate", ], aes(x = mean_auroc, y = mean_auprc,
                                                                                                                              ymin = (mean_auprc - sd_auprc),
                                                                                                                              ymax = (mean_auprc + sd_auprc)), alpha = 0.5,
                width = 0.01) +
  geom_errorbar(data = combo3_benchraptor_shufflestats_df[combo3_benchraptor_shufflestats_df$method == "mlm_estimate", ], aes(x = mean_auroc, y = mean_auprc,
                                                                                                                              xmin = (mean_auroc - sd_auroc),
                                                                                                                              xmax = (mean_auroc + sd_auroc)), alpha = 0.5,
                width = 0.01)
dev.off()

#### COMBO FULL SET RAW ####

combo3_benchraw_files <- bench_out_filelist[str_detect(bench_out_filelist, "benchraw") &
                                                 str_detect(bench_out_filelist, "combo3") &
                                                 !str_detect(bench_out_filelist, "shufflestats")]

combo3_benchraw_shufflestat_files <- bench_out_filelist[str_detect(bench_out_filelist, "benchraw") &
                                                             str_detect(bench_out_filelist, "combo3") &
                                                             str_detect(bench_out_filelist, "shufflestats")]

combo3_benchraw <- read.table(paste0("output/benchmark_out/", combo3_benchraw_files),
                                 sep= '\t',
                                 header = TRUE)

combo3_benchrawmethods_plot <- combo3_benchraw[combo3_benchraw$metric %in% c("auroc", "auprc"), c("net", "method", "metric", "score")]
combo3_benchrawmethods_plot <- tidyr::pivot_wider(combo3_benchrawmethods_plot, names_from = c("metric"), values_from = "score")

combo3_benchrawmethods_plot[, "label"] <- combo3_benchrawmethods_plot$net
combo3_benchrawmethods_plot[str_detect(unlist(combo3_benchrawmethods_plot$net), "shuffle"), "label"] <- ""

combo3_labelvec <- c("weighted", "not weighted")
names(combo3_labelvec) <- unique(combo3_benchrawmethods_plot$net)[1:2]

combo3_benchrawmethods_plot[, "label"] <- combo3_labelvec[combo3_benchrawmethods_plot$net]

combo3_benchraw_shufflestats_lists <- lapply(combo3_benchraw_shufflestat_files, function(file){
  
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
  
  tempout["net"] <- str_remove(file, "_benchraw_shufflestats.tsv")
  
  tempout
  
})

combo3_benchraw_shufflestats_df <- do.call(rbind, combo3_benchraw_shufflestats_lists)

combo3_benchraw_shufflestats_df[, "label"] <- combo3_labelvec[combo3_benchraw_shufflestats_df$net]

saveRDS(combo3_benchrawmethods_plot,
        "plotdata/combo3_benchrawmethods_plot.rds")
# combo3_benchrawmethods_plot <- readRDS("plotdata/combo3_benchrawmethods_plot.rds")

saveRDS(combo3_benchraw_shufflestats_df,
        "plotdata/combo3_benchraw_shufflestats_df.rds")
combo3_benchraw_shufflestats_df <- readRDS("plotdata/combo3_benchraw_shufflestats_df.rds")

pdf("graphics/combo3_benchraw_methodsplot.pdf",
    height = 4,
    width = 5)

ggplot(data = combo3_benchrawmethods_plot[!combo3_benchrawmethods_plot$method %in% c("wsum_estimate", "wsum_corr") &
                                               combo3_benchrawmethods_plot$net != "combo3_unequal", ], aes(x = auroc, y = auprc, colour = net)) + 
  geom_point(size = 2) +
  geom_label_repel(aes(label = label,
                       fill = net),
                   colour = "black",
                   label.size = 0.1,
                   max.overlaps = 15,
                   box.padding = 0.3,
                   label.padding = 0.1,
                   point.padding = 0,
                   size = 2) +
  geom_vline(xintercept = 0.5,
             linetype = "dashed",
             col = "black") +
  scale_colour_manual(values = c("green", "dodger blue")) +
  scale_fill_manual(values = c("green", "dodger blue")) +
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
  geom_point(data = combo3_benchraw_shufflestats_df[!combo3_benchraw_shufflestats_df$method %in% c("wsum_estimate", "wsum_corr"), ], aes(x = mean_auroc, y = mean_auprc), alpha = 0.5) + 
  geom_errorbar(data = combo3_benchraw_shufflestats_df[!combo3_benchraw_shufflestats_df$method %in% c("wsum_estimate", "wsum_corr"), ], aes(x = mean_auroc, y = mean_auprc,
                                                                                                                                                  ymin = (mean_auprc - sd_auprc),
                                                                                                                                                  ymax = (mean_auprc + sd_auprc)), alpha = 0.5, width = 0.01) +
  geom_errorbar(data = combo3_benchraw_shufflestats_df[!combo3_benchraw_shufflestats_df$method %in% c("wsum_estimate", "wsum_corr"), ], aes(x = mean_auroc, y = mean_auprc,
                                                                                                                                                  xmin = (mean_auroc - sd_auroc),
                                                                                                                                                  xmax = (mean_auroc + sd_auroc)), alpha = 0.5, width = 0.01)

dev.off()

pdf("graphics/combo3_benchraw_mlmplot.pdf",
    height = 3.5,
    width = 3.5)

ggplot(data = combo3_benchrawmethods_plot[!str_detect(combo3_benchrawmethods_plot$net, "shuffle") & 
                                               combo3_benchrawmethods_plot$method == "mlm_estimate" &
                                               combo3_benchrawmethods_plot$net != "combo3_unequal", ], aes(x = auroc, y = auprc, colour = net)) + 
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
  scale_colour_manual(values = c("green", "dodger blue")) + 
  coord_cartesian(xlim = c(0.45, 0.75),
                  ylim = c(0.45, 0.8)) +
  theme(strip.text = element_text(colour = 'white'),
        strip.background = element_rect(fill = "black"),
        legend.position = "none",
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_line(colour = "grey", linewidth = 0.2),
        panel.grid.minor = element_line(colour = "grey", linewidth = 0.1),
        axis.text = element_text(colour = "black")) +
  geom_point(data = combo3_benchraw_shufflestats_df[combo3_benchraw_shufflestats_df$method == "mlm_estimate", ], aes(x = mean_auroc, y = mean_auprc), alpha = 0.5) + 
  geom_errorbar(data = combo3_benchraw_shufflestats_df[combo3_benchraw_shufflestats_df$method == "mlm_estimate", ], aes(x = mean_auroc, y = mean_auprc,
                                                                                                                              ymin = (mean_auprc - sd_auprc),
                                                                                                                              ymax = (mean_auprc + sd_auprc)), alpha = 0.5,
                width = 0.01) +
  geom_errorbar(data = combo3_benchraw_shufflestats_df[combo3_benchraw_shufflestats_df$method == "mlm_estimate", ], aes(x = mean_auroc, y = mean_auprc,
                                                                                                                              xmin = (mean_auroc - sd_auroc),
                                                                                                                              xmax = (mean_auroc + sd_auroc)), alpha = 0.5,
                width = 0.01)
dev.off()

