#start with empty workspace

rm(list = ls(all = TRUE))

# turn off scientific notation for plots

options(scipen=10000)

#### set working directory ####

# here create new folder and set working directory within it

dir.create("~/Cel_GRN_manuscript/")
setwd("~/Cel_GRN_manuscript/")

# create subfolders for input, output and graphics

dir.create("input")

# into input folder, add input files 

dir.create("output")

dir.create("output/GRNs")

dir.create("graphics")

#### DEFINE FUNCTIONS ####


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

method_labels <- c("Consensus", "MLM", "ULM", "WSum")
names(method_labels) <- c("consensus_estimate", "mlm_estimate", "ulm_estimate", "wsum_norm")

bench_out_filelist <- list.files("output/benchmark_out/")
shufflestats_files <- bench_out_filelist[str_detect(bench_out_filelist, "_shufflestats")]

allChIP_HOTincl_benchraw_files <- bench_out_filelist[str_detect(bench_out_filelist, "benchRAW") & str_detect(bench_out_filelist, "HOTincl")]

allChIP_HOTincl_benchraw_files <- allChIP_HOTincl_benchraw_files[!str_detect(allChIP_HOTincl_benchraw_files, "dense")]
allChIP_HOTincl_benchraw_files <- allChIP_HOTincl_benchraw_files[!str_detect(allChIP_HOTincl_benchraw_files, "shufflestat")]

allChIP_HOTincl_benchraw_list <- lapply(allChIP_HOTincl_benchraw_files, function(x){
  
  read.table(paste0("output/benchmark_out/", x),
             sep = "\t",
             header = TRUE)
  
})

allChIP_HOTincl_benchraw_all <- do.call(rbind, allChIP_HOTincl_benchraw_list)

allChIP_HOTincl_benchraw_allmethods_plot <- allChIP_HOTincl_benchraw_all[allChIP_HOTincl_benchraw_all$metric %in% c("auroc", "auprc"), c("net", "method", "metric", "score")]
allChIP_HOTincl_benchraw_allmethods_plot <- tidyr::pivot_wider(allChIP_HOTincl_benchraw_allmethods_plot, names_from = c("metric"), values_from = "score")

allChIP_HOTincl_benchraw_allmethods_plot[, "label"] <- allChIP_HOTincl_benchraw_allmethods_plot$net
# allChIP_HOTincl_benchraw_allmethods_plot[str_detect(unlist(allChIP_HOTincl_benchraw_allmethods_plot$net), "shuffle"), "label"] <- ""

allChIP_HOTincl_benchraw_allmethods_plot[, "cutoff"] <- str_remove(allChIP_HOTincl_benchraw_allmethods_plot$net, "_shuffle")
# allChIP_HOTincl_benchraw_allmethods_plot[, "shuffled"] <- str_detect(allChIP_HOTincl_benchraw_allmethods_plot$net, "_shuffle")

benchRAW_HOTincl_shufflestats_files <- shufflestats_files[str_detect(shufflestats_files, "RAW") & str_detect(shufflestats_files, "HOTincl")]

benchRAW_HOTincl_shufflestats_lists <- lapply(benchRAW_HOTincl_shufflestats_files, function(file){
  
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

benchRAW_HOTincl_shufflestats_df <- do.call(rbind, benchRAW_HOTincl_shufflestats_lists)
benchRAW_HOTincl_shufflestats_df[, "cutoff"] <- str_remove(benchRAW_HOTincl_shufflestats_df$net, "_shuffle")

saveRDS(allChIP_HOTincl_benchraw_allmethods_plot,
        "plotdata/allChIP_HOTincl_benchraw_allmethods_plot.rds")

saveRDS(benchRAW_HOTincl_shufflestats_df,
        "plotdata/benchRAW_HOTincl_shufflestats_df.rds")

# pdf("graphics/allChIP_HOTincl_benchRAW_methodsplot.pdf",
#     height = 4,
#     width = 5)
# 
# ggplot(data = allChIP_HOTincl_benchraw_allmethods_plot[!str_detect(allChIP_HOTincl_benchrw_allmethods_plot$net, "shuffle") & 
#                                                          !allChIP_HOTincl_benchraw_allmethods_plot$cutoff %in% c("1500", "3000", "10000") &
#                                                          !allChIP_HOTincl_benchraw_allmethods_plot$method %in% c("wsum_estimate", "wsum_corr"), ], aes(x = auroc, y = auprc, colour = cutoff)) + 
#   geom_point(size = 2) +
#   geom_label_repel(aes(label = label),
#                    label.size = 0.1,
#                    max.overlaps = 15,
#                    box.padding = 0.3,
#                    label.padding = 0.1,
#                    point.padding = 0,
#                    size = 2) +
#   geom_vline(xintercept = 0.5,
#              linetype = "dashed",
#              col = "black") +
#   geom_hline(yintercept = 0.5,
#              linetype = "dashed",
#              col = "black") + 
#   facet_wrap(~ method, labeller = labeller(method = method_labels)) + 
#   xlab("Sensitivity (AUROC)") + 
#   ylab("Precision (AUPRC)") +
#   coord_cartesian(xlim = c(0.45, 0.75),
#                   ylim = c(0.45, 0.8)) +
#   theme(strip.text = element_text(colour = 'white'),
#         strip.background = element_rect(fill = "black"),
#         legend.position = "none",
#         panel.background = element_rect(fill = "white", colour = "black"),
#         panel.grid.major = element_line(colour = "grey", linewidth = 0.2),
#         panel.grid.minor = element_line(colour = "grey", linewidth = 0.1)) +
#   geom_point(data = benchRAW_HOTincl_shufflestats_df[!benchRAW_HOTincl_shufflestats_df$cutoff %in% c("1500", "3000", "10000") &
#                                                           !benchRAW_HOTincl_shufflestats_df$method %in% c("wsum_estimate", "wsum_corr"), ], aes(x = mean_auroc, y = mean_auprc), alpha = 0.5) + 
#   geom_errorbar(data = benchRAW_HOTincl_shufflestats_df[!benchRAW_HOTincl_shufflestats_df$cutoff %in% c("1500", "3000", "10000") &
#                                                              !benchRAW_HOTincl_shufflestats_df$method %in% c("wsum_estimate", "wsum_corr"), ], aes(x = mean_auroc, y = mean_auprc,
#                                                                 ymin = (mean_auprc - sd_auprc),
#                                                                 ymax = (mean_auprc + sd_auprc)), alpha = 0.5) +
#   geom_errorbar(data = benchRAW_HOTincl_shufflestats_df[!benchRAW_HOTincl_shufflestats_df$cutoff %in% c("1500", "3000", "10000") &
#                                                              !benchRAW_HOTincl_shufflestats_df$method %in% c("wsum_estimate", "wsum_corr"), ], aes(x = mean_auroc, y = mean_auprc,
#                                                                 xmin = (mean_auroc - sd_auroc),
#                                                                 xmax = (mean_auroc + sd_auroc)), alpha = 0.5)
# 
# dev.off()


## and with RAPToR-corrected benchmark stats

allChIP_HOTincl_benchraptor_files <- bench_out_filelist[str_detect(bench_out_filelist, "benchRAPToR") & str_detect(bench_out_filelist, "HOTincl")]

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

allChIP_HOTincl_benchraptor_allmethods_plot[, "label"] <- allChIP_HOTincl_benchraptor_allmethods_plot$net
# allChIP_HOTincl_benchraptor_allmethods_plot[str_detect(unlist(allChIP_HOTincl_benchraptor_allmethods_plot$net), "shuffle"), "label"] <- ""

allChIP_HOTincl_benchraptor_allmethods_plot[, "cutoff"] <- str_remove(allChIP_HOTincl_benchraptor_allmethods_plot$net, "_shuffle")
allChIP_HOTincl_benchraptor_allmethods_plot[, "shuffled"] <- str_detect(allChIP_HOTincl_benchraptor_allmethods_plot$net, "_shuffle")

benchraptor_HOTincl_shufflestats_files <- shufflestats_files[str_detect(shufflestats_files, "RAPToR") & str_detect(shufflestats_files, "HOTincl")]

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
  
  tempout["net"] <- paste0(str_extract(file, "[0-9]{3,4}"), "_shuffle")
  
  tempout
  
})

benchraptor_HOTincl_shufflestats_df <- do.call(rbind, benchraptor_HOTincl_shufflestats_lists)
benchraptor_HOTincl_shufflestats_df[, "cutoff"] <- str_remove(benchraptor_HOTincl_shufflestats_df$net, "_shuffle")

allChIP_HOTincl_benchraptor_allmethods_plot[allChIP_HOTincl_benchraptor_allmethods_plot$label == "5000", "label"] <- "none"

saveRDS(allChIP_HOTincl_benchraptor_allmethods_plot,
        "plotdata/allChIP_HOTincl_benchraptor_allmethods_plot.rds")
allChIP_HOTincl_benchraptor_allmethods_plot <- readRDS("plotdata/allChIP_HOTincl_benchraptor_allmethods_plot.rds")

saveRDS(benchraptor_HOTincl_shufflestats_df,
        "plotdata/benchRAPToR_HOTincl_shufflestats_df.rds")
benchraptor_HOTincl_shufflestats_df <- readRDS("plotdata/benchRAPToR_HOTincl_shufflestats_df.rds")

pdf("graphics/allChIP_HOTincl_benchRAPToR_methodsplot.pdf",
    height = 4,
    width = 5)

ggplot(data = allChIP_HOTincl_benchraptor_allmethods_plot[!str_detect(allChIP_HOTincl_benchraptor_allmethods_plot$net, "shuffle") & 
                                                         !allChIP_HOTincl_benchraptor_allmethods_plot$cutoff %in% c("1500", "3000", "10000") &
                                                         !allChIP_HOTincl_benchraptor_allmethods_plot$method %in% c("wsum_estimate", "wsum_corr"), ], aes(x = auroc, y = auprc, colour = cutoff)) + 
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
  geom_point(data = benchraptor_HOTincl_shufflestats_df[!benchraptor_HOTincl_shufflestats_df$cutoff %in% c("1500", "3000", "10000") &
                                                       !benchraptor_HOTincl_shufflestats_df$method %in% c("wsum_estimate", "wsum_corr"), ], aes(x = mean_auroc, y = mean_auprc), alpha = 0.5) + 
  geom_errorbar(data = benchraptor_HOTincl_shufflestats_df[!benchraptor_HOTincl_shufflestats_df$cutoff %in% c("1500", "3000", "10000") &
                                                          !benchraptor_HOTincl_shufflestats_df$method %in% c("wsum_estimate", "wsum_corr"), ], aes(x = mean_auroc, y = mean_auprc,
                                                                                                                                                ymin = (mean_auprc - sd_auprc),
                                                                                                                                                ymax = (mean_auprc + sd_auprc)), alpha = 0.5) +
  geom_errorbar(data = benchraptor_HOTincl_shufflestats_df[!benchraptor_HOTincl_shufflestats_df$cutoff %in% c("1500", "3000", "10000") &
                                                          !benchraptor_HOTincl_shufflestats_df$method %in% c("wsum_estimate", "wsum_corr"), ], aes(x = mean_auroc, y = mean_auprc,
                                                                                                                                                xmin = (mean_auroc - sd_auroc),
                                                                                                                                                xmax = (mean_auroc + sd_auroc)), alpha = 0.5)

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

#### MLM plot for main figure

# pdf("graphics/allChIP_HOTincl_benchRAW_mlmplot.pdf",
#     height = 2.5,
#     width = 2.5)
# 
# ggplot(data = allChIP_HOTincl_benchraw_allmethods_plot[!str_detect(allChIP_HOTincl_benchraw_allmethods_plot$net, "shuffle") & 
#                                                          !allChIP_HOTincl_benchraw_allmethods_plot$cutoff %in% c("1500", "3000", "10000") &
#                                                          allChIP_HOTincl_benchraw_allmethods_plot$method == "mlm_estimate", ], aes(x = auroc, y = auprc, colour = cutoff)) + 
#   geom_point(size = 2) +
#   geom_label_repel(aes(label = label),
#                    label.size = 0.1,
#                    max.overlaps = 10,
#                    box.padding = 0.3,
#                    label.padding = 0.1,
#                    point.padding = 0,
#                    size = 2) +
#   geom_vline(xintercept = 0.5,
#              linetype = "dashed",
#              col = "black") +
#   geom_hline(yintercept = 0.5,
#              linetype = "dashed",
#              col = "black") + 
#   xlab("Sensitivity (AUROC)") + 
#   ylab("Precision (AUPRC)") +
#   coord_cartesian(xlim = c(0.45, 0.75),
#                   ylim = c(0.45, 0.8)) +
#   theme(strip.text = element_text(colour = 'white'),
#         strip.background = element_rect(fill = "black"),
#         legend.position = "none",
#         panel.background = element_rect(fill = "white", colour = "black"),
#         panel.grid.major = element_line(colour = "grey", linewidth = 0.2),
#         panel.grid.minor = element_line(colour = "grey", linewidth = 0.1)) +
#   geom_point(data = benchRAW_HOTincl_shufflestats_df[!benchRAW_HOTincl_shufflestats_df$cutoff %in% c("1500", "3000", "10000") &
#                                                        benchRAW_HOTincl_shufflestats_df$method == "mlm_estimate", ], aes(x = mean_auroc, y = mean_auprc), alpha = 0.5) + 
#   geom_errorbar(data = benchRAW_HOTincl_shufflestats_df[!benchRAW_HOTincl_shufflestats_df$cutoff %in% c("1500", "3000", "10000") &
#                                                           benchRAW_HOTincl_shufflestats_df$method == "mlm_estimate", ], aes(x = mean_auroc, y = mean_auprc,
#                                                                                                                                                 ymin = (mean_auprc - sd_auprc),
#                                                                                                                                                 ymax = (mean_auprc + sd_auprc)), alpha = 0.5) +
#   geom_errorbar(data = benchRAW_HOTincl_shufflestats_df[!benchRAW_HOTincl_shufflestats_df$cutoff %in% c("1500", "3000", "10000") &
#                                                           benchRAW_HOTincl_shufflestats_df$method == "mlm_estimate", ], aes(x = mean_auroc, y = mean_auprc,
#                                                                                                                                                 xmin = (mean_auroc - sd_auroc),
#                                                                                                                                                 xmax = (mean_auroc + sd_auroc)), alpha = 0.5)
# dev.off()
# 
# 
# pdf("graphics/allChIP_HOTincl_benchRAPToR_mlmplot.pdf",
#     height = 2.5,
#     width = 2.5)
# 
# ggplot(data = allChIP_HOTincl_benchraptor_allmethods_plot[!str_detect(allChIP_HOTincl_benchraptor_allmethods_plot$net, "shuffle") & 
#                                                          !allChIP_HOTincl_benchraptor_allmethods_plot$cutoff %in% c("1500", "3000", "10000") &
#                                                          allChIP_HOTincl_benchraptor_allmethods_plot$method == "mlm_estimate", ], aes(x = auroc, y = auprc, colour = cutoff)) + 
#   geom_point(size = 2) +
#   geom_label_repel(aes(label = label),
#                    label.size = 0.1,
#                    max.overlaps = 10,
#                    box.padding = 0.3,
#                    label.padding = 0.1,
#                    point.padding = 0,
#                    size = 2) +
#   geom_vline(xintercept = 0.5,
#              linetype = "dashed",
#              col = "black") +
#   geom_hline(yintercept = 0.5,
#              linetype = "dashed",
#              col = "black") + 
#   xlab("Sensitivity (AUROC)") + 
#   ylab("Precision (AUPRC)") +
#   coord_cartesian(xlim = c(0.45, 0.75),
#                   ylim = c(0.45, 0.8)) +
#   theme(strip.text = element_text(colour = 'white'),
#         strip.background = element_rect(fill = "black"),
#         legend.position = "none",
#         panel.background = element_rect(fill = "white", colour = "black"),
#         panel.grid.major = element_line(colour = "grey", linewidth = 0.2),
#         panel.grid.minor = element_line(colour = "grey", linewidth = 0.1),
#         axis.text = element_text(colour = "black")) +
#   geom_point(data = benchraptor_HOTincl_shufflestats_df[!benchraptor_HOTincl_shufflestats_df$cutoff %in% c("1500", "3000", "10000") &
#                                                        benchraptor_HOTincl_shufflestats_df$method == "mlm_estimate", ], aes(x = mean_auroc, y = mean_auprc), alpha = 0.5) + 
#   geom_errorbar(data = benchraptor_HOTincl_shufflestats_df[!benchraptor_HOTincl_shufflestats_df$cutoff %in% c("1500", "3000", "10000") &
#                                                           benchraptor_HOTincl_shufflestats_df$method == "mlm_estimate", ], aes(x = mean_auroc, y = mean_auprc,
#                                                                                                                             ymin = (mean_auprc - sd_auprc),
#                                                                                                                             ymax = (mean_auprc + sd_auprc)), alpha = 0.5) +
#   geom_errorbar(data = benchraptor_HOTincl_shufflestats_df[!benchraptor_HOTincl_shufflestats_df$cutoff %in% c("1500", "3000", "10000") &
#                                                           benchraptor_HOTincl_shufflestats_df$method == "mlm_estimate", ], aes(x = mean_auroc, y = mean_auprc,
#                                                                                                                             xmin = (mean_auroc - sd_auroc),
#                                                                                                                             xmax = (mean_auroc + sd_auroc)), alpha = 0.5)
# dev.off()

#### benchraptor hot exclusion 

allChIP_HOTexcl_benchraptor_files <- bench_out_filelist[str_detect(bench_out_filelist, "benchRAPToR") & str_detect(bench_out_filelist, "HOTexcl")]

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

allChIP_HOTexcl_benchraptor_allmethods_plot[, "label"] <- allChIP_HOTexcl_benchraptor_allmethods_plot$net
# allChIP_HOTexcl_benchraptor_allmethods_plot[str_detect(unlist(allChIP_HOTexcl_benchraptor_allmethods_plot$net), "shuffle"), "label"] <- ""

allChIP_HOTexcl_benchraptor_allmethods_plot[, "cutoff"] <- str_remove(allChIP_HOTexcl_benchraptor_allmethods_plot$net, "_shuffle")
allChIP_HOTexcl_benchraptor_allmethods_plot[, "shuffled"] <- str_detect(allChIP_HOTexcl_benchraptor_allmethods_plot$net, "_shuffle")

allChIP_HOTexcl_benchraptor_allmethods_plot[allChIP_HOTexcl_benchraptor_allmethods_plot$label == "5000", "label"] <- "none"

benchraptor_HOTexcl_shufflestats_files <- shufflestats_files[str_detect(shufflestats_files, "RAPToR") & str_detect(shufflestats_files, "HOTexcl")]

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
  
  tempout["net"] <- paste0(str_extract(file, "[0-9]{3,4}"), "_shuffle")
  
  tempout
  
})

benchraptor_HOTexcl_shufflestats_df <- do.call(rbind, benchraptor_HOTexcl_shufflestats_lists)
benchraptor_HOTexcl_shufflestats_df[, "cutoff"] <- str_remove(benchraptor_HOTexcl_shufflestats_df$net, "_shuffle")

saveRDS(allChIP_HOTexcl_benchraptor_allmethods_plot,
        "plotdata/allChIP_HOTexcl_benchraptor_allmethods_plot.rds")
allChIP_HOTexcl_benchraptor_allmethods_plot <- readRDS("plotdata/allChIP_HOTexcl_benchraptor_allmethods_plot.rds")

saveRDS(benchraptor_HOTexcl_shufflestats_df,
        "plotdata/benchRAPToR_HOTexcl_shufflestats_df.rds")
benchraptor_HOTexcl_shufflestats_df <- readRDS("plotdata/benchRAPToR_HOTexcl_shufflestats_df.rds")

pdf("graphics/allChIP_HOTexcl_benchRAPToR_mlmplot.pdf",
    height = 2.5,
    width = 2.5)

ggplot(data = allChIP_HOTexcl_benchraptor_allmethods_plot[!str_detect(allChIP_HOTexcl_benchraptor_allmethods_plot$net, "shuffle") & 
                                                            !allChIP_HOTexcl_benchraptor_allmethods_plot$cutoff %in% c("1500", "3000", "10000") &
                                                            allChIP_HOTexcl_benchraptor_allmethods_plot$method == "mlm_estimate", ], aes(x = auroc, y = auprc, colour = cutoff)) + 
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
  geom_point(data = benchraptor_HOTexcl_shufflestats_df[!benchraptor_HOTexcl_shufflestats_df$cutoff %in% c("1500", "3000", "10000") &
                                                          benchraptor_HOTexcl_shufflestats_df$method == "mlm_estimate", ], aes(x = mean_auroc, y = mean_auprc), alpha = 0.5) + 
  geom_errorbar(data = benchraptor_HOTexcl_shufflestats_df[!benchraptor_HOTexcl_shufflestats_df$cutoff %in% c("1500", "3000", "10000") &
                                                             benchraptor_HOTexcl_shufflestats_df$method == "mlm_estimate", ], aes(x = mean_auroc, y = mean_auprc,
                                                                                                                                  ymin = (mean_auprc - sd_auprc),
                                                                                                                                  ymax = (mean_auprc + sd_auprc)), alpha = 0.5) +
  geom_errorbar(data = benchraptor_HOTexcl_shufflestats_df[!benchraptor_HOTexcl_shufflestats_df$cutoff %in% c("1500", "3000", "10000") &
                                                             benchraptor_HOTexcl_shufflestats_df$method == "mlm_estimate", ], aes(x = mean_auroc, y = mean_auprc,
                                                                                                                                  xmin = (mean_auroc - sd_auroc),
                                                                                                                                  xmax = (mean_auroc + sd_auroc)), alpha = 0.5)
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
  geom_label_repel(aes(label = label),
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
                                                                                                                                                      ymax = (mean_auprc + sd_auprc)), alpha = 0.5) +
  geom_errorbar(data = benchraptor_HOTexcl_shufflestats_df[!benchraptor_HOTexcl_shufflestats_df$cutoff %in% c("1500", "3000", "10000") &
                                                             !benchraptor_HOTexcl_shufflestats_df$method %in% c("wsum_estimate", "wsum_corr"), ], aes(x = mean_auroc, y = mean_auprc,
                                                                                                                                                      xmin = (mean_auroc - sd_auroc),
                                                                                                                                                      xmax = (mean_auroc + sd_auroc)), alpha = 0.5)

dev.off()
