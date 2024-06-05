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

bench_out_filelist <- list.files("output/benchmark_out/")

shufflestats_files <- bench_out_filelist[str_detect(bench_out_filelist, "_shufflestats")]

FIMO_nohomo_benchraw_files <- bench_out_filelist[str_detect(bench_out_filelist, "benchRAW") & str_detect(bench_out_filelist, "nohomo")]

# FIMO_nohomo_benchraw_files <- FIMO_nohomo__benchraw_files[!str_detect(FIMO_nohomo__benchraw_files, "shufflestat")]

FIMO_nohomo_benchraw_list <- lapply(FIMO_nohomo_benchraw_files, function(x){
  
  read.table(paste0("output/benchmark_out/", x),
             sep = "\t",
             header = TRUE)
  
})

FIMO_nohomo_benchraw_all <- do.call(rbind, FIMO_nohomo_benchraw_list)

FIMO_nohomo_benchraw_allmethods_plot <- FIMO_nohomo_benchraw_all[FIMO_nohomo_benchraw_all$metric %in% c("auroc", "auprc"), c("net", "method", "metric", "score")]
FIMO_nohomo_benchraw_allmethods_plot <- tidyr::pivot_wider(FIMO_nohomo_benchraw_allmethods_plot, names_from = c("metric"), values_from = "score")

FIMO_nohomo_benchraw_allmethods_plot[, "label"] <- FIMO_nohomo_benchraw_allmethods_plot$net
FIMO_nohomo_benchraw_allmethods_plot[str_detect(unlist(FIMO_nohomo_benchraw_allmethods_plot$net), "shuffle"), "label"] <- ""

FIMO_nohomo_benchraw_allmethods_plot[, "cutoff"] <- str_remove(FIMO_nohomo_benchraw_allmethods_plot$net, "_shuffle")
FIMO_nohomo_benchraw_allmethods_plot[, "shuffled"] <- str_detect(FIMO_nohomo_benchraw_allmethods_plot$net, "_shuffle")

ggplot(data = FIMO_nohomo_benchraw_allmethods_plot, aes(x = auroc, y = auprc, label = label, alpha = shuffled, color = cutoff)) + 
  geom_point() + 
  geom_label_repel(label.size = 0.1,
                   max.overlaps = 50,
                   box.padding = 0.1,
                   label.padding = 0.1,
                   point.padding = 1) + 
  scale_alpha_discrete(range = c(1, 0.2)) +
  theme_classic() + 
  geom_vline(xintercept = 0.5,
             linetype = "dashed",
             col = "grey") +
  geom_hline(yintercept = 0.5,
             linetype = "dashed",
             col = "grey") + 
  facet_wrap(~ method)

ggplot(data = FIMO_nohomo_benchraw_allmethods_plot[FIMO_nohomo_benchraw_allmethods_plot$method == "mlm_estimate", ], aes(x = auroc, y = auprc, label = label, alpha = shuffled, color = cutoff)) + 
  geom_point() + 
  geom_label_repel(label.size = 0.1,
                   max.overlaps = 50,
                   box.padding = 0.1,
                   label.padding = 0.1,
                   point.padding = 1) + 
  scale_alpha_discrete(range = c(1, 0.2)) +
  theme_classic() + 
  geom_vline(xintercept = 0.5,
             linetype = "dashed",
             col = "grey") +
  geom_hline(yintercept = 0.5,
             linetype = "dashed",
             col = "grey")

plot(as.numeric(unlist(FIMO_nohomo_benchraw_allmethods_plot[FIMO_nohomo_benchraw_allmethods_plot$shuffled == FALSE & FIMO_nohomo_benchraw_allmethods_plot$method == "mlm_estimate", "cutoff"])),
     unlist(FIMO_nohomo_benchraw_allmethods_plot[FIMO_nohomo_benchraw_allmethods_plot$shuffled == FALSE & FIMO_nohomo_benchraw_allmethods_plot$method == "mlm_estimate", "auprc"]))

ggplot(data = FIMO_nohomo_benchraw_allmethods_plot[FIMO_nohomo_benchraw_allmethods_plot$shuffled == FALSE, ], aes(x = as.numeric(cutoff), y = auprc, colour = method)) + 
  geom_point() + 
  geom_smooth(se = FALSE) +
  scale_x_log10() + 
  xlab(Log[10]~"cutoff") + 
  ylab("AUPRC")+ 
  coord_cartesian(ylim = c(0.55, 0.7))


#### RAPToR

FIMO_nohomo_benchraptor_files <- bench_out_filelist[str_detect(bench_out_filelist, "benchRAPToR") & str_detect(bench_out_filelist, "nohomo")]

FIMO_nohomo_benchraptor_files <- FIMO_nohomo_benchraptor_files[!str_detect(FIMO_nohomo_benchraptor_files, "shufflestat")]

FIMO_nohomo_benchraptor_list <- lapply(FIMO_nohomo_benchraptor_files, function(x){
  
  read.table(paste0("output/benchmark_out/", x),
             sep = "\t",
             header = TRUE)
  
})

FIMO_nohomo_benchraptor_all <- do.call(rbind, FIMO_nohomo_benchraptor_list)

FIMO_nohomo_benchraptor_allmethods_plot <- FIMO_nohomo_benchraptor_all[FIMO_nohomo_benchraptor_all$metric %in% c("auroc", "auprc"), c("net", "method", "metric", "score")]
FIMO_nohomo_benchraptor_allmethods_plot <- tidyr::pivot_wider(FIMO_nohomo_benchraptor_allmethods_plot, names_from = c("metric"), values_from = "score")

FIMO_nohomo_benchraptor_allmethods_plot[, "label"] <- FIMO_nohomo_benchraptor_allmethods_plot$net
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
  
  tempout["net"] <- paste0(str_extract(file, "[0-9]{3,4}"), "_shuffle")
  
  tempout
  
})

FIMO_nohomo_benchraptor_shufflestats_df <- do.call(rbind, FIMO_nohomo_benchraptor_shufflestats_lists)
FIMO_nohomo_benchraptor_shufflestats_df[, "cutoff"] <- str_remove(FIMO_nohomo_benchraptor_shufflestats_df$net, "_shuffle")

method_labels <- c("Consensus", "MLM", "ULM", "WSum")
names(method_labels) <- c("consensus_estimate", "mlm_estimate", "ulm_estimate", "wsum_norm")

saveRDS(FIMO_nohomo_benchraptor_allmethods_plot,
        "plotdata/FIMO_nohomo_benchraptor_allmethods_plot.rds")
FIMO_nohomo_benchraptor_allmethods_plot <- readRDS("plotdata/FIMO_nohomo_benchraptor_allmethods_plot.rds")

saveRDS(FIMO_nohomo_benchraptor_shufflestats_df,
        "plotdata/FIMO_nohomo_benchraptor_shufflestats_df.rds")
FIMO_nohomo_benchraptor_shufflestats_df <- readRDS("plotdata/FIMO_nohomo_benchraptor_shufflestats_df.rds")

pdf("graphics/FIMO_nohomo_benchraptor_methodsplot.pdf",
    height = 4,
    width = 5)

ggplot(data = FIMO_nohomo_benchraptor_allmethods_plot[!str_detect(FIMO_nohomo_benchraptor_allmethods_plot$net, "shuffle") & 
                                                            !FIMO_nohomo_benchraptor_allmethods_plot$cutoff %in% c("2500", "10000") &
                                                            !FIMO_nohomo_benchraptor_allmethods_plot$method %in% c("wsum_estimate", "wsum_corr"), ], aes(x = auroc, y = auprc, colour = cutoff)) + 
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
  geom_point(data = FIMO_nohomo_benchraptor_shufflestats_df[!FIMO_nohomo_benchraptor_shufflestats_df$cutoff %in% c("2500", "10000") &
                                                          !FIMO_nohomo_benchraptor_shufflestats_df$method %in% c("wsum_estimate", "wsum_corr"), ], aes(x = mean_auroc, y = mean_auprc), alpha = 0.5) + 
  geom_errorbar(data = FIMO_nohomo_benchraptor_shufflestats_df[!FIMO_nohomo_benchraptor_shufflestats_df$cutoff %in% c("2500", "10000") &
                                                             !FIMO_nohomo_benchraptor_shufflestats_df$method %in% c("wsum_estimate", "wsum_corr"), ], aes(x = mean_auroc, y = mean_auprc,
                                                                                                                                                      ymin = (mean_auprc - sd_auprc),
                                                                                                                                                      ymax = (mean_auprc + sd_auprc)), alpha = 0.5) +
  geom_errorbar(data = FIMO_nohomo_benchraptor_shufflestats_df[!FIMO_nohomo_benchraptor_shufflestats_df$cutoff %in% c("2500", "10000") &
                                                             !FIMO_nohomo_benchraptor_shufflestats_df$method %in% c("wsum_estimate", "wsum_corr"), ], aes(x = mean_auroc, y = mean_auprc,
                                                                                                                                                      xmin = (mean_auroc - sd_auroc),
                                                                                                                                                      xmax = (mean_auroc + sd_auroc)), alpha = 0.5)

dev.off()

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

# 
# plot(as.numeric(unlist(FIMO_nohomo_benchraptor_allmethods_plot[FIMO_nohomo_benchraptor_allmethods_plot$shuffled == FALSE & FIMO_nohomo_benchraptor_allmethods_plot$method == "mlm_estimate", "cutoff"])),
#      unlist(FIMO_nohomo_benchraptor_allmethods_plot[FIMO_nohomo_benchraptor_allmethods_plot$shuffled == FALSE & FIMO_nohomo_benchraptor_allmethods_plot$method == "mlm_estimate", "auprc"]))
# 
# unique(FIMO_nohomo_benchraptor_allmethods_plot$method)
# methods_label <- c("MLM", "ULM", "WSUM_est", "WSUM_norm", "WSum", "Consensus")
# names(methods_label) <- unique(FIMO_nohomo_benchraptor_allmethods_plot$method)
# 
# FIMO_nohomo_benchraptor_allmethods_plot$method <- methods_label[FIMO_nohomo_benchraptor_allmethods_plot$method]

saveRDS(FIMO_nohomo_benchraptor_allmethods_plot,
        "plotdata/FIMO_nohomo_benchraptor_allmethods_plot.rds")

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

#### homotypic 

param5 <- 5


FIMO_param5_benchraptor_files <- bench_out_filelist[str_detect(bench_out_filelist, "benchRAPToR") & str_detect(bench_out_filelist, paste0("param", param5))]

# FIMO_nohomo_benchraw_files <- FIMO_nohomo__benchraw_files[!str_detect(FIMO_nohomo__benchraw_files, "shufflestat")]

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
                       colour = param),
                   label.size = 0.1,
                   max.overlaps = 50,
                   box.padding = 0.1,
                   label.padding = 0.1,
                   point.padding = 1,
                   size = 3) +
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
                                                                                                                                          ymax = (mean_auprc + sd_auprc)), alpha = 0.5) +
  geom_errorbar(data = FIMO_param5_benchraptor_shufflestats_df[FIMO_param5_benchraptor_shufflestats_df$method == "mlm_estimate", ], aes(x = mean_auroc, y = mean_auprc,
                                                                                                                                          xmin = (mean_auroc - sd_auroc),
                                                                                                                           xmax = (mean_auroc + sd_auroc)), alpha = 0.5)
                                                                                                                                          
dev.off()

