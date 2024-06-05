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

combo_restrictedobs_labelvec <- c("PWMs", "ChIP", "no weight", "weighted", "uneq. weights")
names(combo_restrictedobs_labelvec) <- unique(combo_benchraptormethods_plot$net)

combo_benchraptormethods_plot[, "label"] <- combo_restrictedobs_labelvec[combo_benchraptormethods_plot$net]

saveRDS(combo_benchraptormethods_plot,
        "plotdata/combo_benchraptormethods_plot.rds")

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

pdf("graphics/combo_restrictedobs_benchRAPToR_methodsplot.pdf",
    height = 4,
    width = 5)

ggplot(data = combo_benchraptormethods_plot[!combo_benchraptormethods_plot$method %in% c("wsum_estimate", "wsum_corr"), ], aes(x = auroc, y = auprc, colour = net)) + 
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
                  ylim = c(0.45, 0.85)) +
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
                                                                                                                      ymax = (mean_auprc + sd_auprc)), alpha = 0.5) +
  geom_errorbar(data = combo_restrictedobs_benchraptor_shufflestats_df[!combo_restrictedobs_benchraptor_shufflestats_df$method %in% c("wsum_estimate", "wsum_corr"), ], aes(x = mean_auroc, y = mean_auprc,
                                                                                                                      xmin = (mean_auroc - sd_auroc),
                                                                                                                      xmax = (mean_auroc + sd_auroc)), alpha = 0.5)

dev.off()

pdf("graphics/combo_restrictedobs_benchraptor_mlmplot.pdf",
    height = 2.5,
    width = 2.5)

ggplot(data = combo_benchraptormethods_plot[!str_detect(combo_benchraptormethods_plot$net, "shuffle") & 
                                                        combo_benchraptormethods_plot$method == "mlm_estimate", ], aes(x = auroc, y = auprc, colour = net)) + 
  geom_point(size = 2) +
  geom_label_repel(aes(label = label),
                   label.size = 0.1,
                   max.overlaps = 10,
                   box.padding = 0,
                   label.padding = 0.1,
                   point.padding = 2,
                   size = 1.5) +
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
  geom_point(data = combo_restrictedobs_benchraptor_shufflestats_df[combo_restrictedobs_benchraptor_shufflestats_df$method == "mlm_estimate", ], aes(x = mean_auroc, y = mean_auprc), alpha = 0.5) + 
  geom_errorbar(data = combo_restrictedobs_benchraptor_shufflestats_df[combo_restrictedobs_benchraptor_shufflestats_df$method == "mlm_estimate", ], aes(x = mean_auroc, y = mean_auprc,
                                                                                                                                          ymin = (mean_auprc - sd_auprc),
                                                                                                                                          ymax = (mean_auprc + sd_auprc)), alpha = 0.5) +
  geom_errorbar(data = combo_restrictedobs_benchraptor_shufflestats_df[combo_restrictedobs_benchraptor_shufflestats_df$method == "mlm_estimate", ], aes(x = mean_auroc, y = mean_auprc,
                                                                                                                                          xmin = (mean_auroc - sd_auroc),
                                                                                                                                          xmax = (mean_auroc + sd_auroc)), alpha = 0.5)
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

pdf("graphics/combo3_benchRAPToR_methodsplot.pdf",
    height = 4,
    width = 5)

ggplot(data = combo3_benchraptormethods_plot[!combo3_benchraptormethods_plot$method %in% c("wsum_estimate", "wsum_corr") &
                                               combo3_benchraptormethods_plot$net != "combo3_unequal", ], aes(x = auroc, y = auprc, colour = net)) + 
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
  geom_point(data = combo3_benchraptor_shufflestats_df[!combo3_benchraptor_shufflestats_df$method %in% c("wsum_estimate", "wsum_corr"), ], aes(x = mean_auroc, y = mean_auprc), alpha = 0.5) + 
  geom_errorbar(data = combo3_benchraptor_shufflestats_df[!combo3_benchraptor_shufflestats_df$method %in% c("wsum_estimate", "wsum_corr"), ], aes(x = mean_auroc, y = mean_auprc,
                                                                                                                                                                            ymin = (mean_auprc - sd_auprc),
                                                                                                                                                                            ymax = (mean_auprc + sd_auprc)), alpha = 0.5) +
  geom_errorbar(data = combo3_benchraptor_shufflestats_df[!combo3_benchraptor_shufflestats_df$method %in% c("wsum_estimate", "wsum_corr"), ], aes(x = mean_auroc, y = mean_auprc,
                                                                                                                                                                            xmin = (mean_auroc - sd_auroc),
                                                                                                                                                                            xmax = (mean_auroc + sd_auroc)), alpha = 0.5)

dev.off()

pdf("graphics/combo3_benchraptor_mlmplot.pdf",
    height = 3.5,
    width = 3.5)

ggplot(data = combo3_benchraptormethods_plot[!str_detect(combo3_benchraptormethods_plot$net, "shuffle") & 
                                              combo3_benchraptormethods_plot$method == "mlm_estimate" &
                                               combo3_benchraptormethods_plot$net != "combo3_unequal", ], aes(x = auroc, y = auprc, colour = net)) + 
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
  geom_point(data = combo3_benchraptor_shufflestats_df[combo3_benchraptor_shufflestats_df$method == "mlm_estimate", ], aes(x = mean_auroc, y = mean_auprc), alpha = 0.5) + 
  geom_errorbar(data = combo3_benchraptor_shufflestats_df[combo3_benchraptor_shufflestats_df$method == "mlm_estimate", ], aes(x = mean_auroc, y = mean_auprc,
                                                                                                                                                        ymin = (mean_auprc - sd_auprc),
                                                                                                                                                        ymax = (mean_auprc + sd_auprc)), alpha = 0.5) +
  geom_errorbar(data = combo3_benchraptor_shufflestats_df[combo3_benchraptor_shufflestats_df$method == "mlm_estimate", ], aes(x = mean_auroc, y = mean_auprc,
                                                                                                                                                        xmin = (mean_auroc - sd_auroc),
                                                                                                                                                        xmax = (mean_auroc + sd_auroc)), alpha = 0.5)
dev.off()


#### NOW FOR COMBO FULL SET ####

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

pdf("graphics/combo3_benchraw_methodsplot.pdf",
    height = 4,
    width = 5)

ggplot(data = combo3_benchrawmethods_plot[!combo3_benchrawmethods_plot$method %in% c("wsum_estimate", "wsum_corr") &
                                               combo3_benchrawmethods_plot$net != "combo3_unequal", ], aes(x = auroc, y = auprc, colour = net)) + 
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
  geom_point(data = combo3_benchraw_shufflestats_df[!combo3_benchraw_shufflestats_df$method %in% c("wsum_estimate", "wsum_corr"), ], aes(x = mean_auroc, y = mean_auprc), alpha = 0.5) + 
  geom_errorbar(data = combo3_benchraw_shufflestats_df[!combo3_benchraw_shufflestats_df$method %in% c("wsum_estimate", "wsum_corr"), ], aes(x = mean_auroc, y = mean_auprc,
                                                                                                                                                  ymin = (mean_auprc - sd_auprc),
                                                                                                                                                  ymax = (mean_auprc + sd_auprc)), alpha = 0.5) +
  geom_errorbar(data = combo3_benchraw_shufflestats_df[!combo3_benchraw_shufflestats_df$method %in% c("wsum_estimate", "wsum_corr"), ], aes(x = mean_auroc, y = mean_auprc,
                                                                                                                                                  xmin = (mean_auroc - sd_auroc),
                                                                                                                                                  xmax = (mean_auroc + sd_auroc)), alpha = 0.5)

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
                                                                                                                              ymax = (mean_auprc + sd_auprc)), alpha = 0.5) +
  geom_errorbar(data = combo3_benchraw_shufflestats_df[combo3_benchraw_shufflestats_df$method == "mlm_estimate", ], aes(x = mean_auroc, y = mean_auprc,
                                                                                                                              xmin = (mean_auroc - sd_auroc),
                                                                                                                              xmax = (mean_auroc + sd_auroc)), alpha = 0.5)
dev.off()

