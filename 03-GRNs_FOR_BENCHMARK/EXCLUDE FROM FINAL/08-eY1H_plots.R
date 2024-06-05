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

eY1H_benchraptormethods_plot <- readRDS("plotdata/eY1H_benchraptormethods_plot.rds")

eY1H_shufflestats_df <- readRDS("plotdata/eY1H_benchraptor_shufflestats_df.rds")

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
                                                                                                                                                              ymax = (mean_auprc + sd_auprc)), alpha = 0.5) +
  geom_errorbar(data = eY1H_shufflestats_df[!eY1H_shufflestats_df$method %in% c("wsum_estimate", "wsum_corr"), ], aes(x = mean_auroc, y = mean_auprc,
                                                                                                                                                              xmin = (mean_auroc - sd_auroc),
                                                                                                                                        xmax = (mean_auroc + sd_auroc)), alpha = 0.5)
                                                                                                                                                              

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
