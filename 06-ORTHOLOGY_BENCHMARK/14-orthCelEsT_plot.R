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

packages <- c("stringr",
              "dplyr",
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

orthCelEsT_benchraptor <- read.table("output/benchmark_out/orthCelEsT_df.tsv",
                            sep = "\t",
                            header = TRUE)

orthCelEsT_benchraptor_shufflestat_files <- paste0(c("orthCelEsT", "CelEsT"), "_benchRAPToR_shufflestats.tsv")

maxCelEsT_benchraptor <- read.table("output/benchmark_out/maxCelEsT_df.tsv",
                                     sep = "\t",
                                     header = TRUE)

maxCelEsT_benchraptor_shufflestat_files <- paste0(c("maxCelEsT", "CelEsT"), "_benchRAPToR_shufflestats.tsv")

#### PREPARE FOR PLOTTING ####

orthCelEsT_benchraptormethods_plot <- orthCelEsT_benchraptor[orthCelEsT_benchraptor$metric %in% c("auroc", "auprc"), c("net", "method", "metric", "score")]
orthCelEsT_benchraptormethods_plot <- tidyr::pivot_wider(orthCelEsT_benchraptormethods_plot, names_from = c("metric"), values_from = "score")

orthCelEsT_benchraptormethods_plot$net <- str_remove(orthCelEsT_benchraptormethods_plot$net, " $")

orthCelEsT_benchraptormethods_plot[, "label"] <- orthCelEsT_benchraptormethods_plot$net
orthCelEsT_benchraptormethods_plot[str_detect(unlist(orthCelEsT_benchraptormethods_plot$net), "shuffle"), "label"] <- ""

orthCelEsT_benchraptor_shufflestats_lists <- lapply(orthCelEsT_benchraptor_shufflestat_files, function(file){

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

orthCelEsT_benchraptor_shufflestats_df <- do.call(rbind, orthCelEsT_benchraptor_shufflestats_lists)

orthCelEsT_benchraptor_shufflestats_df[, "label"] <- orthCelEsT_benchraptor_shufflestats_df$net

saveRDS(orthCelEsT_benchraptormethods_plot,
        "plotdata/orthCelEsT_benchraptormethods_plot.rds")

saveRDS(orthCelEsT_benchraptor_shufflestats_df,
        "plotdata/orthCelEsT_benchraptor_shufflestats_df.rds")

method_labels <- c("Consensus", "MLM", "ULM", "WSum")
names(method_labels) <- c("consensus_estimate", "mlm_estimate", "ulm_estimate", "wsum_norm")

pdf("graphics/orthCelEsT_benchRAPToR_methodsplot.pdf",
    height = 4,
    width = 5)

ggplot(data = orthCelEsT_benchraptormethods_plot[!orthCelEsT_benchraptormethods_plot$method %in% c("wsum_estimate", "wsum_corr") &
                                               orthCelEsT_benchraptormethods_plot$net != "orthCelEsT_unequal", ], aes(x = auroc, y = auprc, colour = as.factor(label))) + 
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
  facet_wrap(~ method,
             labeller = labeller(method = method_labels)
  ) + 
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
  geom_point(data = orthCelEsT_benchraptor_shufflestats_df[!orthCelEsT_benchraptor_shufflestats_df$method %in% c("wsum_estimate", "wsum_corr"), ], aes(x = mean_auroc, y = mean_auprc), alpha = 0.5) + 
  geom_errorbar(data = orthCelEsT_benchraptor_shufflestats_df[!orthCelEsT_benchraptor_shufflestats_df$method %in% c("wsum_estimate", "wsum_corr"), ], aes(x = mean_auroc, y = mean_auprc,
                                                                                                                                                  ymin = (mean_auprc - sd_auprc),
                                                                                                                                                  ymax = (mean_auprc + sd_auprc)), alpha = 0.5) +
  geom_errorbar(data = orthCelEsT_benchraptor_shufflestats_df[!orthCelEsT_benchraptor_shufflestats_df$method %in% c("wsum_estimate", "wsum_corr"), ], aes(x = mean_auroc, y = mean_auprc,
                                                                                                                                                  xmin = (mean_auroc - sd_auroc),
                                                                                                                                                  xmax = (mean_auroc + sd_auroc),
                                                                                                                                                  color = label), alpha = 0.5)

dev.off()

pdf("graphics/orthCelEsT_benchraptor_mlmplot.pdf",
    height = 3.5,
    width = 3.5)

ggplot(data = orthCelEsT_benchraptormethods_plot[!str_detect(orthCelEsT_benchraptormethods_plot$net, "shuffle") & 
                                               orthCelEsT_benchraptormethods_plot$method == "mlm_estimate" &
                                               orthCelEsT_benchraptormethods_plot$net != "orthCelEsT_unequal", ], aes(x = auroc, y = auprc, colour = net)) + 
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
  geom_point(data = orthCelEsT_benchraptor_shufflestats_df[orthCelEsT_benchraptor_shufflestats_df$method == "mlm_estimate", ], aes(x = mean_auroc, y = mean_auprc), alpha = 0.5) + 
  geom_errorbar(data = orthCelEsT_benchraptor_shufflestats_df[orthCelEsT_benchraptor_shufflestats_df$method == "mlm_estimate", ], aes(x = mean_auroc, y = mean_auprc,
                                                                                                                              ymin = (mean_auprc - sd_auprc),
                                                                                                                              ymax = (mean_auprc + sd_auprc)), alpha = 0.5) +
  geom_errorbar(data = orthCelEsT_benchraptor_shufflestats_df[orthCelEsT_benchraptor_shufflestats_df$method == "mlm_estimate", ], aes(x = mean_auroc, y = mean_auprc,
                                                                                                                              xmin = (mean_auroc - sd_auroc),
                                                                                                                              xmax = (mean_auroc + sd_auroc)), alpha = 0.5)
dev.off()

#### maxCelEsT ####

maxCelEsT_benchraptormethods_plot <- maxCelEsT_benchraptor[maxCelEsT_benchraptor$metric %in% c("auroc", "auprc"), c("net", "method", "metric", "score")]
maxCelEsT_benchraptormethods_plot <- tidyr::pivot_wider(maxCelEsT_benchraptormethods_plot, names_from = c("metric"), values_from = "score")

maxCelEsT_benchraptormethods_plot$net <- str_remove(maxCelEsT_benchraptormethods_plot$net, " $")

maxCelEsT_benchraptormethods_plot[, "label"] <- maxCelEsT_benchraptormethods_plot$net
maxCelEsT_benchraptormethods_plot[str_detect(unlist(maxCelEsT_benchraptormethods_plot$net), "shuffle"), "label"] <- ""

maxCelEsT_benchraptor_shufflestats_lists <- lapply(maxCelEsT_benchraptor_shufflestat_files, function(file){
  
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

maxCelEsT_benchraptor_shufflestats_df <- do.call(rbind, maxCelEsT_benchraptor_shufflestats_lists)

maxCelEsT_benchraptor_shufflestats_df[, "label"] <- maxCelEsT_benchraptor_shufflestats_df$net

saveRDS(maxCelEsT_benchraptormethods_plot,
        "plotdata/maxCelEsT_benchraptormethods_plot.rds")

saveRDS(maxCelEsT_benchraptor_shufflestats_df,
        "plotdata/maxCelEsT_benchraptor_shufflestats_df.rds")

method_labels <- c("Consensus", "MLM", "ULM", "WSum")
names(method_labels) <- c("consensus_estimate", "mlm_estimate", "ulm_estimate", "wsum_norm")

pdf("graphics/maxCelEsT_benchRAPToR_methodsplot.pdf",
    height = 4,
    width = 5)

ggplot(data = maxCelEsT_benchraptormethods_plot[!maxCelEsT_benchraptormethods_plot$method %in% c("wsum_estimate", "wsum_corr") &
                                                   maxCelEsT_benchraptormethods_plot$net != "maxCelEsT_unequal", ], aes(x = auroc, y = auprc, colour = as.factor(label))) + 
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
  facet_wrap(~ method,
             labeller = labeller(method = method_labels)
  ) + 
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
  geom_point(data = maxCelEsT_benchraptor_shufflestats_df[!maxCelEsT_benchraptor_shufflestats_df$method %in% c("wsum_estimate", "wsum_corr"), ], aes(x = mean_auroc, y = mean_auprc), alpha = 0.5) + 
  geom_errorbar(data = maxCelEsT_benchraptor_shufflestats_df[!maxCelEsT_benchraptor_shufflestats_df$method %in% c("wsum_estimate", "wsum_corr"), ], aes(x = mean_auroc, y = mean_auprc,
                                                                                                                                                          ymin = (mean_auprc - sd_auprc),
                                                                                                                                                          ymax = (mean_auprc + sd_auprc)), alpha = 0.5) +
  geom_errorbar(data = maxCelEsT_benchraptor_shufflestats_df[!maxCelEsT_benchraptor_shufflestats_df$method %in% c("wsum_estimate", "wsum_corr"), ], aes(x = mean_auroc, y = mean_auprc,
                                                                                                                                                          xmin = (mean_auroc - sd_auroc),
                                                                                                                                                          xmax = (mean_auroc + sd_auroc),
                                                                                                                                                          color = label), alpha = 0.5)

dev.off()

