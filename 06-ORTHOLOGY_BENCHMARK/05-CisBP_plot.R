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
  tempGRN <- do.call(rbind, lapply(names(TF_orthology_prob_list), function(thisTF){
    # thisTF <- names(TF_orthology_prob_list)[5]
    x <- TF_orthology_prob_list[[thisTF]]
    
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
              "ggplot2",
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

#### MAIN CisBP orthology PLOT ####

TForth_benchraptor_shufflestat_files <- paste0(c("TForth", "CelEsT", "control", "FIMOfull"), "_benchRAPToR_shufflestats.tsv")

TForth_benchraptor <- read.table("output/benchmark_out/FIMO_TForth_df.tsv",
                                 sep= '\t',
                                 header = TRUE)

TForth_benchraptormethods_plot <- TForth_benchraptor[TForth_benchraptor$metric %in% c("auroc", "auprc"), c("net", "method", "metric", "score")]
TForth_benchraptormethods_plot <- tidyr::pivot_wider(TForth_benchraptormethods_plot, names_from = c("metric"), values_from = "score")

TForth_benchraptormethods_plot[, "label"] <- TForth_benchraptormethods_plot$net
TForth_benchraptormethods_plot[str_detect(unlist(TForth_benchraptormethods_plot$net), "shuffle"), "label"] <- ""

TForth_labelvec <- c("Conservation filtered", "control", "full 1500", "CelEsT*")
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

method_labels <- c("Consensus", "MLM", "ULM", "WSum")
names(method_labels) <- c("consensus_estimate", "mlm_estimate", "ulm_estimate", "wsum_norm")

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
