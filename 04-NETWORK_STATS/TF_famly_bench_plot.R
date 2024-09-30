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

#### LOAD PACKAGES & FUNCTIONS ####

## First specify the packages of interest

packages <- c("openxlsx",
              "dplyr",
              "stringr",
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

github_packages <- c()

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

wTF3 <- readRDS("output/wTF3_modified.rds")

benchmark_obs <- read.table("output/benchmark_observations.txt",
                            sep = "\t",
                            header = TRUE)

CelEsT <- read.table("output/GRNs/allthree_equalweights.txt",
                     sep = "\t",
                     header = TRUE)

#### PLOT BY FAMILY ####

wTF3[wTF3$family == "AT Hook" & wTF3$bench_no != 0, ]

main_families <- c("ZF - NHR",
                   "ZF - C2H2",
                   "HD"
                   ,
                   "bHLH",
                   "WH",
                   "bZIP"
                   )

family_stats_list <- lapply(main_families, function(fam){
  
  temp_out <- read.table(paste0("output/benchmark_out/FAMILY_", fam,"_benchRAW.tsv"),
                         header = TRUE,
                         sep = "\t")
  
  temp_shuffle <- read.table(paste0("output/benchmark_out/", fam,"_benchRAW_shufflestats.tsv"),
                             header = TRUE,
                             sep = "\t")
  
  temp_shuffle <- temp_shuffle[str_detect(temp_shuffle$net, "shuffle"), 4:ncol(temp_shuffle)]
  
  temp_shuffle_auprc <- temp_shuffle[temp_shuffle$method == "mlm_estimate" & temp_shuffle$metric == "auprc", ]
  temp_shuffle_auroc <- temp_shuffle[temp_shuffle$method == "mlm_estimate" & temp_shuffle$metric == "auroc", ]
  
  output_vec <- c("auprc" = temp_out[temp_out$method == "mlm_estimate" & temp_out$metric == "auprc", "score"],
                  "random_mean_auprc" = mean(temp_shuffle_auprc$score),
                  "sd_auprc" = sd(temp_shuffle_auprc$score),
                  "auroc" = temp_out[temp_out$method == "mlm_estimate" & temp_out$metric == "auroc", "score"],
                  "random_mean_auroc" = mean(temp_shuffle_auroc$score),
                  "sd_auroc" = sd(temp_shuffle_auroc$score))

  output_vec
  
})


names(family_stats_list) <- main_families

family_stats_plot <- do.call(rbind, family_stats_list)
family_stats_plot <- as.data.frame(family_stats_plot)

family_stats_plot[, "group"] <- row.names(family_stats_plot)

family_stats_plot[, "auprc_95CI_margin"] <- qt(0.975, df = 99) * family_stats_plot$sd_auprc / sqrt(99)
family_stats_plot[, "auroc_95CI_margin"] <- qt(0.975, df = 99) * family_stats_plot$sd_auroc / sqrt(99)

family_stats_plot[, "label"] <- row.names(family_stats_plot)
family_stats_plot[, "label"] <- str_remove_all(family_stats_plot[, "label"], " ")

family_stats_plot$label  <- factor(family_stats_plot$label, levels = c("ZF-NHR", "ZF-C2H2", "HD", "bHLH", "WH", "bZIP"))

# annotate with number of benchmarking experiments for point size
family_bench_no <- sapply(main_families, function(thisfam){
  
  sum(wTF3[wTF3$family == thisfam, "bench_no"])
  
})

names(family_bench_no) <- str_remove_all(names(family_bench_no), " ")

family_stats_plot[, "bench_no"] <- family_bench_no[family_stats_plot$label]

family_TF_no <- sapply(main_families, function(thisfam){

  sum((wTF3[wTF3$family == thisfam, "Sequence.name"] %in% benchmark_obs$target_gseq) & (wTF3[wTF3$family == thisfam, "Sequence.name"] %in% unique(CelEsT$source)))
  
})

names(family_TF_no) <- str_remove_all(names(family_TF_no), " ")

family_stats_plot[, "TF_no"] <- family_TF_no[family_stats_plot$label]

saveRDS(family_stats_plot,
        "plotdata/family_bench_plot.pdf")

pdf("graphics/bench_by_family.pdf",
    width = 2.5,
    height = 2.5)

ggplot(family_stats_plot, aes(x = auroc, y = auprc, colour = label)) + 
  geom_point(aes(size = bench_no)) + 
  geom_point(aes(size = TF_no), 
               colour = "white") +
  xlab("") + 
  ylab("") +
  coord_cartesian(xlim = c(0.2, 0.85),
                  ylim = c(0.4, 0.9)) +
  theme_classic() +
  geom_label_repel(aes(label = label,
                       fill = label),
                   colour = "black",
                   label.size = 0.1,
                   max.overlaps = 15,
                   box.padding = 0.3,
                   label.padding = 0.1,
                   point.padding = 0,
                   size = 2) +
  scale_fill_manual(values = c("magenta", "olivedrab3", "dodger blue", "yellow3", "lightsalmon1", "orchid1")) + 
  scale_color_manual(values = c("magenta", "olivedrab3", "dodger blue", "yellow3", "lightsalmon1", "orchid1")) +
  theme(strip.text = element_text(colour = 'white'),
        strip.background = element_rect(fill = "black"),
        legend.position = "none",
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_line(colour = "grey", linewidth = 0.2),
        panel.grid.minor = element_line(colour = "grey", linewidth = 0.1),
        axis.text = element_text(colour = "black")) +
geom_errorbar(data = family_stats_plot, aes(x = random_mean_auroc, y = random_mean_auprc,
                                            ymin = (random_mean_auprc - sd_auprc),
                                            ymax = (random_mean_auprc + sd_auprc)), alpha = 0.5, width = 0.01) +
  geom_errorbar(data = family_stats_plot, aes(x = random_mean_auroc, y = random_mean_auprc,
                                              xmin = (random_mean_auroc - sd_auroc),
                                              xmax = (random_mean_auroc + sd_auroc)), alpha = 0.5, width = 0.01) +
  geom_vline(xintercept = 0.5,
             linetype = "dashed",
             col = "black") +
  geom_hline(yintercept = 0.5,
             linetype = "dashed",
             col = "black")
  
dev.off()

family_statsRAPToR_list <- lapply(main_families, function(fam){
  
  temp_out <- read.table(paste0("output/benchmark_out/FAMILY_", fam,"_benchRAPToR.tsv"),
                         header = TRUE,
                         sep = "\t")
  
  temp_shuffle <- read.table(paste0("output/benchmark_out/", fam,"_benchRAPToR_shufflestats.tsv"),
                             header = TRUE,
                             sep = "\t")
  
  temp_shuffle <- temp_shuffle[str_detect(temp_shuffle$net, "shuffle"), 4:ncol(temp_shuffle)]
  
  temp_shuffle_auprc <- temp_shuffle[temp_shuffle$method == "mlm_estimate" & temp_shuffle$metric == "auprc", ]
  temp_shuffle_auroc <- temp_shuffle[temp_shuffle$method == "mlm_estimate" & temp_shuffle$metric == "auroc", ]
  
  output_vec <- c("auprc" = temp_out[temp_out$method == "mlm_estimate" & temp_out$metric == "auprc", "score"],
                  "random_mean_auprc" = mean(temp_shuffle_auprc$score),
                  "sd_auprc" = sd(temp_shuffle_auprc$score),
                  "auroc" = temp_out[temp_out$method == "mlm_estimate" & temp_out$metric == "auroc", "score"],
                  "random_mean_auroc" = mean(temp_shuffle_auroc$score),
                  "sd_auroc" = sd(temp_shuffle_auroc$score))
  
  output_vec
  
})


names(family_statsRAPToR_list) <- main_families

family_statsRAPToR_plot <- do.call(rbind, family_statsRAPToR_list)
family_statsRAPToR_plot <- as.data.frame(family_statsRAPToR_plot)

family_statsRAPToR_plot[, "group"] <- row.names(family_statsRAPToR_plot)

family_statsRAPToR_plot[, "auprc_95CI_margin"] <- qt(0.975, df = 99) * family_statsRAPToR_plot$sd_auprc / sqrt(99)
family_statsRAPToR_plot[, "auroc_95CI_margin"] <- qt(0.975, df = 99) * family_statsRAPToR_plot$sd_auroc / sqrt(99)

family_statsRAPToR_plot[, "label"] <- row.names(family_statsRAPToR_plot)
family_statsRAPToR_plot[, "label"] <- str_remove_all(family_statsRAPToR_plot[, "label"], " ")

family_statsRAPToR_plot$label  <- factor(family_statsRAPToR_plot$label, levels = c("ZF-NHR", "ZF-C2H2", "HD", "bHLH", "WH", "bZIP"))

# annotate with number of benchmarking experiments for point size

family_statsRAPToR_plot[, "bench_no"] <- family_bench_no[family_statsRAPToR_plot$label]

family_statsRAPToR_plot[, "TF_no"] <- family_TF_no[family_statsRAPToR_plot$label]

saveRDS(family_statsRAPToR_plot,
        "plotdata/family_benchRAPToR_plot.pdf")

pdf("graphics/bench_by_family_RAPToR.pdf",
    width = 2.5,
    height = 2.5)

ggplot(family_statsRAPToR_plot, aes(x = auroc, y = auprc, colour = label)) + 
  geom_point(aes(size = bench_no)) + 
  geom_point(aes(size = TF_no), 
             colour = "white") + 
  xlab("") + 
  ylab("") +
  coord_cartesian(xlim = c(0.2, 0.85),
                  ylim = c(0.4, 0.9)) +
  theme_classic() +
  geom_label_repel(aes(label = label,
                       fill = label),
                   colour = "black",
                   label.size = 0.1,
                   max.overlaps = 15,
                   box.padding = 0.3,
                   label.padding = 0.1,
                   point.padding = 0,
                   size = 2) +
  scale_fill_manual(values = c("magenta", "olivedrab3", "dodger blue", "yellow3", "lightsalmon1", "orchid1")) + 
  scale_color_manual(values = c("magenta", "olivedrab3", "dodger blue", "yellow3", "lightsalmon1", "orchid1")) +
  theme(strip.text = element_text(colour = 'white'),
        strip.background = element_rect(fill = "black"),
        legend.position = "none",
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_line(colour = "grey", linewidth = 0.2),
        panel.grid.minor = element_line(colour = "grey", linewidth = 0.1),
        axis.text = element_text(colour = "black")) +
  geom_errorbar(data = family_statsRAPToR_plot, aes(x = random_mean_auroc, y = random_mean_auprc,
                                              ymin = (random_mean_auprc - sd_auprc),
                                              ymax = (random_mean_auprc + sd_auprc)), alpha = 0.5, width = 0.01) +
  geom_errorbar(data = family_statsRAPToR_plot, aes(x = random_mean_auroc, y = random_mean_auprc,
                                              xmin = (random_mean_auroc - sd_auroc),
                                              xmax = (random_mean_auroc + sd_auroc)), alpha = 0.5, width = 0.01) +
  geom_vline(xintercept = 0.5,
             linetype = "dashed",
             col = "black") +
  geom_hline(yintercept = 0.5,
             linetype = "dashed",
             col = "black")

dev.off()

