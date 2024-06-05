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
              "ggplot2")

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

allChIP_HOTexcl_benchraptor_allmethods_plot <- readRDS("plotdata/allChIP_HOTexcl_benchraptor_allmethods_plot.rds")
allChIP_HOTincl_benchraptor_allmethods_plot <- readRDS("plotdata/allChIP_HOTincl_benchraptor_allmethods_plot.rds")

#### PREPARE FOR PLOTTING ####

## process networks with STREME-based orthology profiling

ChIPorth_HOTexcl_benchraptor_files <- bench_out_filelist[str_detect(bench_out_filelist, "benchRAPToR") &
                                                           str_detect(bench_out_filelist, 'HOTexcl') &
                                                              str_detect(bench_out_filelist, "ChIPorth")]

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

## process networks with CisBP motifs instead of STREME motifs where applicable

STREMECisBP_HOTexcl_benchraptor_files <- bench_out_filelist[str_detect(bench_out_filelist, "benchRAPToR") &
                                                              !str_detect(bench_out_filelist, "HOTincl") &
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

## add network with signal-based order for comparison

allChIP_HOTexcl_benchraptor_allmethods_plot <- allChIP_HOTexcl_benchraptor_allmethods_plot[allChIP_HOTexcl_benchraptor_allmethods_plot$method == "mlm_estimate", ]

allChIP_HOTexcl_benchraptor_allmethods_plot <- allChIP_HOTexcl_benchraptor_allmethods_plot[, colnames(allChIP_HOTexcl_benchraptor_allmethods_plot) %in% colnames(STREMECisBP_HOTexcl_benchraptor_allmethods_plot)]

ChIPorth_HOTexcl_benchraptor_allmethods_plot[, "order"] <- "STREME motif conservation"
allChIP_HOTexcl_benchraptor_allmethods_plot[, "order"] <- "ChIP signal"
STREMECisBP_HOTexcl_benchraptor_allmethods_plot[, "order"] <- "STREME+CisBP motif conservation"

three_HOTexcl_benchraptor_allmethods_plot <- rbind(allChIP_HOTexcl_benchraptor_allmethods_plot,
                                                   STREMECisBP_HOTexcl_benchraptor_allmethods_plot,
                                                   ChIPorth_HOTexcl_benchraptor_allmethods_plot)

three_HOTexcl_benchraptor_allmethods_plot$order <- factor(three_HOTexcl_benchraptor_allmethods_plot$order, 
                                                             levels = c("ChIP signal", "STREME+CisBP motif conservation", "STREME motif conservation"))

saveRDS(three_HOTexcl_benchraptor_allmethods_plot,
        "plotdata/three_HOTexcl_benchraptor_allmethods_plot.rds")

three_HOTexcl_benchraptor_allmethods_plot <- readRDS("plotdata/three_HOTexcl_benchraptor_allmethods_plot.rds")

pdf("graphics/ChIPorth_allnetlines_HOTexcl_AUPRC.pdf",
    width = 2.5,
    height =2.5)

ggplot(data = three_HOTexcl_benchraptor_allmethods_plot, aes(x = as.numeric(net), y = auprc, colour = order)) + 
  geom_point() + 
  geom_line() + 
  scale_color_manual(values = c("dark green", "magenta", "blue")) +
  scale_x_log10(breaks = c(100, 500, 1000, 1500, 2000, 2500, 3000, 5000)) +
  ylab("Precision (AUPRC)") +
  xlab("Cutoff (log"[10]~"scale)") +
  coord_cartesian(ylim = c(0.55, 0.75)) +
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

pdf("graphics/ChIPorth_allnetlines_HOTexcl_AUROC.pdf",
    width = 2.5,
    height =2.5)

ggplot(data = three_HOTexcl_benchraptor_allmethods_plot, aes(x = as.numeric(net), y = auroc, colour = order)) + 
  geom_point() + 
  geom_line() + 
  scale_color_manual(values = c("dark green", "magenta", "blue")) +
  scale_x_log10(breaks = c(100, 500, 1000, 1500, 2000, 2500, 3000, 5000)) +
  ylab("Sensitivity (AUROC)") +
  xlab("Cutoff (log"[10]~"scale)") +
  coord_cartesian(ylim = c(0.55, 0.75)) +
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

pdf("graphics/ChIPorth_allnetlines_HOTexcl_AUPRC_LEGEND.pdf",
    width = 5,
    height = 5)

ggplot(data = three_HOTexcl_benchraptor_allmethods_plot, aes(x = as.numeric(net), y = auprc, colour = order)) + 
  geom_point() + 
  geom_line() + 
  scale_color_manual(values = c("dark green", "magenta", "blue")) +
  scale_x_log10(breaks = c(100, 500, 1000, 1500, 2000, 2500, 3000, 5000)) +
  ylab("Precision (AUPRC)") +
  xlab("Cutoff (log"[10]~"scale)") +
  coord_cartesian(ylim = c(0.55, 0.75)) +
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

dev.off()

pdf("graphics/ChIPorth_twonetlines_HOTexcl_AUPRC.pdf",
    width = 2.5,
    height =2.5)

ggplot(data = three_HOTexcl_benchraptor_allmethods_plot[!str_detect(three_HOTexcl_benchraptor_allmethods_plot$order, "CisBP"), ], aes(x = as.numeric(net), y = auprc, colour = order)) + 
  geom_point() + 
  geom_line() + 
  scale_color_manual(values = c("dark green", "magenta")) +
  scale_x_log10(breaks = c(100, 500, 1000, 1500, 2000, 2500, 3000, 5000)) +
  ylab("Precision (AUPRC)") +
  xlab("Cutoff (log"[10]~"scale)") +
  coord_cartesian(ylim = c(0.55, 0.75)) +
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

pdf("graphics/ChIPorth_twonetlines_HOTexcl_AUROC.pdf",
    width = 2.5,
    height =2.5)

ggplot(data = three_HOTexcl_benchraptor_allmethods_plot[!str_detect(three_HOTexcl_benchraptor_allmethods_plot$order, "CisBP"), ], aes(x = as.numeric(net), y = auroc, colour = order)) + 
  geom_point() + 
  geom_line() + 
  scale_color_manual(values = c("dark green", "magenta")) +
  scale_x_log10(breaks = c(100, 500, 1000, 1500, 2000, 2500, 3000, 5000)) +
  ylab("Sensitivity (AUROC)") +
  xlab("Cutoff (log"[10]~"scale)") +
  coord_cartesian(ylim = c(0.55, 0.75)) +
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

#### without HOT region filtering ####

ChIPorth_HOTincl_benchraptor_files <- bench_out_filelist[str_detect(bench_out_filelist, "benchRAPToR") &
                                                           str_detect(bench_out_filelist, 'HOTincl') &
                                                           str_detect(bench_out_filelist, "ChIPorth")]

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

## process networks with CisBP motifs instead of STREME motifs where applicable

STREMECisBP_HOTincl_benchraptor_files <- bench_out_filelist[str_detect(bench_out_filelist, "benchRAPToR") &
                                                              str_detect(bench_out_filelist, "HOTincl") &
                                                              str_detect(bench_out_filelist, "STRCisBP")]

STREMECisBP_HOTincl_benchraptor_files <- STREMECisBP_HOTincl_benchraptor_files[!str_detect(STREMECisBP_HOTincl_benchraptor_files, "dense")]
STREMECisBP_HOTincl_benchraptor_files <- STREMECisBP_HOTincl_benchraptor_files[!str_detect(STREMECisBP_HOTincl_benchraptor_files, "shufflestat")]

STREMECisBP_HOTincl_benchraptor_list <- lapply(STREMECisBP_HOTincl_benchraptor_files, function(x){
  
  read.table(paste0("output/benchmark_out/", x),
             sep = "\t",
             header = TRUE)
  
})

STREMECisBP_HOTincl_benchraptor_all <- do.call(rbind, STREMECisBP_HOTincl_benchraptor_list)

STREMECisBP_HOTincl_benchraptor_allmethods_plot <- STREMECisBP_HOTincl_benchraptor_all[STREMECisBP_HOTincl_benchraptor_all$metric %in% c("auroc", "auprc"), c("net", "method", "metric", "score")]
STREMECisBP_HOTincl_benchraptor_allmethods_plot <- tidyr::pivot_wider(STREMECisBP_HOTincl_benchraptor_allmethods_plot, names_from = c("metric"), values_from = "score")

STREMECisBP_HOTincl_benchraptor_allmethods_plot[, "label"] <- as.character(STREMECisBP_HOTincl_benchraptor_allmethods_plot$net)
# STREMECisBP_HOTincl_benchraptor_allmethods_plot[str_detect(unlist(STREMECisBP_HOTincl_benchraptor_allmethods_plot$net), "shuffle"), "label"] <- ""

STREMECisBP_HOTincl_benchraptor_allmethods_plot <- STREMECisBP_HOTincl_benchraptor_allmethods_plot[STREMECisBP_HOTincl_benchraptor_allmethods_plot$method == "mlm_estimate", ]

## add network with signal-based order for comparison

allChIP_HOTincl_benchraptor_allmethods_plot <- allChIP_HOTincl_benchraptor_allmethods_plot[allChIP_HOTincl_benchraptor_allmethods_plot$method == "mlm_estimate", ]

allChIP_HOTincl_benchraptor_allmethods_plot <- allChIP_HOTincl_benchraptor_allmethods_plot[, colnames(allChIP_HOTincl_benchraptor_allmethods_plot) %in% colnames(STREMECisBP_HOTincl_benchraptor_allmethods_plot)]

ChIPorth_HOTincl_benchraptor_allmethods_plot[, "order"] <- "STREME motif conservation"
allChIP_HOTincl_benchraptor_allmethods_plot[, "order"] <- "ChIP signal"
STREMECisBP_HOTincl_benchraptor_allmethods_plot[, "order"] <- "STREME+CisBP motif conservation"

three_HOTincl_benchraptor_allmethods_plot <- rbind(allChIP_HOTincl_benchraptor_allmethods_plot,
                                                   STREMECisBP_HOTincl_benchraptor_allmethods_plot,
                                                   ChIPorth_HOTincl_benchraptor_allmethods_plot)

three_HOTincl_benchraptor_allmethods_plot$order <- factor(three_HOTincl_benchraptor_allmethods_plot$order, 
                                                          levels = c("ChIP signal", "STREME+CisBP motif conservation", "STREME motif conservation"))

saveRDS(three_HOTincl_benchraptor_allmethods_plot,
        "plotdata/three_HOTincl_benchraptor_allmethods_plot.rds")

pdf("graphics/ChIPorth_allnetlines_HOTincl_AUPRC.pdf",
    width = 2.5,
    height =2.5)

ggplot(data = three_HOTincl_benchraptor_allmethods_plot, aes(x = as.numeric(net), y = auprc, colour = order)) + 
  geom_point() + 
  geom_line() + 
  scale_color_manual(values = c("dark green", "blue", "magenta")) +
  scale_x_log10(breaks = c(100, 500, 1000, 1500, 2000, 2500, 3000, 5000)) +
  ylab("Precision (AUPRC)") +
  xlab("Cutoff (log"[10]~"scale)") +
  coord_cartesian(ylim = c(0.55, 0.75)) +
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

pdf("graphics/ChIPorth_allnetlines_HOTincl_AUROC.pdf",
    width = 2.5,
    height =2.5)

ggplot(data = three_HOTincl_benchraptor_allmethods_plot, aes(x = as.numeric(net), y = auroc, colour = order)) + 
  geom_point() + 
  geom_line() + 
  scale_color_manual(values = c("dark green", "blue", "magenta")) +
  scale_x_log10(breaks = c(100, 500, 1000, 1500, 2000, 2500, 3000, 5000)) +
  ylab("Sensitivity (AUROC)") +
  xlab("Cutoff (log"[10]~"scale)") +
  coord_cartesian(ylim = c(0.55, 0.75)) +
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

pdf("graphics/ChIPorth_allnetlines_HOTincl_AUPRC_LEGEND.pdf",
    width = 5,
    height = 5)

ggplot(data = three_HOTincl_benchraptor_allmethods_plot, aes(x = as.numeric(net), y = auprc, colour = order)) + 
  geom_point() + 
  geom_line() + 
  scale_color_manual(values = c("dark green", "blue", "magenta")) +
  scale_x_log10(breaks = c(100, 500, 1000, 1500, 2000, 2500, 3000, 5000)) +
  ylab("Precision (AUPRC)") +
  xlab("Cutoff (log"[10]~"scale)") +
  coord_cartesian(ylim = c(0.55, 0.75)) +
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

dev.off()

pdf("graphics/ChIPorth_twonetlines_HOTincl_AUPRC.pdf",
    width = 2.5,
    height =2.5)

ggplot(data = three_HOTincl_benchraptor_allmethods_plot[!str_detect(three_HOTincl_benchraptor_allmethods_plot$order, "CisBP"), ], aes(x = as.numeric(net), y = auprc, colour = order)) + 
  geom_point() + 
  geom_line() + 
  scale_color_manual(values = c("dark green", "magenta")) +
  scale_x_log10(breaks = c(100, 500, 1000, 1500, 2000, 2500, 3000, 5000)) +
  ylab("Precision (AUPRC)") +
  xlab("Cutoff (log"[10]~"scale)") +
  coord_cartesian(ylim = c(0.55, 0.75)) +
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

pdf("graphics/ChIPorth_twonetlines_HOTincl_AUROC.pdf",
    width = 2.5,
    height =2.5)

ggplot(data = three_HOTincl_benchraptor_allmethods_plot[!str_detect(three_HOTincl_benchraptor_allmethods_plot$order, "CisBP"), ], aes(x = as.numeric(net), y = auroc, colour = order)) + 
  geom_point() + 
  geom_line() + 
  scale_color_manual(values = c("dark green", "magenta")) +
  scale_x_log10(breaks = c(100, 500, 1000, 1500, 2000, 2500, 3000, 5000)) +
  ylab("Sensitivity (AUROC)") +
  xlab("Cutoff (log"[10]~"scale)") +
  coord_cartesian(ylim = c(0.55, 0.75)) +
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

