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

packages <- c("ggplot2")

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

allmodERN_TFsONLY_BM <- readRDS("output/allmodERN_TFsONLY_BM.rds")

allTFS_manualtoptargets <- readRDS("output/MODernENCODE_manualtoptargets_operonexcluded.rds")
allTFS_manualtoptargets_HOTexcl <- readRDS("output/MODernENCODE_manualtoptargets_operon&HOTexcluded.rds")

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

#### PLOT HOTCUT ####

HOTcut_TF_numbers <- sapply(hotcutfiles, function(thisfile){
  
  Hotcut <- str_extract(thisfile, "HOTcut[0-9]{2}")
  
  thisGRN <- read.table(paste0("output/GRNs/allChIP_3000_", Hotcut, ".txt"),
                        header = TRUE)
  
  length(unique(thisGRN$source))
  
})

names(HOTcut_TF_numbers) <- str_remove(str_remove(names(HOTcut_TF_numbers), "MODernENCODE_manualtoptargets_operon&HOTexcluded_HOTcut"), "\\.rds")

#### INPUT DATA ####

method_labels <- c("Consensus", "MLM", "ULM", "WSum")
names(method_labels) <- c("consensus_estimate", "mlm_estimate", "ulm_estimate", "wsum_norm")

bench_out_filelist <- list.files("output/benchmark_out/")

HOTcut_benchraptor_files <- bench_out_filelist[str_detect(bench_out_filelist, "HotCut")]

HOTcut_benchraptor_list <- lapply(HOTcut_benchraptor_files, function(x){
  
  read.table(paste0("output/benchmark_out/", x),
             sep = "\t",
             header = TRUE)
  
})

HOTcut_benchraptor_all <- do.call(rbind, HOTcut_benchraptor_list)

HOTcut_benchraptor_allmethods_plot <- HOTcut_benchraptor_all[HOTcut_benchraptor_all$metric %in% c("auroc", "auprc"), c("net", "method", "metric", "score")]
HOTcut_benchraptor_allmethods_plot <- tidyr::pivot_wider(HOTcut_benchraptor_allmethods_plot, names_from = c("metric"), values_from = "score")

HOTcut_benchraptor_allmethods_plot[, "label"] <- HOTcut_benchraptor_allmethods_plot$net


pdf("graphics/allChIP_HOTcut__benchraptor_cutoffs_methodlines_AUPRC.pdf",
    height = 2,
    width = 2)

ggplot(data = HOTcut_benchraptor_allmethods_plot[!HOTcut_benchraptor_allmethods_plot$method %in% c("wsum_estimate", "wsum_corr"), ], aes(x = as.numeric(cutoff), y = auprc, colour = method)) + 
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
        axis.title = element_text(size = 10) ) + 
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

pdf("graphics/allChIP_HOTcut_methodlines_AUROC.pdf",
    height = 2,
    width = 2)

ggplot(data = HOTcut_benchraptor_allmethods_plot[!HOTcut_benchraptor_allmethods_plot$method %in% c("wsum_estimate", "wsum_corr"), ], aes(x = as.numeric(cutoff), y = auroc, colour = method)) + 
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

