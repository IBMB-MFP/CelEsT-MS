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

allChIP_HOTincl_benchraw_files <- bench_out_filelist[str_detect(bench_out_filelist, "benchRAW") & str_detect(bench_out_filelist, "HOTincl")]

allChIP_HOTincl_benchraw_files_dense <- allChIP_HOTincl_benchraw_files[str_detect(allChIP_HOTincl_benchraw_files, "dense")]
allChIP_HOTincl_benchraw_files <- allChIP_HOTincl_benchraw_files[!str_detect(allChIP_HOTincl_benchraw_files, "dense")]

allChIP_HOTincl_benchraw_list <- lapply(allChIP_HOTincl_benchraw_files, function(x){
  
  read.table(paste0("output/benchmark_out/", x),
             sep = "\t",
             header = TRUE)
  
})

allChIP_HOTincl_benchraw_all <- do.call(rbind, allChIP_HOTincl_benchraw_list)

# allChIP_HOTincl_benchraw_plot <- allChIP_HOTincl_benchraw_all[allChIP_HOTincl_benchraw_all$metric %in% c("auroc", "auprc") & allChIP_HOTincl_benchraw_all$method == "mlm_estimate", c("net", "metric", "score")]
# 
# allChIP_HOTincl_benchraw_plot <- tidyr::pivot_wider(allChIP_HOTincl_benchraw_plot, names_from = "metric", values_from = "score")
# allChIP_HOTincl_benchraw_plot[, "colour"] <- "red"
# allChIP_HOTincl_benchraw_plot[str_detect(allChIP_HOTincl_benchraw_plot$net, "shuffle"), "colour"] <- "grey"
# 
# ggplot(data = allChIP_HOTincl_benchraw_plot, aes(x = auroc, y = auprc, label = net)) + 
#   geom_point(col = allChIP_HOTincl_benchraw_plot$colour) + 
#   geom_label_repel() + 
#   theme_classic() + 
#   geom_vline(xintercept = 0.5,
#              linetype = "dashed",
#              col = "grey") +
#   geom_hline(yintercept = 0.5,
#              linetype = "dashed",
#              col = "grey")

allChIP_HOTincl_benchraw_allmethods_plot <- allChIP_HOTincl_benchraw_all[allChIP_HOTincl_benchraw_all$metric %in% c("auroc", "auprc"), c("net", "method", "metric", "score")]
allChIP_HOTincl_benchraw_allmethods_plot <- tidyr::pivot_wider(allChIP_HOTincl_benchraw_allmethods_plot, names_from = c("metric"), values_from = "score")

allChIP_HOTincl_benchraw_allmethods_plot[, "colour"] <- "red"
allChIP_HOTincl_benchraw_allmethods_plot[str_detect(allChIP_HOTincl_benchraw_allmethods_plot$net, "shuffle"), "colour"] <- "grey"

allChIP_HOTincl_benchraw_allmethods_plot[, "label"] <- allChIP_HOTincl_benchraw_allmethods_plot$net
allChIP_HOTincl_benchraw_allmethods_plot[str_detect(unlist(allChIP_HOTincl_benchraw_allmethods_plot$net), "shuffle"), "label"] <- ""

allChIP_HOTincl_benchraw_allmethods_plot[, "cutoff"] <- str_remove(allChIP_HOTincl_benchraw_allmethods_plot$net, "_shuffle")
allChIP_HOTincl_benchraw_allmethods_plot[, "shuffled"] <- str_detect(allChIP_HOTincl_benchraw_allmethods_plot$net, "_shuffle")

ggplot(data = allChIP_HOTincl_benchraw_allmethods_plot, aes(x = auroc, y = auprc, label = label, alpha = shuffled, color = cutoff)) + 
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

ggplot(data = allChIP_HOTincl_benchraw_allmethods_plot[allChIP_HOTincl_benchraw_allmethods_plot$method == "mlm_estimate" & allChIP_HOTincl_benchraw_allmethods_plot$cutoff %in% c("100", "500", "1000", "2000", "3000", "5000"), ], aes(x = auroc, y = auprc, label = label, alpha = shuffled, color = cutoff)) + 
  geom_point() + 
  geom_label_repel(label.size = 0.1,
                   max.overlaps = 15,
                   box.padding = 0.1,
                   label.padding = 0.1,
                   point.padding = 1,
                   color = "black") + 
  theme_light() + 
  scale_alpha_discrete(range = c(1, 0.4)) +
  geom_vline(xintercept = 0.5,
             linetype = "dashed",
             col = "grey") +
  geom_hline(yintercept = 0.5,
             linetype = "dashed",
             col = "grey") + 
  coord_cartesian(xlim = c(0.4, 0.8),
                  ylim = c(0.4, 0.8))


##### DENSE #####


allChIP_HOTincl_benchraw_dense_list <- lapply(allChIP_HOTincl_benchraw_files_dense, function(x){
  
  read.table(paste0("output/benchmark_out/", x),
             sep = "\t",
             header = TRUE)
  
})

allChIP_HOTincl_benchraw_dense_all <- do.call(rbind, allChIP_HOTincl_benchraw_dense_list)

# allChIP_HOTincl_benchraw_dense_plot <- allChIP_HOTincl_benchraw_dense_all[allChIP_HOTincl_benchraw_dense_all$metric %in% c("auroc", "auprc") & allChIP_HOTincl_benchraw_dense_all$method == "mlm_estimate", c("net", "metric", "score")]
# 
# allChIP_HOTincl_benchraw_dense_plot <- tidyr::pivot_wider(allChIP_HOTincl_benchraw_dense_plot, names_from = "metric", values_from = "score")
# allChIP_HOTincl_benchraw_dense_plot[, "colour"] <- "red"
# allChIP_HOTincl_benchraw_dense_plot[str_detect(allChIP_HOTincl_benchraw_dense_plot$net, "shuffle"), "colour"] <- "grey"
# 
# ggplot(data = allChIP_HOTincl_benchraw_dense_plot, aes(x = auroc, y = auprc, label = net)) + 
#   geom_point(col = allChIP_HOTincl_benchraw_dense_plot$colour) + 
#   geom_label_repel() + 
#   theme_classic() + 
#   geom_vline(xintercept = 0.5,
#              linetype = "dashed",
#              col = "grey") +
#   geom_hline(yintercept = 0.5,
#              linetype = "dashed",
#              col = "grey")

allChIP_HOTincl_benchraw_dense_allmethods_plot <- allChIP_HOTincl_benchraw_dense_all[allChIP_HOTincl_benchraw_dense_all$metric %in% c("auroc", "auprc"), c("net", "method", "metric", "score")]
allChIP_HOTincl_benchraw_dense_allmethods_plot <- tidyr::pivot_wider(allChIP_HOTincl_benchraw_dense_allmethods_plot, names_from = c("metric"), values_from = "score")

allChIP_HOTincl_benchraw_dense_allmethods_plot[, "colour"] <- "red"
allChIP_HOTincl_benchraw_dense_allmethods_plot[str_detect(allChIP_HOTincl_benchraw_dense_allmethods_plot$net, "shuffle"), "colour"] <- "grey"

allChIP_HOTincl_benchraw_dense_allmethods_plot[, "label"] <- allChIP_HOTincl_benchraw_dense_allmethods_plot$net
allChIP_HOTincl_benchraw_dense_allmethods_plot[str_detect(unlist(allChIP_HOTincl_benchraw_dense_allmethods_plot$net), "shuffle"), "label"] <- ""

allChIP_HOTincl_benchraw_dense_allmethods_plot[, "cutoff"] <- str_remove(allChIP_HOTincl_benchraw_dense_allmethods_plot$net, "_shuffle")
allChIP_HOTincl_benchraw_dense_allmethods_plot[, "shuffled"] <- str_detect(allChIP_HOTincl_benchraw_dense_allmethods_plot$net, "_shuffle")

ggplot(data = allChIP_HOTincl_benchraw_dense_allmethods_plot, aes(x = auroc, y = auprc, label = label, alpha = shuffled, color = cutoff)) + 
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

ggplot(data = allChIP_HOTincl_benchraw_dense_allmethods_plot[allChIP_HOTincl_benchraw_dense_allmethods_plot$method == "mlm_estimate" & allChIP_HOTincl_benchraw_dense_allmethods_plot$cutoff %in% c("100", "500", "1000", "2000", "3000", "5000"), ], aes(x = auroc, y = auprc, label = label, alpha = shuffled, color = cutoff)) + 
  geom_point() + 
  geom_label_repel(label.size = 0.1,
                   max.overlaps = 15,
                   box.padding = 0.1,
                   label.padding = 0.1,
                   point.padding = 1,
                   color = "black") + 
  theme_light() + 
  scale_alpha_discrete(range = c(1, 0.4)) +
  geom_vline(xintercept = 0.5,
             linetype = "dashed",
             col = "grey") +
  geom_hline(yintercept = 0.5,
             linetype = "dashed",
             col = "grey") + 
  coord_cartesian(xlim = c(0.4, 0.8),
                  ylim = c(0.4, 0.8))

# 
# ggplot(data = allChIP_HOTincl_benchraw_allmethods_plot[!str_detect(allChIP_HOTincl_benchraw_allmethods_plot$net, "shuffle"), ], aes(x = as.numeric(net), y = auroc)) + 
#   geom_point() + 
#   facet_wrap(~ method) + 
#   scale_x_log10()
#   

allChIP_HOTexcl_benchraw_files <- bench_out_filelist[str_detect(bench_out_filelist, "benchRAW") & str_detect(bench_out_filelist, "HOTexcl")]
allChIP_HOTexcl_benchraw_files <- allChIP_HOTexcl_benchraw_files[!str_detect(allChIP_HOTexcl_benchraw_files, "absolute")]
allChIP_HOTexcl_benchraw_files <- allChIP_HOTexcl_benchraw_files[!str_detect(allChIP_HOTexcl_benchraw_files, "absoliute")]


allChIP_HOTexcl_benchraw_files_dense <- allChIP_HOTexcl_benchraw_files[str_detect(allChIP_HOTexcl_benchraw_files, "dense")]
allChIP_HOTexcl_benchraw_files <- allChIP_HOTexcl_benchraw_files[!str_detect(allChIP_HOTexcl_benchraw_files, "dense")]

allChIP_HOTexcl_benchraw_list <- lapply(allChIP_HOTexcl_benchraw_files, function(x){
  
  read.table(paste0("output/benchmark_out/", x),
             sep = "\t",
             header = TRUE)
  
})

allChIP_HOTexcl_benchraw_all <- do.call(rbind, allChIP_HOTexcl_benchraw_list)

# allChIP_HOTexcl_benchraw_plot <- allChIP_HOTexcl_benchraw_all[allChIP_HOTexcl_benchraw_all$metric %in% c("auroc", "auprc") & allChIP_HOTexcl_benchraw_all$method == "mlm_estimate", c("net", "metric", "score")]
# 
# allChIP_HOTexcl_benchraw_plot <- tidyr::pivot_wider(allChIP_HOTexcl_benchraw_plot, names_from = "metric", values_from = "score")
# allChIP_HOTexcl_benchraw_plot[, "colour"] <- "red"
# allChIP_HOTexcl_benchraw_plot[str_detect(allChIP_HOTexcl_benchraw_plot$net, "shuffle"), "colour"] <- "grey"
# 
# ggplot(data = allChIP_HOTexcl_benchraw_plot, aes(x = auroc, y = auprc, label = net)) + 
#   geom_point(col = allChIP_HOTexcl_benchraw_plot$colour) + 
#   geom_label_repel() + 
#   theme_classic() + 
#   geom_vline(xintercept = 0.5,
#              linetype = "dashed",
#              col = "grey") +
#   geom_hline(yintercept = 0.5,
#              linetype = "dashed",
#              col = "grey")

allChIP_HOTexcl_benchraw_allmethods_plot <- allChIP_HOTexcl_benchraw_all[allChIP_HOTexcl_benchraw_all$metric %in% c("auroc", "auprc"), c("net", "method", "metric", "score")]
allChIP_HOTexcl_benchraw_allmethods_plot <- tidyr::pivot_wider(allChIP_HOTexcl_benchraw_allmethods_plot, names_from = "metric", values_from = "score")

allChIP_HOTexcl_benchraw_allmethods_plot[, "colour"] <- "red"
allChIP_HOTexcl_benchraw_allmethods_plot[str_detect(allChIP_HOTexcl_benchraw_allmethods_plot$net, "shuffle"), "colour"] <- "grey"

allChIP_HOTexcl_benchraw_allmethods_plot[, "label"] <- allChIP_HOTexcl_benchraw_allmethods_plot$net
allChIP_HOTexcl_benchraw_allmethods_plot[str_detect(unlist(allChIP_HOTexcl_benchraw_allmethods_plot$net), "shuffle"), "label"] <- ""
allChIP_HOTexcl_benchraw_allmethods_plot[, "cutoff"] <- str_remove(allChIP_HOTexcl_benchraw_allmethods_plot$net, "_shuffle")
allChIP_HOTexcl_benchraw_allmethods_plot[, "shuffled"] <- str_detect(allChIP_HOTexcl_benchraw_allmethods_plot$net, "_shuffle")

ggplot(data = allChIP_HOTexcl_benchraw_allmethods_plot, aes(x = auroc, y = auprc, label = label, alpha = shuffled, color = cutoff)) + 
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

ggplot(data = allChIP_HOTexcl_benchraw_allmethods_plot[allChIP_HOTexcl_benchraw_allmethods_plot$method == "mlm_estimate" & allChIP_HOTexcl_benchraw_allmethods_plot$cutoff %in% c("100", "500", "1000", "2000", "3000", "5000"), ], aes(x = auroc, y = auprc, label = label, alpha = shuffled, color = cutoff)) + 
  geom_point() + 
  geom_label_repel(label.size = 0.1,
                   max.overlaps = 15,
                   box.padding = 0.1,
                   label.padding = 0.1,
                   point.padding = 1,
                   color = "black") + 
  theme_light() + 
  scale_alpha_discrete(range = c(1, 0.4)) +
  geom_vline(xintercept = 0.5,
             linetype = "dashed",
             col = "grey") +
  geom_hline(yintercept = 0.5,
             linetype = "dashed",
             col = "grey") + 
  coord_cartesian(xlim = c(0.4, 0.8),
                  ylim = c(0.4, 0.8))

#### DENSE

allChIP_HOTexcl_benchraw_dense_list <- lapply(allChIP_HOTexcl_benchraw_files_dense, function(x){
  
  read.table(paste0("output/benchmark_out/", x),
             sep = "\t",
             header = TRUE)
  
})

allChIP_HOTexcl_benchraw_dense_all <- do.call(rbind, allChIP_HOTexcl_benchraw_dense_list)

# allChIP_HOTexcl_benchraw_dense_plot <- allChIP_HOTexcl_benchraw_dense_all[allChIP_HOTexcl_benchraw_dense_all$metric %in% c("auroc", "auprc") & allChIP_HOTexcl_benchraw_dense_all$method == "mlm_estimate", c("net", "metric", "score")]
# 
# allChIP_HOTexcl_benchraw_dense_plot <- tidyr::pivot_wider(allChIP_HOTexcl_benchraw_dense_plot, names_from = "metric", values_from = "score")
# allChIP_HOTexcl_benchraw_dense_plot[, "colour"] <- "red"
# allChIP_HOTexcl_benchraw_dense_plot[str_detect(allChIP_HOTexcl_benchraw_dense_plot$net, "shuffle"), "colour"] <- "grey"
# 
# ggplot(data = allChIP_HOTexcl_benchraw_dense_plot, aes(x = auroc, y = auprc, label = net)) + 
#   geom_point(col = allChIP_HOTexcl_benchraw_dense_plot$colour) + 
#   geom_label_repel() + 
#   theme_classic() + 
#   geom_vline(xintercept = 0.5,
#              linetype = "dashed",
#              col = "grey") +
#   geom_hline(yintercept = 0.5,
#              linetype = "dashed",
#              col = "grey")

allChIP_HOTexcl_benchraw_dense_allmethods_plot <- allChIP_HOTexcl_benchraw_dense_all[allChIP_HOTexcl_benchraw_dense_all$metric %in% c("auroc", "auprc"), c("net", "method", "metric", "score")]
allChIP_HOTexcl_benchraw_dense_allmethods_plot <- tidyr::pivot_wider(allChIP_HOTexcl_benchraw_dense_allmethods_plot, names_from = c("metric"), values_from = "score")

allChIP_HOTexcl_benchraw_dense_allmethods_plot[, "colour"] <- "red"
allChIP_HOTexcl_benchraw_dense_allmethods_plot[str_detect(allChIP_HOTexcl_benchraw_dense_allmethods_plot$net, "shuffle"), "colour"] <- "grey"

allChIP_HOTexcl_benchraw_dense_allmethods_plot[, "label"] <- allChIP_HOTexcl_benchraw_dense_allmethods_plot$net
allChIP_HOTexcl_benchraw_dense_allmethods_plot[str_detect(unlist(allChIP_HOTexcl_benchraw_dense_allmethods_plot$net), "shuffle"), "label"] <- ""

allChIP_HOTexcl_benchraw_dense_allmethods_plot[, "cutoff"] <- str_remove(allChIP_HOTexcl_benchraw_dense_allmethods_plot$net, "_shuffle")
allChIP_HOTexcl_benchraw_dense_allmethods_plot[, "shuffled"] <- str_detect(allChIP_HOTexcl_benchraw_dense_allmethods_plot$net, "_shuffle")

ggplot(data = allChIP_HOTexcl_benchraw_dense_allmethods_plot, aes(x = auroc, y = auprc, label = label, alpha = shuffled, color = cutoff)) + 
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

ggplot(data = allChIP_HOTexcl_benchraw_dense_allmethods_plot[allChIP_HOTexcl_benchraw_dense_allmethods_plot$method == "mlm_estimate" & allChIP_HOTexcl_benchraw_dense_allmethods_plot$cutoff %in% c("100", "500", "1000", "2000", "3000", "5000"), ], aes(x = auroc, y = auprc, label = label, alpha = shuffled, color = cutoff)) + 
  geom_point() + 
  geom_label_repel(label.size = 0.1,
                   max.overlaps = 15,
                   box.padding = 0.1,
                   label.padding = 0.1,
                   point.padding = 1,
                   color = "black") + 
  theme_light() + 
  scale_alpha_discrete(range = c(1, 0.4)) +
  geom_vline(xintercept = 0.5,
             linetype = "dashed",
             col = "grey") +
  geom_hline(yintercept = 0.5,
             linetype = "dashed",
             col = "grey") + 
  coord_cartesian(xlim = c(0.4, 0.8),
                  ylim = c(0.4, 0.8))



ggplot(data = allChIP_HOTexcl_benchraw_allmethods_plot[!str_detect(allChIP_HOTexcl_benchraw_allmethods_plot$net, "shuffle"), ], aes(x = as.numeric(net), y = auroc)) +
  geom_point() +
  facet_wrap(~ method) +
  scale_x_log10()


ggplot(data = allChIP_HOTexcl_benchraw_allmethods_plot[!str_detect(allChIP_HOTexcl_benchraw_allmethods_plot$net, "shuffle"), ], aes(x = as.numeric(net), y = auprc)) +
  geom_point() +
  facet_wrap(~ method) +
  scale_x_log10()






allChIP_HOTexcl_benchraptor_files <- bench_out_filelist[str_detect(bench_out_filelist, "benchRAPToR") & str_detect(bench_out_filelist, "HOTexcl")]
allChIP_HOTexcl_benchraptor_files <- allChIP_HOTexcl_benchraptor_files[!str_detect(allChIP_HOTexcl_benchraptor_files, "absolute")]
allChIP_HOTexcl_benchraptor_files <- allChIP_HOTexcl_benchraptor_files[!str_detect(allChIP_HOTexcl_benchraptor_files, "absoliute")]


allChIP_HOTexcl_benchraptor_files_dense <- allChIP_HOTexcl_benchraptor_files[str_detect(allChIP_HOTexcl_benchraptor_files, "dense")]
allChIP_HOTexcl_benchraptor_files <- allChIP_HOTexcl_benchraptor_files[!str_detect(allChIP_HOTexcl_benchraptor_files, "dense")]

allChIP_HOTexcl_benchraptor_list <- lapply(allChIP_HOTexcl_benchraptor_files, function(x){
  
  read.table(paste0("output/benchmark_out/", x),
             sep = "\t",
             header = TRUE)
  
})

allChIP_HOTexcl_benchraptor_all <- do.call(rbind, allChIP_HOTexcl_benchraptor_list)

# allChIP_HOTexcl_benchraptor_plot <- allChIP_HOTexcl_benchraptor_all[allChIP_HOTexcl_benchraptor_all$metric %in% c("auroc", "auprc") & allChIP_HOTexcl_benchraptor_all$method == "mlm_estimate", c("net", "metric", "score")]
# 
# allChIP_HOTexcl_benchraptor_plot <- tidyr::pivot_wider(allChIP_HOTexcl_benchraptor_plot, names_from = "metric", values_from = "score")
# allChIP_HOTexcl_benchraptor_plot[, "colour"] <- "red"
# allChIP_HOTexcl_benchraptor_plot[str_detect(allChIP_HOTexcl_benchraptor_plot$net, "shuffle"), "colour"] <- "grey"
# 
# ggplot(data = allChIP_HOTexcl_benchraptor_plot, aes(x = auroc, y = auprc, label = net)) + 
#   geom_point(col = allChIP_HOTexcl_benchraptor_plot$colour) + 
#   geom_label_repel() + 
#   theme_classic() + 
#   geom_vline(xintercept = 0.5,
#              linetype = "dashed",
#              col = "grey") +
#   geom_hline(yintercept = 0.5,
#              linetype = "dashed",
#              col = "grey")

allChIP_HOTexcl_benchraptor_allmethods_plot <- allChIP_HOTexcl_benchraptor_all[allChIP_HOTexcl_benchraptor_all$metric %in% c("auroc", "auprc"), c("net", "method", "metric", "score")]
allChIP_HOTexcl_benchraptor_allmethods_plot <- tidyr::pivot_wider(allChIP_HOTexcl_benchraptor_allmethods_plot, names_from = "metric", values_from = "score")

allChIP_HOTexcl_benchraptor_allmethods_plot[, "colour"] <- "red"
allChIP_HOTexcl_benchraptor_allmethods_plot[str_detect(allChIP_HOTexcl_benchraptor_allmethods_plot$net, "shuffle"), "colour"] <- "grey"

allChIP_HOTexcl_benchraptor_allmethods_plot[, "label"] <- allChIP_HOTexcl_benchraptor_allmethods_plot$net
allChIP_HOTexcl_benchraptor_allmethods_plot[str_detect(unlist(allChIP_HOTexcl_benchraptor_allmethods_plot$net), "shuffle"), "label"] <- ""
allChIP_HOTexcl_benchraptor_allmethods_plot[, "cutoff"] <- str_remove(allChIP_HOTexcl_benchraptor_allmethods_plot$net, "_shuffle")
allChIP_HOTexcl_benchraptor_allmethods_plot[, "shuffled"] <- str_detect(allChIP_HOTexcl_benchraptor_allmethods_plot$net, "_shuffle")

ggplot(data = allChIP_HOTexcl_benchraptor_allmethods_plot, aes(x = auroc, y = auprc, label = label, alpha = shuffled, color = cutoff)) + 
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

ggplot(data = allChIP_HOTexcl_benchraptor_allmethods_plot[allChIP_HOTexcl_benchraptor_allmethods_plot$method == "mlm_estimate" & allChIP_HOTexcl_benchraptor_allmethods_plot$cutoff %in% c("100", "500", "1000", "2000", "3000", "5000"), ], aes(x = auroc, y = auprc, label = label, alpha = shuffled, color = cutoff)) + 
  geom_point() + 
  geom_label_repel(label.size = 0.1,
                   max.overlaps = 15,
                   box.padding = 0.1,
                   label.padding = 0.1,
                   point.padding = 1,
                   color = "black") + 
  theme_light() + 
  scale_alpha_discrete(range = c(1, 0.4)) +
  geom_vline(xintercept = 0.5,
             linetype = "dashed",
             col = "grey") +
  geom_hline(yintercept = 0.5,
             linetype = "dashed",
             col = "grey") + 
  coord_cartesian(xlim = c(0.4, 0.8),
                  ylim = c(0.4, 0.8))


#### try with shuffle stats

library(dplyr)
shufflestats_files <- bench_out_filelist[str_detect(bench_out_filelist, "_shufflestats")]

benchRAPToR_HOTincl_shufflestats_files <- shufflestats_files[str_detect(shufflestats_files, "RAPToR") & str_detect(shufflestats_files, "HOTincl")]

benchRAPToR_HOTincl_shufflestats_lists <- lapply(benchRAPToR_HOTincl_shufflestats_files, function(file){

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

benchRAPToR_HOTincl_shufflestats_df <- do.call(rbind, benchRAPToR_HOTincl_shufflestats_lists)
benchRAPToR_HOTincl_shufflestats_df[, "cutoff"] <- str_remove(benchRAPToR_HOTincl_shufflestats_df$net, "_shuffle")

allChIP_HOTincl_benchraptor_allmethods_plot[, "cutoff"] <- str_remove(allChIP_HOTincl_benchraptor_allmethods_plot$net, "_shuffle")
allChIP_HOTincl_benchraptor_allmethods_plot[, "shuffled"] <- str_detect(allChIP_HOTincl_benchraptor_allmethods_plot$net, "_shuffle")

ggplot(data = allChIP_HOTincl_benchraptor_allmethods_plot[!str_detect(allChIP_HOTincl_benchraptor_allmethods_plot$net, "shuffle"), ], aes(x = auroc, y = auprc, colour = cutoff)) + 
  geom_point() +
  geom_label_repel(aes(label = label),
                   label.size = 0.1,
                   max.overlaps = 50,
                   box.padding = 0.1,
                   label.padding = 0.1,
                   point.padding = 1,
                   size = 2) +
  theme_classic() + 
  geom_vline(xintercept = 0.5,
             linetype = "dashed",
             col = "grey") +
  geom_hline(yintercept = 0.5,
             linetype = "dashed",
             col = "grey") + 
  facet_wrap(~ method) + 
  geom_point(data = benchRAPToR_HOTincl_shufflestats_df, aes(x = mean_auroc, y = mean_auprc), alpha = 0.5) + 
  geom_errorbar(data = benchRAPToR_HOTincl_shufflestats_df, aes(x = mean_auroc, y = mean_auprc,
                                                                ymin = (mean_auprc - sd_auprc),
                                                                ymax = (mean_auprc + sd_auprc)), alpha = 0.5) +
  geom_errorbar(data = benchRAPToR_HOTincl_shufflestats_df, aes(x = mean_auroc, y = mean_auprc,
                                                                xmin = (mean_auroc - sd_auroc),
                                                                xmax = (mean_auroc + sd_auroc)), alpha = 0.5)


