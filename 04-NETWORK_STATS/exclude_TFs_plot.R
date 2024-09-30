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
              "stringr")

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

benchmark_obs <- read.table("output/benchmark_observations.txt",
                            sep = "\t",
                            header = TRUE)

exclude_files <- list.files("output/benchmark_out/")

CelEsT_bench_RAW <- read.table("output/benchmark_out/combo3_benchRAW.tsv",
                               header = TRUE,
                               sep = "\t")

CelEsT_bench <- read.table("output/benchmark_out/combo3_benchRAPToR.tsv",
                           header = TRUE,
                           sep = "\t")

#### Benchmark for excluded TFs RAPToR-correction ####

exclude_files_RAPToR <- exclude_files[str_detect(exclude_files, "EXCLUDEsample") & !str_detect(exclude_files, "RAW")]

exclude_stats <- sapply(exclude_files_RAPToR, function(thisfile){

  tempin <- read.table(paste0("output/benchmark_out/", thisfile),
                       header = TRUE,
                       sep = "\t")
  
  c("auprc" = tempin[tempin$metric == "auprc", "score"],
    "auroc" = tempin[tempin$metric == "auroc", "score"])
  
})

CelEsT_bench_df <- data.frame("auprc" = CelEsT_bench[CelEsT_bench$method == "mlm_estimate" & CelEsT_bench$metric == "auprc" & CelEsT_bench$net == "combo3equal", "score"],
  "auroc" = CelEsT_bench[CelEsT_bench$method == "mlm_estimate" & CelEsT_bench$metric == "auroc" & CelEsT_bench$net == "combo3equal", "score"])
 
exclude_stats_df <- data.frame(t(exclude_stats))

exclude_benchraptor_shufflestats <- lapply("EXCLUDEsample_benchRAPToR_shufflestats.tsv", function(file){
  
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
  
  tempout
  
})

exclude_benchraptor_shufflestats_df <- as.data.frame(t(unlist(exclude_benchraptor_shufflestats)))
row.names(exclude_benchraptor_shufflestats_df) <- exclude_benchraptor_shufflestats_df[, 1]
exclude_benchraptor_shufflestats_df <- exclude_benchraptor_shufflestats_df[, 2:ncol(exclude_benchraptor_shufflestats_df)]
exclude_benchraptor_shufflestats_df[] <- lapply(exclude_benchraptor_shufflestats_df, as.numeric)

pdf("graphics/CelEsT_exclude.pdf",
    height = 2.5,
    width = 2.5)

ggplot(exclude_stats_df, aes(x = auroc, y = auprc)) + 
  geom_point(alpha = 0.5,
             colour = "dark grey",
             size = 0.5) +
  theme_classic() +
  theme(strip.text = element_text(colour = 'white'),
        strip.background = element_rect(fill = "black"),
        legend.position = "none",
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_line(colour = "grey", linewidth = 0.2),
        panel.grid.minor = element_line(colour = "grey", linewidth = 0.1),
        axis.text = element_text(colour = "black")) +
  coord_cartesian(xlim = c(0.45, 0.8),
                  ylim = c(0.45, 0.8)) + 
  geom_point(data = CelEsT_bench_df,
             colour = "red",
             size = 2) +
  xlab("") + 
  ylab("") +
  geom_vline(xintercept = 0.5,
             linetype = "dashed",
             col = "black") +
  geom_hline(yintercept = 0.5,
             linetype = "dashed",
             col = "black") + 
  geom_errorbar(data = exclude_benchraptor_shufflestats_df, aes(x = mean_auroc, y = mean_auprc,
                                                                ymin = (mean_auprc - sd_auprc),
                                                                ymax = (mean_auprc + sd_auprc)),
                alpha = 0.5, width = 0.01) +
  geom_errorbar(data = exclude_benchraptor_shufflestats_df, aes(x = mean_auroc, y = mean_auprc,
                                                                xmin = (mean_auroc - sd_auroc),
                                                                xmax = (mean_auroc + sd_auroc)), 
                alpha = 0.5, width = 0.01)


dev.off()

#### Benchmark for excluded TFs RAW ####

exclude_files_RAW <- exclude_files[str_detect(exclude_files, "EXCLUDEsample") & str_detect(exclude_files, "RAW")]

exclude_stats_RAW <- sapply(exclude_files_RAW, function(thisfile){
  
  tempin <- read.table(paste0("output/benchmark_out/", thisfile),
                       header = TRUE,
                       sep = "\t")
  
  c("auprc" = tempin[tempin$metric == "auprc", "score"],
    "auroc" = tempin[tempin$metric == "auroc", "score"])
  
})

CelEsT_bench_RAW_df <- data.frame("auprc" = CelEsT_bench_RAW[CelEsT_bench_RAW$method == "mlm_estimate" & CelEsT_bench_RAW$metric == "auprc" & CelEsT_bench_RAW$net == "combo3equal", "score"],
                              "auroc" = CelEsT_bench_RAW[CelEsT_bench_RAW$method == "mlm_estimate" & CelEsT_bench_RAW$metric == "auroc" & CelEsT_bench_RAW$net == "combo3equal", "score"])

exclude_stats_RAW_df <- data.frame(t(exclude_stats_RAW))

exclude_benchRAW_shufflestats <- lapply("EXCLUDEsample_benchRAW_shufflestats.tsv", function(file){
  
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
  
  tempout
  
})

exclude_benchRAW_shufflestats_df <- as.data.frame(t(unlist(exclude_benchRAW_shufflestats)))
row.names(exclude_benchRAW_shufflestats_df) <- exclude_benchRAW_shufflestats_df[, 1]
exclude_benchRAW_shufflestats_df <- exclude_benchRAW_shufflestats_df[, 2:ncol(exclude_benchRAW_shufflestats_df)]
exclude_benchRAW_shufflestats_df[] <- lapply(exclude_benchRAW_shufflestats_df, as.numeric)

pdf("graphics/CelEsT_exclude_RAW.pdf",
    height = 2.5,
    width = 2.5)

ggplot(exclude_stats_RAW_df, aes(x = auroc, y = auprc)) + 
  geom_point(alpha = 0.5,
             colour = "dark grey",
             size = 0.5) +
  theme_classic() +
  theme(strip.text = element_text(colour = 'white'),
        strip.background = element_rect(fill = "black"),
        legend.position = "none",
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_line(colour = "grey", linewidth = 0.2),
        panel.grid.minor = element_line(colour = "grey", linewidth = 0.1),
        axis.text = element_text(colour = "black")) +
  coord_cartesian(xlim = c(0.45, 0.8),
                  ylim = c(0.45, 0.8)) + 
  geom_point(data = CelEsT_bench_RAW_df,
             colour = "red",
             size = 2) +
  xlab("") + 
  ylab("") +
  geom_vline(xintercept = 0.5,
             linetype = "dashed",
             col = "black") +
  geom_hline(yintercept = 0.5,
             linetype = "dashed",
             col = "black") +
  geom_errorbar(data = exclude_benchRAW_shufflestats_df, aes(x = mean_auroc, y = mean_auprc,
                                                                ymin = (mean_auprc - sd_auprc),
                                                                ymax = (mean_auprc + sd_auprc)),
                alpha = 0.5, width = 0.01) +
  geom_errorbar(data = exclude_benchRAW_shufflestats_df, aes(x = mean_auroc, y = mean_auprc,
                                                                xmin = (mean_auroc - sd_auroc),
                                                                xmax = (mean_auroc + sd_auroc)), 
                alpha = 0.5, width = 0.01)


dev.off()





