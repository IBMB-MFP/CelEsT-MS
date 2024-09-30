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

packages <- c("ggplot2",
              "stringr",
              "ggrepel",
              "dplyr",
              "openxlsx")

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

wTF3 <- read.xlsx("input/WTF3.xlsx")

wTF3[, "family"] <- str_remove(wTF3$DBD, " - [0-9]{1,2} finger.*")
wTF3[, "family"] <- str_remove(wTF3$family, " - [0-9]{1,2} domain.*")

wTF3[, "family"] <- str_remove(wTF3$family, " x[0-9]{1}$")

wTF3[, "family"] <- str_remove(wTF3$family, " $")

wTF3[str_detect(wTF3$DBD, "^WH"), "family"] <- "WH"

wTF3[str_detect(wTF3$DBD, "^AT Hook"), "family"] <- "AT Hook"

wTF3[str_detect(wTF3$DBD, "^HD"), "family"] <- "HD"

wTF3[str_detect(wTF3$DBD, "^ZF C2H2"), "family"] <- "ZF - C2H2"

main_families <- names(table(wTF3$family)[order(table(wTF3$family), decreasing = TRUE)])[1:9]

wTF3[!wTF3$family %in% main_families, "family"] <- "Other"

benchmark_obs <- read.table("output/benchmark_observations.txt",
                  sep = '\t')

#### by TFs RAPToR-corrected DE ####

bysource <- read.table("output/benchmark_out/recall_bysource.tsv",
                       header = T,
                       fill = T)

bysource[, 3:7] <- bysource[, 1:5]
bysource <- bysource[, 3:7]

bysource <- bysource[bysource$metric %in% c("auroc", "auprc", "recall"), ]
bysource[, "name"] <- wormRef::Cel_genes[match(bysource$source, Cel_genes$sequence_name), "public_name"]
bysource[bysource$name == "atf-5", "name"] <- "atf-4"

bysource[, "bench_no"] <- as.numeric(unlist(table(benchmark_obs$target_name)[bysource$name]))

# bysource_auroc <- bysource[bysource$metric == "auroc", ]
bysource_auprc <- bysource[bysource$metric == "auprc", ]
bysource_auprc[, "auroc"] <- bysource[bysource$metric == "auroc", "score"]

bysource_auprc[, "family"] <- wTF3[match(bysource_auprc$source, wTF3$Sequence.name), "family"]

bysource_auprc[, "label"] <- toupper(bysource_auprc$name)
bysource_auprc[!(quantile(bysource_auprc$score, 0.1) & bysource_auprc$auroc < quantile(bysource_auprc$auroc, 0.1) | bysource_auprc$score > quantile(bysource_auprc$score, 0.9) & bysource_auprc$auroc > quantile(bysource_auprc$auroc, 0.9)), "label"] <- ""

bysource_auprc$family  <- factor(bysource_auprc$family, levels = c("ZF - NHR", "ZF - C2H2", "HD", "bHLH", "WH", "bZIP", "AT Hook", "Other"))

#### by TFs RAW DE ####

bysource_RAW <- read.table("output/benchmark_out/recall_bysource_RAW.tsv",
                           header = T,
                           fill = T)

bysource_RAW[, 3:7] <- bysource_RAW[, 1:5]
bysource_ <- bysource_RAW[, 3:7]

bysource_RAW <- bysource_RAW[bysource_RAW$metric %in% c("auroc", "auprc", "recall"), ]
bysource_RAW[, "name"] <- wormRef::Cel_genes[match(bysource_RAW$source, Cel_genes$sequence_name), "public_name"]
bysource_RAW[bysource_RAW$name == "atf-5", "name"] <- "atf-4"

bysource_RAW[, "bench_no"] <- unlist(table(benchmark_obs$target_name)[bysource_RAW$name])

# bysource_auroc <- bysource_RAW[bysource_RAW$metric == "auroc", ]
bysourceRAW_auprc <- bysource_RAW[bysource_RAW$metric == "auprc", ]
bysourceRAW_auprc[, "auroc"] <- bysource_RAW[bysource_RAW$metric == "auroc", "score"]

bysourceRAW_auprc[, "family"] <- wTF3[match(bysourceRAW_auprc$source, wTF3$Sequence.name), "family"]

bysourceRAW_auprc[, "label"] <- toupper(bysourceRAW_auprc$name)
bysourceRAW_auprc[!(quantile(bysourceRAW_auprc$score, 0.1) & bysourceRAW_auprc$auroc < quantile(bysourceRAW_auprc$auroc, 0.1) | bysourceRAW_auprc$score > quantile(bysourceRAW_auprc$score, 0.9) & bysourceRAW_auprc$auroc > quantile(bysourceRAW_auprc$auroc, 0.9)), "label"] <- ""

bysourceRAW_auprc$family  <- factor(bysourceRAW_auprc$family, levels = c("ZF - NHR", "ZF - C2H2", "HD", "bHLH", "WH", "bZIP", "AT Hook", "Other"))

bysourceRAW_auprc$bench_no <- as.numeric(unlist(bysourceRAW_auprc$bench_no))

# Mark TFs labelled in other plot

bysourceRAW_auprc[bysource_auprc$label != "", "label"] <- bysource_auprc[bysource_auprc$label != "", "label"] 
bysource_auprc[bysourceRAW_auprc$label != "", "label"] <- bysourceRAW_auprc[bysourceRAW_auprc$label != "", "label"] 

saveRDS(bysource_auprc, "plotdata/bench_byTF_RAPToR.rds")
saveRDS(bysourceRAW_auprc, "plotdata/bench_byTF_raw.rds")

pdf("graphics/byTF_RAPToR.pdf",
    height = 2.5,
    width = 2.5)

ggplot(bysource_auprc, aes(x = auroc, y = score, colour = family))+
  geom_point(aes(size = bench_no)) + 
  theme_classic() + 
  xlab("") + 
  ylab("") + 
  geom_label_repel(aes(label = label,
                       fill = family),
                   colour = "black",
                   label.size = 0.1,
                   max.overlaps = 20,
                   box.padding = 0.3,
                   label.padding = 0.1,
                   point.padding = 0,
                   size = 2) +
    scale_size_continuous(range = c(0.5, 3)) +
  scale_fill_manual(values = c("mediumpurple1", "olivedrab3", "dodger blue", "yellow3", "lightsalmon1", "orchid1", "burlywood4", "gray70")) + 
  scale_color_manual(values = c("mediumpurple1", "olivedrab3", "dodger blue", "yellow3", "lightsalmon1", "orchid1", "burlywood4", "gray70")) +
  theme(strip.text = element_text(colour = 'white'),
        strip.background = element_rect(fill = "black"),
        legend.position = "none",
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_line(colour = "grey", linewidth = 0.2),
        panel.grid.minor = element_line(colour = "grey", linewidth = 0.1),
        axis.text = element_text(colour = "black")) + 
  geom_vline(xintercept = 0.5,
             linetype = "dashed",
             col = "black") +
  geom_hline(yintercept = 0.5,
             linetype = "dashed",
             col = "black")

dev.off()


pdf("graphics/byTF_RAW.pdf",
    height = 2.5,
    width = 2.5)

ggplot(bysourceRAW_auprc, aes(x = auroc, y = score, colour = family))+
  geom_point(aes(size = bench_no)) + 
  theme_classic() + 
  xlab("") + 
  ylab("") + 
  geom_label_repel(aes(label = label,
                       fill = family),
                   colour = "black",
                   label.size = 0.1,
                   max.overlaps = 15,
                   box.padding = 0.3,
                   label.padding = 0.1,
                   point.padding = 0,
                   size = 2) +
  scale_size_continuous(range = c(0.5, 3)) +
  scale_fill_manual(values = c("mediumpurple1", "olivedrab3", "dodger blue", "yellow3", "lightsalmon1", "orchid1", "burlywood4", "gray70")) + 
  scale_color_manual(values = c("mediumpurple1", "olivedrab3", "dodger blue", "yellow3", "lightsalmon1", "orchid1", "burlywood4", "gray70")) +
  theme(strip.text = element_text(colour = 'white'),
        strip.background = element_rect(fill = "black"),
        legend.position = "none",
        panel.background = element_rect(fill = "white", colour = "black"),
        panel.grid.major = element_line(colour = "grey", linewidth = 0.2),
        panel.grid.minor = element_line(colour = "grey", linewidth = 0.1),
        axis.text = element_text(colour = "black")) + 
  geom_vline(xintercept = 0.5,
             linetype = "dashed",
             col = "black") +
  geom_hline(yintercept = 0.5,
             linetype = "dashed",
             col = "black")

dev.off()
