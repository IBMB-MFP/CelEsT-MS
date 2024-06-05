#start with empty workspace

rm(list = ls(all = TRUE))

# turn off scientific notation for plots

options(scipen=10000)

#### set working directory ####

# here create new folder and set working directory within it

dir.create("~/Cel_GRN_manuscript")
setwd("~/Cel_GRN_manuscript")

# create subfolders for input, output and graphics

dir.create("input")

# into input folder, add input files 

dir.create("output")
dir.create("graphics")

#### LOAD PACKAGES & FUNCTIONS ####

## First specify the packages of interest

packages <- c("biomaRt",
              "dplyr",
              "tidyr",
              "stringr",
              "gplots",
              "Hmisc",
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

biocmanager_packages <- c("decoupleR"
) 

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

#### DEFINE FUNCTIONS ####

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

meta.analysis.decouple <- function(DEdata,
                                   use_network = CelEsT){
  # x <- "daf2daf18"
  #   DEdata <- allstudiesDE_RAPToR_df[, str_detect(colnames(allstudiesDE_RAPToR_df), paste0(x, "$")), drop = FALSE]
  DEdata[is.na(DEdata)] <- 0
  
  DEdata_decouple <- decoupleR::decouple(mat = DEdata, 
                                         network = use_network,
                                         .source = "source",
                                         .target = "target",
                                         statistics = "mlm",
                                         args = list(mlm = list(.mor = "weight")),
                                         consensus_score = FALSE)
  
  DEdata_decouple_wide <- pivot_wider(data = DEdata_decouple[, c("source", "condition", "score")], names_from = "source", values_from = "score")
  
  DEdata_decouple_wide <- data.frame(t(DEdata_decouple_wide))
  
  colnames(DEdata_decouple_wide) <- DEdata_decouple_wide[1, ]
  DEdata_decouple_wide <- DEdata_decouple_wide[2:nrow(DEdata_decouple_wide), , drop = FALSE]
  
  DEdata_decouple_wide[] <- lapply(DEdata_decouple_wide, as.numeric)
  
  DEdata_decouple_zscores <- data.frame(matrix(nrow = nrow(DEdata_decouple_wide), ncol = ncol(DEdata_decouple_wide)))
  
  row.names(DEdata_decouple_zscores) <- row.names(DEdata_decouple_wide)
  colnames(DEdata_decouple_zscores) <- colnames(DEdata_decouple_wide)
  
  DEdata_decouple_zscores[] <- lapply(DEdata_decouple_wide, function(x){
    
    (x - mean(x)) / sd(x)
    
  })
  
  DEdata_zscores <- data.frame(matrix(nrow = nrow(DEdata), ncol = ncol(DEdata)))
  
  row.names(DEdata_zscores) <- row.names(DEdata)
  colnames(DEdata_zscores) <- colnames(DEdata)
  
  DEdata_zscores[] <- lapply(DEdata, function(x){
    
    (x - mean(x)) / sd(x)
    
  })
  
  # get geometric mean of p-values 
  
  DEdata_decouple_wide_p <- as.data.frame(pivot_wider(data = DEdata_decouple[, c("source", "condition", "p_value")], names_from = "source", values_from = "p_value"))
  row.names(DEdata_decouple_wide_p) <- unlist(DEdata_decouple_wide_p[, 1])
  DEdata_decouple_wide_p <- DEdata_decouple_wide_p[, 2:ncol(DEdata_decouple_wide_p)]
  
  DEdata_decouple_wide_p <- data.frame(t(DEdata_decouple_wide_p))
  
  DEdata_decouple_wide_p[] <- lapply(DEdata_decouple_wide_p, as.numeric)
  
  ## output to list
  
  list("decouple_z" = DEdata_decouple_zscores,
       "DE_z" = DEdata_zscores,
       "decouple_p" = DEdata_decouple_wide_p)
  
}

#### INPUT DATA ####

CelEsT <- read.table("output/GRNs/allthree_equalweights.txt",
                     sep = "\t",
                     header = TRUE)

orthCelEsT <- read.table("output/GRNs/orthCelEsT_equalweights.txt",
                     sep = "\t",
                     header = TRUE)

maxCelEsT <- read.table("output/GRNs/orthCelEsTMAXcov_equalweights.txt",
                         sep = "\t",
                         header = TRUE)

HS_DE <- read.table("output/HS_DE_RAPToR.txt",
                    sep = "\t",
                    header = TRUE)

HS_DE[is.na(HS_DE)] <- 0

MALE_DE <- read.table("output/MALE_DE_RAPToR.txt",
                      sep = "\t",
                      header = TRUE)

MALE_DE[is.na(MALE_DE)] <- 0


PA14_DE <- read.table("output/PA14_DE_RAPToR.txt",
                      sep = "\t",
                      header = TRUE)

PA14_DE[is.na(PA14_DE)] <- 0


fullset_TFs_BM <- readRDS("output/fullset_TFs_BM.rds")
fullset_TFs_BM[, "label"] <- fullset_TFs_BM$wormbase_locus
fullset_TFs_BM[fullset_TFs_BM$label == "", "label"] <- fullset_TFs_BM[fullset_TFs_BM$label == "", "wormbase_gseq"]

mainpal <- (colorRampPalette(c("blue", "white", "red"))(100))

#### HEAT SHOCK ####

HS_decouple <- meta.analysis.decouple(DEdata = HS_DE,
                                           use_network = CelEsT)

# HS_decouple_orth <- meta.analysis.decouple(DEdata = HS_DE,
#                                            network_use = orthCelEsT)
# 
# HS_decouple_max <- meta.analysis.decouple(DEdata = HS_DE,
#                                                 network_use = maxCelEsT)

# HS_orthmax_corrs <- lapply(1:ncol(HS_decouple_mlm_wide), function(i){
# 
#   shared <- base::intersect(base::intersect(row.names(HS_decouple_mlm_wide), row.names(HS_decouple_orth_mlm_wide)), row.names(HS_decouple_max_mlm_wide))
#   
#   temptogether <- cbind(HS_decouple_mlm_wide[shared, i], HS_decouple_orth_mlm_wide[shared, i], HS_decouple_max_mlm_wide[shared, i])
#   colnames(temptogether) <- c("CelEsT", "orthCelEsT", "maxCelEsT")
#   
#   temp_cormat <- rcorr(temptogether, type = "pearson")
# 
#   temp_cormat$r
# 
#   })

# check for outliers, rerun without outliers if needed

heatmap.2(as.matrix(HS_decouple[["decouple_z"]]),
          col = mainpal,
          dendrogram = "column",
          density = "none",
          key = FALSE,
          trace = "none")

pdf("graphics/HS_samples_heatmap.pdf",
    height = 4, 
    width = 4)

heatmap.2(as.matrix(HS_decouple[["decouple_z"]]),
          col = mainpal,
          dendrogram = "column",
          density = "none",
          key = FALSE,
          trace = "none",
          breaks = seq(from = -10, to = 10, length.out = 101))

dev.off()

pdf("graphics/HS_samples_heatmap_LEGEND.pdf")

heatmap.2(as.matrix(HS_decouple_zscores),
          col = mainpal,
          dendrogram = "column",
          density = "none",
          breaks = seq(from = -10, to = 10, length.out = 101))

dev.off()

rcorr(as.matrix(HS_decouple_zscores))

# weird! 0.00 correlation? Yes. Wow.

cor.test(HS_decouple_zscores$GSE122015,HS_decouple_zscores$GSE162064 )

# exclude outlier sample; GSE122015

HS_decouple_zscores <- HS_decouple_zscores[, colnames(HS_decouple_zscores) != "GSE122015"]

# compare to expression

HS_DE_zscores <- data.frame(matrix(nrow = nrow(HS_DE), ncol = ncol(HS_DE)))

row.names(HS_DE_zscores) <- row.names(HS_DE)
colnames(HS_DE_zscores) <- colnames(HS_DE)

HS_DE_zscores[] <- lapply(HS_DE, function(x){
  
  (x - mean(x)) / sd(x)
  
})

HS_TF_DE_meanz <- rowMeans(HS_DE_zscores)

HS_TF_DE_meanz <- HS_TF_DE_meanz[row.names(HS_decouple_zscores)]

HS_DE_vs_activity_plot <- data.frame("DE" = HS_TF_DE_meanz,
                                     "act" = HS_meanz[names(HS_TF_DE_meanz)],
                                     "label" = fullset_TFs_BM[match(names(HS_TF_DE_meanz), fullset_TFs_BM$wormbase_gseq), "label"]
)

HS_DE_vs_activity_plot <- HS_DE_vs_activity_plot[!is.na(HS_DE_vs_activity_plot$DE), ]

HS_DE_vs_activity_plot[abs(HS_DE_vs_activity_plot$DE) < 1.96 & abs(HS_DE_vs_activity_plot$act)  < 1.96, "label"] <- ""

cor.test(HS_DE_vs_activity_plot$DE,
         HS_DE_vs_activity_plot$act)

pdf("graphics/HS_act_vs_DE.pdf",
    height = 2.5,
    width = 2.5)

ggplot(HS_DE_vs_activity_plot, aes(x = DE, y = act, label = label)) + 
  geom_point(size = 0.3) + 
  geom_smooth(method = "lm",
              se = FALSE,
              col = "red",
              size = 0.5) + 
  geom_label_repel(size = 1.8,
                   label.padding = 0.1) + 
  xlab("Differential expression z-score") + 
  ylab("TF activity z-score") +
  geom_hline(yintercept = 0,
             linetype = "dashed",
             col = "grey",
             size = 0.3) +
  geom_vline(xintercept = 0,
             linetype = "dashed",
             col = "grey",
             size = 0.3) +
  theme_classic() + 
  theme(axis.text = element_text(colour = "black"),
        axis.title = element_text(size = 10))

dev.off()

# get geometric mean of p-values 

HS_decouple_mlm_wide_p <- as.data.frame(pivot_wider(data = HS_decouple[HS_decouple$statistic == "mlm", c("source", "condition", "p_value")], names_from = "source", values_from = "p_value"))

HS_decouple_mlm_wide_p <- HS_decouple_mlm_wide_p[, 2:ncol(HS_decouple_mlm_wide_p)]

HS_decouple_mlm_wide_p <- HS_decouple_mlm_wide_p[row.names(HS_decouple_mlm_wide_p) != "GSE122015", ]

HS_decouple_mlm_wide_p <- data.frame(t(HS_decouple_mlm_wide_p))

HS_decouple_mlm_wide_p[] <- lapply(HS_decouple_mlm_wide_p, as.numeric)

HS_decouple_mlm_wide_pgeomean <- apply(HS_decouple_mlm_wide_p, 1, gm_mean)

HS_decouple_mlm_wide_pgeomean[order(HS_decouple_mlm_wide_pgeomean)]

HS_forbubbleplot <- data.frame(mean_z = HS_meanz[names(HS_decouple_mlm_wide_pgeomean)],
                                 p_val = -log10(HS_decouple_mlm_wide_pgeomean),
                                 label = fullset_TFs_BM[match(names(HS_decouple_mlm_wide_pgeomean), fullset_TFs_BM$wormbase_gseq), "label"])

HS_forbubbleplot[abs(HS_forbubbleplot$p_val) < 2.2, "label"] <- ""

saveRDS(HS_forbubbleplot,
     "plotdata/HS_forbubbleplot.rds")

pdf("graphics/HS_bubbleplot.pdf",
    height = 2.5,
    width = 2.5)

ggplot(aes(x = mean_z, y = p_val, size = 2^p_val, label = toupper(label)), data = HS_forbubbleplot) + 
  geom_point(col = "red", alpha = 0.7) + 
  theme_classic() + 
  geom_label_repel(size = 1.8, max.overlaps = 40,
                   label.padding = 0.1) +
  geom_vline(xintercept = 0, linetype = "dashed", col = "grey") + 
  theme(legend.position = "none",
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10)) + 
  ylab(substitute("-log"[10]~"(mean p-value)")) + 
  xlab("TF activity (mean z-score)") +
  coord_cartesian(xlim = c(-3, 3),
                  ylim = c(0, 6))

dev.off()

#### MALES ####

MALE_decouple <- decoupleR::decouple(mat = MALE_DE, 
                                   network = CelEsT,
                                   .source = "source",
                                   .target = "target",
                                   statistics = "mlm",
                                   args = list(mlm = list(.mor = "weight")))

MALE_decouple_orth <- decoupleR::decouple(mat = MALE_DE, 
                                     network = orthCelEsT,
                                     .source = "source",
                                     .target = "target",
                                     statistics = "mlm",
                                     args = list(mlm = list(.mor = "weight")))

MALE_decouple_max <- decoupleR::decouple(mat = MALE_DE, 
                                          network = maxCelEsT,
                                          .source = "source",
                                          .target = "target",
                                          statistics = "mlm",
                                          args = list(mlm = list(.mor = "weight")))

MALE_decouple_mlm_wide <- make.decouple.wide(MALE_decouple)
MALE_decouple_orth_mlm_wide <- make.decouple.wide(MALE_decouple_orth)
MALE_decouple_max_mlm_wide <- make.decouple.wide(MALE_decouple_max)

lapply(1:ncol(MALE_decouple_mlm_wide), function(i){
  
  shared <- base::intersect(base::intersect(row.names(MALE_decouple_mlm_wide), row.names(MALE_decouple_orth_mlm_wide)), row.names(MALE_decouple_max_mlm_wide))
  
  temptogether <- cbind(MALE_decouple_mlm_wide[shared, i], MALE_decouple_orth_mlm_wide[shared, i], MALE_decouple_max_mlm_wide[shared, i])
  colnames(temptogether) <- c("CelEsT", "orthCelEsT", "maxCelEsT")
  
  temp_cormat <- rcorr(temptogether, type = "pearson")
  
  temp_cormat$r
  
})

MALE_decouple_zscores <- data.frame(matrix(nrow = nrow(MALE_decouple_mlm_wide), ncol = ncol(MALE_decouple_mlm_wide)))

row.names(MALE_decouple_zscores) <- row.names(MALE_decouple_mlm_wide)
colnames(MALE_decouple_zscores) <- colnames(MALE_decouple_mlm_wide)

MALE_decouple_zscores[] <- lapply(MALE_decouple_mlm_wide, function(x){
  
  (x - mean(x)) / sd(x)
  
})

# names(MALE_meanz) <- fullset_TFs_BM[match(names(MALE_meanz) , fullset_TFs_BM$wormbase_gseq), "label"]

heatmap.2(as.matrix(MALE_decouple_zscores),
          col = mainpal)
rcorr(as.matrix(MALE_decouple_mlm_wide))

pdf("graphics/MALE_samples_heatmap.pdf",
    height = 4, 
    width = 4)

heatmap.2(as.matrix(MALE_decouple_zscores),
          col = mainpal,
          dendrogram = "column",
          density = "none",
          key = FALSE,
          trace = "none",
          breaks = seq(from = -10, to = 10, length.out = 101))

dev.off()

pdf("graphics/MALE_samples_heatmap_LEGEND.pdf",
    height = 4, 
    width = 4)

heatmap.2(as.matrix(MALE_decouple_zscores),
          col = mainpal,
          dendrogram = "column",
          density = "none",
          breaks = seq(from = -10, to = 10, length.out = 101))

dev.off()

# ONE AGREES POORLY WITH OTHERS. 447
MALE_decouple_zscores <- MALE_decouple_zscores[, colnames(MALE_decouple_zscores) != "GSE222447"]

MALE_meanz <- rowMeans(MALE_decouple_zscores)
MALE_meanz <- MALE_meanz[order(MALE_meanz, decreasing = TRUE)]

MALE_DE_zscores <- data.frame(matrix(nrow = nrow(MALE_DE), ncol = ncol(MALE_DE)))

row.names(MALE_DE_zscores) <- row.names(MALE_DE)
colnames(MALE_DE_zscores) <- colnames(MALE_DE)

MALE_DE_zscores[] <- lapply(MALE_DE, function(x){
  
  (x - mean(x)) / sd(x)
  
})

MALE_TF_DE_meanz <- rowMeans(MALE_DE_zscores)

MALE_TF_DE_meanz <- MALE_TF_DE_meanz[row.names(MALE_decouple_zscores)]

MALE_DE_vs_activity_plot <- data.frame("DE" = MALE_TF_DE_meanz,
                                     "act" = MALE_meanz[names(MALE_TF_DE_meanz)],
                                     "label" = fullset_TFs_BM[match(names(MALE_TF_DE_meanz), fullset_TFs_BM$wormbase_gseq), "label"]
)

MALE_DE_vs_activity_plot <- MALE_DE_vs_activity_plot[!is.na(MALE_DE_vs_activity_plot$DE), ]

MALE_DE_vs_activity_plot[abs(MALE_DE_vs_activity_plot$DE) < 1.9 & abs(MALE_DE_vs_activity_plot$act)  < 2.3, "label"] <- ""

cor.test(MALE_DE_vs_activity_plot$DE,
         MALE_DE_vs_activity_plot$act)

pdf("graphics/MALE_act_vs_DE.pdf",
    height = 2.5,
    width = 2.5)

ggplot(MALE_DE_vs_activity_plot, aes(x = DE, y = act, label = label)) + 
  geom_point(size = 0.3) + 
  geom_smooth(method = "lm",
              se = FALSE,
              col = "red",
              size = 0.5) + 
  geom_label_repel(size = 1.8,
                   max.overlaps = 40,
                   label.padding = 0.1) + 
  xlab("Differential expression z-score") + 
  ylab("TF activity z-score") +
  geom_hline(yintercept = 0,
             linetype = "dashed",
             col = "grey",
             size = 0.3) +
  geom_vline(xintercept = 0,
             linetype = "dashed",
             col = "grey",
             size = 0.3) +
  theme_classic() + 
  theme(axis.text = element_text(colour = "black"),
        axis.title = element_text(size = 10))

dev.off()

# get geometric mean of p-values 

MALE_decouple_mlm_wide_p <- as.data.frame(pivot_wider(data = MALE_decouple[MALE_decouple$statistic == "mlm", c("source", "condition", "p_value")], names_from = "source", values_from = "p_value"))
row.names(MALE_decouple_mlm_wide_p) <- unlist(MALE_decouple_mlm_wide_p[, 1])
MALE_decouple_mlm_wide_p <- MALE_decouple_mlm_wide_p[, 2:ncol(MALE_decouple_mlm_wide_p)]

MALE_decouple_mlm_wide_p <- data.frame(t(MALE_decouple_mlm_wide_p))

MALE_decouple_mlm_wide_p[] <- lapply(MALE_decouple_mlm_wide_p, as.numeric)

MALE_decouple_mlm_wide_p <- MALE_decouple_mlm_wide_p[, colnames(MALE_decouple_mlm_wide_p) != "GSE222447"]

MALE_decouple_mlm_wide_pgeomean <- apply(MALE_decouple_mlm_wide_p, 1, gm_mean)

MALE_decouple_mlm_wide_pgeomean[order(MALE_decouple_mlm_wide_pgeomean)]

MALE_forbubbleplot <- data.frame(mean_z = MALE_meanz[names(MALE_decouple_mlm_wide_pgeomean)],
                                 p_val = -log10(MALE_decouple_mlm_wide_pgeomean),
                                 label = fullset_TFs_BM[match(names(MALE_decouple_mlm_wide_pgeomean), fullset_TFs_BM$wormbase_gseq), "label"])

MALE_forbubbleplot[abs(MALE_forbubbleplot$p_val) < 2.7, "label"] <- ""

# # mark out genes in DREAM
# MALE_forbubbleplot[, "pointcol"] <- "red"
# MALE_forbubbleplot[MALE_forbubbleplot$label %in% c("lin-15B", "efl-1", "efl-2", "dpl-1"), "pointcol"] <- "blue"
# MALE_forbubbleplot[, "textcol"] <- "black"
# MALE_forbubbleplot[MALE_forbubbleplot$label %in% c("lin-15B", "efl-1", "efl-2", "dpl-1"), "textcol"] <- "blue"

saveRDS(MALE_forbubbleplot,
        "plotdata/MALE_forbubbleplot.rds")

pdf("graphics/MALE_bubbleplot.pdf",
    height = 2.5,
    width = 2.5)

ggplot(aes(x = mean_z, y = p_val, size = 2^p_val, label = toupper(label)), data = MALE_forbubbleplot) + 
  geom_point(colour = "red",  alpha = 0.7) + 
  theme_classic() + 
  geom_label_repel(size = 1.8, max.overlaps = 40, 
                   label.padding = 0.1) +
  geom_vline(xintercept = 0, linetype = "dashed", col = "grey") + 
  theme(legend.position = "none",
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10)) + 
  ylab(substitute("-log"[10]~"(mean p-value)")) + 
  xlab("TF activity (mean z-score)") +
  coord_cartesian(xlim = c(-4, 4),
                  ylim = c(0, 10))

dev.off()

#### PA14 ####

PA14_decouple <- decoupleR::decouple(mat = PA14_DE, 
                                     network = CelEsT,
                                     .source = "source",
                                     .target = "target",
                                     statistics = "mlm",
                                     args = list(mlm = list(.mor = "weight")))

PA14_decouple_orth <- decoupleR::decouple(mat = PA14_DE, 
                                     network = orthCelEsT,
                                     .source = "source",
                                     .target = "target",
                                     statistics = "mlm",
                                     args = list(mlm = list(.mor = "weight")))

PA14_decouple_max <- decoupleR::decouple(mat = PA14_DE, 
                                          network = maxCelEsT,
                                          .source = "source",
                                          .target = "target",
                                          statistics = "mlm",
                                          args = list(mlm = list(.mor = "weight")))

PA14_decouple_mlm_wide <- make.decouple.wide(PA14_decouple)
PA14_decouple_max_mlm_wide <- make.decouple.wide(PA14_decouple_max)
PA14_decouple_orth_mlm_wide <- make.decouple.wide(PA14_decouple_orth)

lapply(1:ncol(PA14_decouple_mlm_wide), function(i){
  
  shared <- base::intersect(base::intersect(row.names(PA14_decouple_mlm_wide), row.names(PA14_decouple_orth_mlm_wide)), row.names(PA14_decouple_max_mlm_wide))
  
  temptogether <- cbind(PA14_decouple_mlm_wide[shared, i], PA14_decouple_orth_mlm_wide[shared, i], PA14_decouple_max_mlm_wide[shared, i])
  colnames(temptogether) <- c("CelEsT", "orthCelEsT", "maxCelEsT")
  
  temp_cormat <- rcorr(temptogether, type = "pearson")
  
  temp_cormat$r
  
})

## there are two that are obviously the same dataset

# so exclude GSE146406

PA14_decouple_zscores <- data.frame(matrix(nrow = nrow(PA14_decouple_mlm_wide[, colnames(PA14_decouple_mlm_wide) != "GSE146406"]), ncol = ncol(PA14_decouple_mlm_wide[, colnames(PA14_decouple_mlm_wide) != "GSE146406"])))

row.names(PA14_decouple_zscores) <- row.names(PA14_decouple_mlm_wide)
colnames(PA14_decouple_zscores) <- colnames(PA14_decouple_mlm_wide)[colnames(PA14_decouple_mlm_wide) != "GSE146406"]

PA14_decouple_zscores[] <- lapply(PA14_decouple_mlm_wide[, colnames(PA14_decouple_mlm_wide) != "GSE146406"], function(x){
  
  (x - mean(x)) / sd(x)
  
})

heatmap.2(as.matrix(PA14_decouple_zscores),
          col = mainpal)

PA14_meanz <- rowMeans(PA14_decouple_zscores)
PA14_meanz <- PA14_meanz[order(PA14_meanz, decreasing = TRUE)]

Hmisc::rcorr(as.matrix(PA14_DE))
rowMeans(Hmisc::rcorr(as.matrix(PA14_decouple_mlm_wide))[[1]])



# names(PA14_meanz) <- fullset_TFs_BM[match(names(PA14_meanz) , fullset_TFs_BM$wormbase_gseq), "label"]

PA14_DE_zscores <- data.frame(matrix(nrow = nrow(PA14_DE), ncol = ncol(PA14_DE) -1))

row.names(PA14_DE_zscores) <- row.names(PA14_DE)
colnames(PA14_DE_zscores) <- colnames(PA14_DE)[colnames(PA14_DE) != "GSE146406"]

PA14_DE_zscores[] <- lapply(PA14_DE[, colnames(PA14_DE) != "GSE146406"], function(x){
  
  (x - mean(x)) / sd(x)
  
})

PA14_TF_DE_meanz <- rowMeans(PA14_DE_zscores)


PA14_TF_DE_meanz <- PA14_TF_DE_meanz[row.names(PA14_decouple_zscores)]

pdf("graphics/PA14_samples_heatmap.pdf",
    height = 4, 
    width = 4)

heatmap.2(as.matrix(PA14_decouple_zscores),
          col = mainpal,
          trace = "none",
          density = "none",
          dendrogram = "column",
          key = FALSE,
          breaks = seq(from = -10, to = 10, length.out = 101)
          )

dev.off()

# get geometric mean of p-values 

PA14_decouple_mlm_wide_p <- as.data.frame(pivot_wider(data = PA14_decouple_mlm[, c("source", "condition", "p_value")], names_from = "source", values_from = "p_value"))
row.names(PA14_decouple_mlm_wide_p) <- unlist(PA14_decouple_mlm_wide_p[, 1])
PA14_decouple_mlm_wide_p <- PA14_decouple_mlm_wide_p[, 2:ncol(PA14_decouple_mlm_wide_p)]

PA14_decouple_mlm_wide_p <- data.frame(t(PA14_decouple_mlm_wide_p))

PA14_decouple_mlm_wide_p[] <- lapply(PA14_decouple_mlm_wide_p, as.numeric)

PA14_decouple_mlm_wide_pgeomean <- apply(PA14_decouple_mlm_wide_p[, colnames(PA14_decouple_mlm_wide_p) != "GSE146406"], 1, gm_mean)

PA14_forbubbleplot <- data.frame(mean_z = PA14_meanz[names(PA14_decouple_mlm_wide_pgeomean)],
                                 p_val = -log10(PA14_decouple_mlm_wide_pgeomean),
                                 label = fullset_TFs_BM[match(names(PA14_decouple_mlm_wide_pgeomean), fullset_TFs_BM$wormbase_gseq), "label"])

PA14_forbubbleplot[abs(PA14_forbubbleplot$p_val) < 3.5, "label"] <- ""

saveRDS(PA14_forbubbleplot,
        "plotdata/PA14_forbubbleplot.rds")

pdf("graphics/PA14_bubbleplot.pdf",
    height = 2.5,
    width = 2.5)

ggplot(aes(x = mean_z, y = p_val, size = 2^p_val, label = toupper(label)), data = PA14_forbubbleplot) + 
  geom_point(col = "red", alpha = 0.7) + 
  theme_classic() + 
  geom_label_repel(size = 1.8, max.overlaps = 40,
                   label.padding = 0.1) +
  geom_vline(xintercept = 0, linetype = "dashed", col = "grey") + 
  theme(legend.position = "none",
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10)) + 
  ylab(substitute("-log"[10]~"(mean p-value)")) + 
  xlab("TF activity (mean z-score)") +
  coord_cartesian(xlim = c(-7, 7),
                  ylim = c(0, 32))

dev.off()

PA14_DE_vs_activity_plot <- data.frame("DE" = PA14_TF_DE_meanz,
                                       "act" = PA14_meanz[names(PA14_TF_DE_meanz)],
                                       "label" = fullset_TFs_BM[match(names(PA14_TF_DE_meanz), fullset_TFs_BM$wormbase_gseq), "label"]
)

PA14_DE_vs_activity_plot <- PA14_DE_vs_activity_plot[!is.na(PA14_DE_vs_activity_plot$DE), ]

PA14_DE_vs_activity_plot[abs(PA14_DE_vs_activity_plot$DE) < 2.1 & abs(PA14_DE_vs_activity_plot$act)  < 2.1, "label"] <- ""

cor.test(PA14_DE_vs_activity_plot$DE,
         PA14_DE_vs_activity_plot$act)$estimate

pdf("graphics/PA14_act_vs_DE.pdf",
    height = 2.5,
    width = 2.5)

ggplot(PA14_DE_vs_activity_plot, aes(x = DE, y = act, label = label)) + 
  geom_point(size = 0.3) + 
  geom_smooth(method = "lm",
              se = FALSE,
              col = "red",
              size = 0.5) + 
  geom_label_repel(size = 1.8,
                   label.padding = 0.1,
                   max.overlaps = 40) + 
  xlab("Differential expression z-score") + 
  ylab("TF activity z-score") +
  geom_hline(yintercept = 0,
             linetype = "dashed",
             col = "grey",
             size = 0.3) +
  geom_vline(xintercept = 0,
             linetype = "dashed",
             col = "grey",
             size = 0.3) +
  theme_classic() + 
  theme(axis.text = element_text(colour = "black"),
        axis.title = element_text(size = 10))

dev.off()
