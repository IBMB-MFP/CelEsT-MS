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

dir.create("output/STREME_out")

#### LOAD PACKAGES & FUNCTIONS ####

## First specify the packages of interest

packages <- c("dplyr",
              "stringr",
              "openxlsx",
              "ggrepel",
              "gplots",
              "Hmisc")

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

gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

map2color <- function(x, 
                      pal,
                      limits = NULL){
  
  if(is.null(limits)){
    limits = range(x)
  }
  
  pal[findInterval(x, seq(limits[1], limits[2], length.out = length(pal) + 1), all.inside = TRUE)]
  
}

#### INPUT DATA ####

mainpal <- (colorRampPalette(c("blue", "white", "red"))(100))

daf16_SRA <- read.xlsx("input/IIS_SRA.xlsx",
                       sheet = 1)

daf18_SRA <- read.xlsx("input/IIS_SRA.xlsx",
                       sheet = 2)

daf2_SRA <- read.xlsx("input/IIS_SRA.xlsx",
                      sheet = 3)

daf2daf16_SRA <- read.xlsx("input/IIS_SRA.xlsx",
                           sheet = 4)

daf2daf18_SRA <- read.xlsx("input/IIS_SRA.xlsx",
                           sheet = 5)

all_SRA <- do.call(rbind, lapply(c("daf16_SRA",
                                   "daf18_SRA",
                                   "daf2_SRA",
                                   "daf2daf16_SRA",
                                   "daf2daf18_SRA"), function(thisgene){
                                     
                                     thisdata <- get(thisgene)
                                     
                                     thisdata[, "genotype"] <- str_remove(thisgene, "_SRA")
                                     
                                     thisdata
                                     
                                   }))

allstudiesDE_RAPToR_df <- read.table("output/IIS_DE_RAPToR.txt",
            sep = "\t",
            header = TRUE)

fullset_TFs_BM <- readRDS("output/fullset_TFs_BM.rds")
fullset_TFs_BM[, "label"] <- fullset_TFs_BM$wormbase_locus
fullset_TFs_BM[fullset_TFs_BM$label == "", "label"] <- fullset_TFs_BM[fullset_TFs_BM$label == "", "wormbase_gseq"]

CelEsT <- read.table("output/GRNs/allthree_equalweights.txt",
                     sep = "\t",
                     header = TRUE)

orthCelEsT <- read.table("output/GRNs/orthCelEsT_equalweights.txt",
                         sep = "\t",
                         header = TRUE)

maxCelEsT <- read.table("output/GRNs/orthCelEsTMAXcov_equalweights.txt",
                        sep = "\t",
                        header = TRUE)


#### DAF2 ####

genotypelist <- unique(all_SRA$genotype)

decouplelist <- lapply(genotypelist[c(1, 3:4)], function(x){print(x)
  meta.analysis.decouple(allstudiesDE_RAPToR_df[, str_remove(colnames(allstudiesDE_RAPToR_df), "GSE[0-9]+") == x, drop = FALSE])})
names(decouplelist) <- genotypelist[c(1, 3:4)]

decouplelist2 <- lapply(genotypelist[c(2, 5)], function(x){print(x)
  meta.analysis.decouple(allstudiesDE_RAPToR_df[, str_remove(colnames(allstudiesDE_RAPToR_df), "GSE[0-9]+") == x, drop = FALSE])})
names(decouplelist2) <- genotypelist[c(2, 5)]

# check for outliers, rerun without outliers if needed

pdf("graphics/daf2_samples_heatmap.pdf",
    height = 4, 
    width = 4)

heatmap.2(as.matrix(decouplelist[["daf2"]][["decouple_z"]]),
          col = mainpal,
          dendrogram = "column",
          density = "none",
          key = FALSE,
          trace = "none",
          breaks = seq(from = -10, to = 10, length.out = 101))

dev.off()

# exclude "GSE36041daf2"

rcorr(as.matrix(decouplelist[["daf2"]][["decouple_z"]]))
rowMeans(rcorr(as.matrix(decouplelist[["daf2"]][["decouple_z"]]))[[1]])

daf2_forbubbleplot <- data.frame(mean_z = rowMeans(decouplelist[["daf2"]][["decouple_z"]][, colnames(decouplelist[["daf2"]][["decouple_z"]]) != "GSE36041daf2"]),
           p_val = -log10(apply(decouplelist[["daf2"]][["decouple_p"]][, colnames(decouplelist[["daf2"]][["decouple_z"]]) != "GSE36041daf2"], 1, gm_mean)),
           label = toupper(fullset_TFs_BM[match(row.names(decouplelist[["daf2"]][["decouple_p"]][, colnames(decouplelist[["daf2"]][["decouple_z"]]) != "GSE36041daf2"]), fullset_TFs_BM$wormbase_gseq), "label"]))

daf2_forbubbleplot[abs(daf2_forbubbleplot$p_val) < 3.5, "label"] <- ""

saveRDS(daf2_forbubbleplot,
        "plotdata/daf2_forbubbleplot.rds")
# daf2_forbubbleplot <- readRDS("plotdata/daf2_forbubbleplot.rds")

daf2colours_for_plot <- map2color(daf2_forbubbleplot$p_val, pal = colorRampPalette(c("grey", "grey", "grey", "red", "red", "red"))(100))

pdf("graphics/daf2_bubbleplot.pdf",
    height = 2.5,
    width = 2.5)

ggplot(aes(x = mean_z, y = p_val, size = 2^p_val, label = label), data = daf2_forbubbleplot) + 
  geom_point(col = daf2colours_for_plot, alpha = 0.7) + 
  theme_classic() + 
  geom_label_repel(size = 1.8, max.overlaps = 40, label.padding = 0.1) +
  geom_vline(xintercept = 0, linetype = "dashed", col = "grey") + 
  theme(legend.position = "none",
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10)) + 
  ylab(substitute("-log"[10]~"(mean p-value)")) + 
  xlab("TF activity (mean z-score)")

dev.off()


daf2_DE_vs_activity_plot <- data.frame("DE" = rowMeans(decouplelist[["daf2"]][["DE_z"]][, colnames(decouplelist[["daf2"]][["DE_z"]]) != "GSE36041daf2"])[row.names(decouplelist[["daf2"]][["decouple_z"]])],
                                     "act" = rowMeans(decouplelist[["daf2"]][["decouple_z"]][, colnames(decouplelist[["daf2"]][["decouple_z"]]) != "GSE36041daf2"]),
                                     "label" = fullset_TFs_BM[match(row.names(decouplelist[["daf2"]][["decouple_z"]]), fullset_TFs_BM$wormbase_gseq), "label"])

daf2_DE_vs_activity_plot <- daf2_DE_vs_activity_plot[!is.na(daf2_DE_vs_activity_plot$DE), ]

daf2_DE_vs_activity_plot[abs(daf2_DE_vs_activity_plot$DE) < 1.96 & abs(daf2_DE_vs_activity_plot$act)  < 1.96, "label"] <- ""

cor.test(daf2_DE_vs_activity_plot$DE,
     daf2_DE_vs_activity_plot$act)

saveRDS(daf2_DE_vs_activity_plot,
        "plotdata/daf2_DE_vs_activity_plot.rds")

pdf("graphics/daf2_act_vs_DE.pdf",
    height = 2.5,
    width = 2.5)

ggplot(daf2_DE_vs_activity_plot, aes(x = DE, y = act, label = label)) + 
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

#### DAF2 DAF16 ####

pdf("graphics/daf2daf16_samples_heatmap.pdf",
    height = 4, 
    width = 4)

heatmap.2(as.matrix(decouplelist[["daf2daf16"]][["decouple_z"]]),
          col = mainpal,
          dendrogram = "column",
          density = "none",
          key = FALSE,
          trace = "none",
          breaks = seq(from = -10, to = 10, length.out = 101))

dev.off()

daf2daf16_forbubbleplot <- data.frame(mean_z = rowMeans(decouplelist[["daf2daf16"]][["decouple_z"]]),
                                 p_val = -log10(apply(decouplelist[["daf2daf16"]][["decouple_p"]], 1, gm_mean)),
                                 label = toupper(fullset_TFs_BM[match(row.names(decouplelist[["daf2daf16"]][["decouple_p"]]), fullset_TFs_BM$wormbase_gseq), "label"]))

daf2daf16_forbubbleplot[abs(daf2daf16_forbubbleplot$p_val) < 5, "label"] <- ""

saveRDS(daf2daf16_forbubbleplot,
        "plotdata/daf2daf16_forbubbleplot.rds")

# daf2daf16_forbubbleplot <- readRDS("plotdata/daf2daf16_forbubbleplot.rds")

daf2daf16colours_for_plot <- map2color(daf2daf16_forbubbleplot$p_val, pal = colorRampPalette(c("grey", "grey", "grey", "red", "red", "red"))(100))


pdf("graphics/daf2daf16_bubbleplot.pdf",
    height = 2.5,
    width = 2.5)

ggplot(aes(x = mean_z, y = p_val, size = 2^p_val, label = label), data = daf2daf16_forbubbleplot) + 
  geom_point(col = daf2daf16colours_for_plot, alpha = 0.7) + 
  theme_classic() + 
  geom_label_repel(size = 1.8, max.overlaps = 40, label.padding = 0.1) +
  geom_vline(xintercept = 0, linetype = "dashed", col = "grey") + 
  theme(legend.position = "none",
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10)) + 
  ylab(substitute("-log"[10]~"(mean p-value)")) + 
  xlab("TF activity (mean z-score)") 

dev.off()

daf2daf16_DE_vs_activity_plot <- data.frame("DE" = rowMeans(decouplelist[["daf2daf16"]][["DE_z"]][, colnames(decouplelist[["daf2daf16"]][["DE_z"]]) != "GSE36041daf2daf16"])[row.names(decouplelist[["daf2daf16"]][["decouple_z"]])],
                                       "act" = rowMeans(decouplelist[["daf2daf16"]][["decouple_z"]][, colnames(decouplelist[["daf2daf16"]][["decouple_z"]]) != "GSE36041daf2daf16"]),
                                       "label" = fullset_TFs_BM[match(row.names(decouplelist[["daf2daf16"]][["decouple_z"]]), fullset_TFs_BM$wormbase_gseq), "label"])

daf2daf16_DE_vs_activity_plot <- daf2daf16_DE_vs_activity_plot[!is.na(daf2daf16_DE_vs_activity_plot$DE), ]

daf2daf16_DE_vs_activity_plot[abs(daf2daf16_DE_vs_activity_plot$DE) < 1.96 & abs(daf2daf16_DE_vs_activity_plot$act)  < 1.96, "label"] <- ""

saveRDS(daf2daf16_DE_vs_activity_plot,
        "plotdata/daf2daf16_DE_vs_activity_plot.rds")

pdf("graphics/daf2daf16_act_vs_DE.pdf",
    height = 2.5,
    width = 2.5)

ggplot(daf2daf16_DE_vs_activity_plot, aes(x = DE, y = act, label = label)) + 
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

#### DAF-16 ####

pdf("graphics/daf16_samples_heatmap.pdf",
    height = 4, 
    width = 4)

heatmap.2(as.matrix(decouplelist[["daf16"]][["decouple_z"]]),
          col = mainpal,
          dendrogram = "column",
          density = "none",
          key = FALSE,
          trace = "none",
          breaks= seq(from = -10, to = 10, length.out = 101))

dev.off()

daf16_forbubbleplot <- data.frame(mean_z = rowMeans(decouplelist[["daf16"]][["decouple_z"]][colnames(decouplelist[["daf16"]][["decouple_z"]]) %in% c("GSE240821daf16", "GSE108848daf16")]),
                                 p_val = -log10(apply(decouplelist[["daf16"]][["decouple_p"]][colnames(decouplelist[["daf16"]][["decouple_z"]]) %in% c("GSE240821daf16", "GSE108848daf16")], 1, gm_mean)),
                                 label = toupper(fullset_TFs_BM[match(row.names(decouplelist[["daf16"]][["decouple_p"]][colnames(decouplelist[["daf16"]][["decouple_z"]]) %in% c("GSE240821daf16", "GSE108848daf16")]), fullset_TFs_BM$wormbase_gseq), "label"]))

daf16_forbubbleplot[abs(daf16_forbubbleplot$p_val) < 3, "label"] <- ""

saveRDS(daf16_forbubbleplot,
        "plotdata/daf16_forbubbleplot.rds")

daf16_forbubbleplot <- readRDS("plotdata/daf16_forbubbleplot.rds")

daf16colours_for_plot <- map2color(daf16_forbubbleplot$p_val, pal = colorRampPalette(c("grey", "grey", "grey", "red", "red", "red"))(100))

pdf("graphics/daf16_bubbleplot.pdf",
    height = 2.5,
    width = 2.5)

ggplot(aes(x = mean_z, y = p_val, size = 2^p_val, label = label), data = daf16_forbubbleplot) + 
  geom_point(col = daf16colours_for_plot, alpha = 0.7) + 
  theme_classic() + 
  geom_label_repel(size = 1.8, max.overlaps = 40, label.padding = 0.1) +
  geom_vline(xintercept = 0, linetype = "dashed", col = "grey") + 
  theme(legend.position = "none",
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10)) + 
  ylab(substitute("-log"[10]~"(mean p-value)")) + 
  xlab("TF activity (mean z-score)") 

dev.off()

#### DAF-18 ####

daf18_forbubbleplot <- data.frame(mean_z = decouplelist2[["daf18"]][["decouple_z"]][[1]],
                                 p_val = -log10(decouplelist2[["daf18"]][["decouple_p"]][[1]]),
                                 label = toupper(fullset_TFs_BM[match(row.names(decouplelist[["daf16"]][["decouple_p"]]), fullset_TFs_BM$wormbase_gseq), "label"]))

daf18_forbubbleplot[abs(daf18_forbubbleplot$p_val) < 2.2, "label"] <- ""

saveRDS(daf18_forbubbleplot,
        "plotdata/daf18_forbubbleplot.rds")

# daf18_forbubbleplot <- readRDS("plotdata/daf18_forbubbleplot.rds")

daf18colours_for_plot <- map2color(daf18_forbubbleplot$p_val, pal = colorRampPalette(c("grey", "grey", "grey", "red", "red", "red"))(100))


pdf("graphics/daf18_bubbleplot.pdf",
    height = 2.5,
    width = 2.5)

ggplot(aes(x = mean_z, y = p_val, size = 2^p_val, label = label), data = daf18_forbubbleplot) + 
  geom_point(col = daf18colours_for_plot, alpha = 0.7) + 
  theme_classic() + 
  geom_label_repel(size = 1.8, max.overlaps = 40, label.padding = 0.1) +
  geom_vline(xintercept = 0, linetype = "dashed", col = "grey") + 
  theme(legend.position = "none",
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10)) + 
  ylab(substitute("-log"[10]~"(mean p-value)")) + 
  xlab("TF activity (mean z-score)") 

dev.off()

#### DAF-2 DAF-18 ####

daf2daf18_forbubbleplot <- data.frame(mean_z = decouplelist2[["daf2daf18"]][["decouple_z"]][[1]],
                                  p_val = -log10(decouplelist2[["daf2daf18"]][["decouple_p"]][[1]]),
                                  label = toupper(fullset_TFs_BM[match(row.names(decouplelist[["daf16"]][["decouple_p"]]), fullset_TFs_BM$wormbase_gseq), "label"]))

daf2daf18_forbubbleplot[abs(daf2daf18_forbubbleplot$p_val) < 3, "label"] <- ""

saveRDS(daf2daf18_forbubbleplot,
        "plotdata/daf2daf18_forbubbleplot.rds")

# daf2daf18_forbubbleplot <- readRDS("plotdata/daf2daf18_forbubbleplot.rds")

daf2daf18colours_for_plot <- map2color(daf2daf18_forbubbleplot$p_val, pal = colorRampPalette(c("grey", "grey", "grey", "red", "red", "red"))(100))

pdf("graphics/daf2daf18_bubbleplot.pdf",
    height = 2.5,
    width = 2.5)

ggplot(aes(x = mean_z, y = p_val, size = 2^p_val, label = label), data = daf2daf18_forbubbleplot) + 
  geom_point(col = daf2daf18colours_for_plot, alpha = 0.7) + 
  theme_classic() + 
  geom_label_repel(size = 1.8, max.overlaps = 40, label.padding = 0.1) +
  geom_vline(xintercept = 0, linetype = "dashed", col = "grey") + 
  theme(legend.position = "none",
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        axis.title.x = element_text(size = 10),
        axis.title.y = element_text(size = 10)) + 
  ylab(substitute("-log"[10]~"(p-value)")) + 
  xlab("TF activity (z-score)") 

dev.off()
