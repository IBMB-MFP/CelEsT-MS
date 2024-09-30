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

setwd("~/Cel_GRN_revisions")

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
              "biomaRt",
              "readxl", # for reading xls files (NOTE not xlsx)
              "ggplot2",
              "ggrepel",
              "gplots") 

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
GRN_filelist <- list.files("output/GRNs") 

# GRN_filelist <- list.files("~/Cel_GRN_manuscript/output/GRNs") 

#### SCATTERPLOTS ####

TForthprobs_files <- bench_out_filelist[str_detect(bench_out_filelist, "TF_orth")]
TForthprobs_GRNs <- GRN_filelist[str_detect(GRN_filelist, "TF_orth")]

cutoffs_vec <- c(500, 1000, 1500, 2000, 2500, 3000, 5000, 10000)
fdr_vec <- c(0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9)

method_labels <- c("Consensus", "MLM", "ULM", "WSum")
names(method_labels) <- c("consensus_estimate", "mlm_estimate", "ulm_estimate", "wsum_norm")

TFnumbers <- sapply(cutoffs_vec, function(thiscutoff){

  sapply(fdr_vec, function(thisfdr){

    # tempGRN <- read.table(paste0("output/GRNs/TF_orthprobs_cut", thiscutoff, "_fdr", thisfdr, ".txt"),
    #                       header = TRUE)
    tempGRN <- read.table(paste0("~/Cel_GRN_manuscript/output/GRNs/TF_orthprobs_cut", thiscutoff, "_fdr", thisfdr, ".txt"),
                          header = TRUE)
    
    length(unique(tempGRN$source))
    
  })
  
})

colnames(TFnumbers) <- cutoffs_vec
row.names(TFnumbers) <- fdr_vec

TF_meantarget_numbers <- sapply(cutoffs_vec, function(thiscutoff){
  
  sapply(fdr_vec, function(thisfdr){
    # 
    # tempGRN <- read.table(paste0("output/GRNs/TF_orthprobs_cut", thiscutoff, "_fdr", thisfdr, ".txt"),
    #                       header = TRUE)

    tempGRN <- read.table(paste0("~/Cel_GRN_manuscript/output/GRNs/TF_orthprobs_cut", thiscutoff, "_fdr", thisfdr, ".txt"),
                          header = TRUE)
    
    mean(table(tempGRN$source))
    
  })
  
})

colnames(TF_meantarget_numbers) <- cutoffs_vec
row.names(TF_meantarget_numbers) <- fdr_vec

TF_mediantarget_numbers <- sapply(cutoffs_vec, function(thiscutoff){
  
  sapply(fdr_vec, function(thisfdr){
    
    # tempGRN <- read.table(paste0("output/GRNs/TF_orthprobs_cut", thiscutoff, "_fdr", thisfdr, ".txt"),
    #                       header = TRUE)
    
    tempGRN <- read.table(paste0("~/Cel_GRN_manuscript/output/GRNs/TF_orthprobs_cut", thiscutoff, "_fdr", thisfdr, ".txt"),
                          header = TRUE)
    
    median(table(tempGRN$source))
    
  })
  
})

colnames(TF_mediantarget_numbers) <- cutoffs_vec
row.names(TF_mediantarget_numbers) <- fdr_vec

for(i in 1:length(fdr_vec)){

TForthprobs_tempfdr_benchraptor_list <- lapply(TForthprobs_files[str_detect(TForthprobs_files, paste0("fdr", fdr_vec[i], "_bench"))], function(x){
         
read.table(paste0("output/benchmark_out/", x),
                    sep = "\t",
                    header = TRUE)
         
})
       
TForthprobs_tempfdr_benchraptor_all <- do.call(rbind, TForthprobs_tempfdr_benchraptor_list)
       
TForthprobs_tempfdr_benchraptor_allmethods_plot <- TForthprobs_tempfdr_benchraptor_all[TForthprobs_tempfdr_benchraptor_all$metric %in% c("auroc", "auprc"), c("net", "method", "metric", "score")]
TForthprobs_tempfdr_benchraptor_allmethods_plot <- tidyr::pivot_wider(TForthprobs_tempfdr_benchraptor_allmethods_plot, names_from = c("metric"), values_from = "score")
       
TForthprobs_tempfdr_benchraptor_allmethods_plot[, "label"] <- as.character(TForthprobs_tempfdr_benchraptor_allmethods_plot$net)

pdf(paste0("graphics/TForth_methodsplot_fdr", fdr_vec[i], ".pdf"),
    width = 5,
    height = 4)

print(
  ggplot(data = TForthprobs_tempfdr_benchraptor_allmethods_plot[!TForthprobs_tempfdr_benchraptor_allmethods_plot$method %in% c("wsum_estimate", "wsum_corr"), ], aes(x = auroc, y = auprc, colour = net)) + 
  geom_point(size = 2) +
  geom_label_repel(aes(label = net),
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
  facet_wrap(~ method, labeller = labeller(method = method_labels)) + 
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
        axis.text = element_text(colour = "black"))
)
  
dev.off()

}

for(i in 1:length(cutoffs_vec)){

  TForthprobs_tempcut_benchraptor_list <- lapply(TForthprobs_files[str_detect(TForthprobs_files, paste0("cut", cutoffs_vec[i], "_"))], function(x){
    
    read.table(paste0("output/benchmark_out/", x),
               sep = "\t",
               header = TRUE)
    
  })
  
  TForthprobs_tempcut_benchraptor_all <- do.call(rbind, TForthprobs_tempcut_benchraptor_list)
  
  TForthprobs_tempcut_benchraptor_allmethods_plot <- TForthprobs_tempcut_benchraptor_all[TForthprobs_tempcut_benchraptor_all$metric %in% c("auroc", "auprc"), c("net", "method", "metric", "score")]
  TForthprobs_tempcut_benchraptor_allmethods_plot <- tidyr::pivot_wider(TForthprobs_tempcut_benchraptor_allmethods_plot, names_from = c("metric"), values_from = "score")
  
  TForthprobs_tempcut_benchraptor_allmethods_plot[, "label"] <- as.character(TForthprobs_tempcut_benchraptor_allmethods_plot$net)
  
  pdf(paste0("graphics/TForth_methodsplot_fdr", cutoffs_vec[i], ".pdf"),
      width = 5,
      height = 4)
  
  print(
    ggplot(data = TForthprobs_tempcut_benchraptor_allmethods_plot[!TForthprobs_tempcut_benchraptor_allmethods_plot$method %in% c("wsum_estimate", "wsum_corr"), ], aes(x = auroc, y = auprc, colour = net)) + 
      geom_point(size = 2) +
      geom_label_repel(aes(label = net),
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
      facet_wrap(~ method, labeller = labeller(method = method_labels)) + 
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
            axis.text = element_text(colour = "black"))
  )
  
  dev.off()
  
}

#### HEATMAPS FOR CUTOFF / FDR ####

TF_auroc_mlm <- sapply(cutoffs_vec, function(thiscutoff){
  
  sapply(fdr_vec, function(thisfdr){

    tempbench <- read.table(paste0("output/benchmark_out/TF_orthprobs_cut", thiscutoff, "_fdr", thisfdr, "_benchRAPToR.tsv"),
                          sep = "\t",
                          header = TRUE)
    
    tempbench[tempbench$method == "mlm_estimate" & tempbench$metric == "auroc", "score"]
    
  })
  
})

colnames(TF_auroc_mlm) <- cutoffs_vec
row.names(TF_auroc_mlm) <- fdr_vec

TF_auprc_mlm <- sapply(cutoffs_vec, function(thiscutoff){
  
  sapply(fdr_vec, function(thisfdr){
    
    tempbench <- read.table(paste0("output/benchmark_out/TF_orthprobs_cut", thiscutoff, "_fdr", thisfdr, "_benchRAPToR.tsv"),
                            sep = "\t",
                            header = TRUE)
    
    tempbench[tempbench$method == "mlm_estimate" & tempbench$metric == "auprc", "score"]
    
  })
  
})

colnames(TF_auprc_mlm) <- cutoffs_vec
row.names(TF_auprc_mlm) <- fdr_vec

mainpal <- (colorRampPalette(c("red", "orange", "yellow", "white"))(100))
colors = seq(0.6, 0.7, length.out=101)

saveRDS(t(TF_auprc_mlm),
        "plotdata/TF_auprc_mlm.rds")

pdf("graphics/TForth_cutoffheatmap_auprc_mlm.pdf",
    width = 5,
    height = 5)

gplots::heatmap.2(x = t(TF_auprc_mlm),
                  dendrogram = "none",
                  trace = "none",
                  Rowv = FALSE,
                  Colv = FALSE,
                  density.info = "none",
                  col = mainpal,
                  breaks = colors,
                  key = FALSE,
                  xlab = "TF-target conservation FDR",
                  ylab = "Max. targets cutoff")

dev.off()

saveRDS(t(TF_auroc_mlm),
        "plotdata/TF_auroc_mlm.rds")

pdf("graphics/TForth_cutoffheatmap_auroc_mlm.pdf",
    width = 5,
    height = 5)

gplots::heatmap.2(x = t(TF_auroc_mlm),
                  dendrogram = "none",
                  trace = "none",
                  Rowv = FALSE,
                  Colv = FALSE,
                  density.info = "none",
                  col = mainpal,
                  breaks = colors,
                  key = FALSE)

dev.off()

pdf("graphics/TForth_cutoffheatmap_auprc_mlm_KEY.pdf",
    width = 5,
    height = 5)

gplots::heatmap.2(x = t(TF_auprc_mlm),
                  dendrogram = "none",
                  trace = "none",
                  Rowv = FALSE,
                  Colv = FALSE,
                  density.info = "none",
                  col = mainpal,
                  breaks = colors,
                  key = TRUE,
                  xlab = "TF-target conservation FDR",
                  ylab = "Max. targets cutoff")

dev.off()

