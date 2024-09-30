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

combine.GRNs <- function(GRNlist,
                         weighted = TRUE,
                         weights = NULL,
                         filename = NULL,
                         min_targets = 15){

  if(is.null(filename)){message("Error: you must suply a filename")
                               stop()}
  
  if(!is.null(weights) & length(GRNlist) != length(weights)){message("Error: weight must be a vector of length equal to the input list")
    stop()}

  combinedGRN_temp <- do.call(rbind, GRNlist)
  
  temp_unique_targets <- unique(combinedGRN_temp$target)[order(unique(combinedGRN_temp$target))]
  temp_unique_TFs <- unique(combinedGRN_temp$source)[order(unique(combinedGRN_temp$source))]
  
  tempGRNcombo_list <- lapply(GRNlist, function(thisGRN){
    
    templist <- sapply(temp_unique_TFs, function(thissource){
      
      thisone_targets <- thisGRN[thisGRN$source == thissource, "target"]
      
      tempout <- temp_unique_targets %in% thisone_targets
      
      names(tempout) <- temp_unique_targets
      
      tempout
      
    })
    
    t(templist)
    
  })
  
  
  if(weighted == TRUE & is.null(weights)){
    
    tempGRNcombo_array <- array(unlist(tempGRNcombo_list), dim = c(length(temp_unique_TFs), length(temp_unique_targets), length(tempGRNcombo_list)))
    
    dimnames(tempGRNcombo_array) <- list(temp_unique_TFs, temp_unique_targets, names(GRNlist))
    
    tempGRNcombo_df <- reshape2::melt(dimSums(tempGRNcombo_array, 3))
    
    colnames(tempGRNcombo_df) <- c("source", "target", "weight")
    
    tempGRNcombo_df <- tempGRNcombo_df[!tempGRNcombo_df$weight == 0, ]
    
    for(i in 1:length(GRNlist)){
    tempGRNcombo_df[tempGRNcombo_df$weight == i, "weight"] <- i / length(GRNlist)
    }
    
    tempGRNcombo_df <- tempGRNcombo_df[!tempGRNcombo_df$source %in% names(table(tempGRNcombo_df$source))[table(tempGRNcombo_df$source) < min_targets], ]
    
    write.table(tempGRNcombo_df,
                file = paste0("output/GRNs/", filename, ".txt"),
                row.names = FALSE,
                sep = "\t")
    
  }
  
  if(weighted == TRUE & !is.null(weights)){
    
    tempGRNcombo_list_weight <- lapply(1:length(tempGRNcombo_list), function(i){
      
      tempGRNcombo_list[[i]] * weights[i]
      
    })
    
    names(tempGRNcombo_list_weight) <- names(tempGRNcombo_list)
    
    tempGRNcombo_array <- array(unlist(tempGRNcombo_list_weight), dim = c(length(temp_unique_TFs), length(temp_unique_targets), length(tempGRNcombo_list_weight)))
    
    dimnames(tempGRNcombo_array) <- list(temp_unique_TFs, temp_unique_targets, names(GRNlist))
    
    tempGRNcombo_df <- reshape2::melt(dimSums(tempGRNcombo_array, 3))
    
    colnames(tempGRNcombo_df) <- c("source", "target", "weight")
    
    tempGRNcombo_df <- tempGRNcombo_df[!tempGRNcombo_df$weight == 0, ]
    
    tempGRNcombo_df <- tempGRNcombo_df[!tempGRNcombo_df$source %in% names(table(tempGRNcombo_df$source))[table(tempGRNcombo_df$source) < min_targets], ]
    
    write.table(tempGRNcombo_df,
                file = paste0("output/GRNs/", filename, ".txt"),
                row.names = FALSE,
                sep = "\t")
    
  }
  
  if(weighted == FALSE){
    
    tempGRNcombo_array <- array(unlist(tempGRNcombo_list), dim = c(length(temp_unique_TFs), length(temp_unique_targets), length(tempGRNcombo_list)))
    
    dimnames(tempGRNcombo_array) <- list(temp_unique_TFs, temp_unique_targets, names(GRNlist))
    
    tempGRNcombo_df <- reshape2::melt(dimSums(tempGRNcombo_array, 3))
    
    colnames(tempGRNcombo_df) <- c("source", "target", "weight")
    
    tempGRNcombo_df <- tempGRNcombo_df[!tempGRNcombo_df$weight == 0, ]
    
    tempGRNcombo_df[, "weight"] <- 1
    
    tempGRNcombo_df <- tempGRNcombo_df[!tempGRNcombo_df$source %in% names(table(tempGRNcombo_df$source))[table(tempGRNcombo_df$source) < min_targets], ]
    
    write.table(tempGRNcombo_df,
                file = paste0("output/GRNs/", filename, ".txt"),
                row.names = FALSE,
                sep = "\t")
    
  }
  
}

do.all.GRN.combinations <- function(GRNlist,
                                    prefix,
                                    weight_vec){
  
  combine.GRNs(GRNlist = GRNlist,
               weighted = FALSE,
               filename = paste0(prefix, "_noweights"))
  
  combine.GRNs(GRNlist = GRNlist,
               weighted = TRUE,
               weights = NULL,
               filename = paste0(prefix, "_equalweights"))
  
  combine.GRNs(GRNlist = GRNlist,
               weights = weight_vec,
               filename = paste0(prefix, "_weighted_unequal"))
  
}

#### LOAD PACKAGES & FUNCTIONS ####

## First specify the packages of interest

packages <- c("CHNOSZ")

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

github_packages <- c("r-lib/conflicted",
                     "LBMC/wormRef")

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
# allmodERN_TFsONLY_BM <- readRDS("~/Cel_GRN_manuscript/output/allmodERN_TFsONLY_BM.rds")


CisBP_TFinfo_withmotif <- readRDS("output/CisBP_TFinfo_withmotif.rds")
# CisBP_TFinfo_withmotif <- readRDS("~/Cel_GRN_manuscript/output/CisBP_TFinfo_withmotif.rds")

observations <- read.table("output/benchmark_observations.txt",
                           sep = "\t",
                           header = TRUE)

FIMO_nohomo_1000 <- read.table("output/GRNs/FIMO_nohomo_1000_unfiltered.txt",
                               sep = "\t",
                               header = TRUE)

fullsetTFs_BM <- readRDS("output/fullset_TFs_BM.rds")
# fullsetTFs_BM <- readRDS("~/Cel_GRN_manuscript/output/fullset_TFs_BM.rds")

# load unfiltered GRNs and filter after combination

ChIP_HOTexcl_nocut <- read.table("output/GRNs/allChIP_10000_HOTexcl_unfiltered.txt",
                                 sep = "\t",
                                 header = TRUE)
# ChIP_HOTexcl_nocut <- read.table("~/Cel_GRN_manuscript/output/GRNs/allChIP_10000_HOTexcl_unfiltered.txt",
#                                 sep = "\t",
#                                 header = TRUE)

eY1H_net <- read.table("output/GRNs/walhoutGRN_unfiltered.txt",
                       sep = "\t",
                       header = TRUE)
# eY1H_net <- read.table("~/Cel_GRN_manuscript/output/GRNs/walhoutGRN_unfiltered.txt",
#                        sep = "\t",
#                        header = TRUE)

#### COMBINE NETWORKS ####

# combine to form simple network

do.all.GRN.combinations(GRNlist = list(FIMO_nohomo_1000,
                                       ChIP_HOTexcl_nocut,
                                       eY1H_net),
                        prefix = "allthree",
                        weight_vec = c(0.2, 0.5, 0.3))

# do FIMO and MODern overlap restricted set with restricted observations for direct comparison

genes_present_in_ChIP_and_FIMO <- base::intersect(base::intersect(unique(FIMO_nohomo_1000$source), unique(ChIP_HOTexcl_nocut$source)), observations$target_gseq)

write.table(observations[observations$target_gseq %in% genes_present_in_ChIP_and_FIMO, ], 
            "output/FIMO1000ChIPnocutexcl_observations.txt",
            col.names = TRUE,
            row.names = TRUE,
            sep = "\t")

# combine FIMO and ChIP network with weights

combine.GRNs(GRNlist = list(FIMO_nohomo_1000,
                            ChIP_HOTexcl_nocut),
             weighted = TRUE,
             filename = "FIMO1000ChIPnocutexcl_weighted")

combine.GRNs(GRNlist = list(FIMO_nohomo_1000,
                            ChIP_HOTexcl_nocut),
             weighted = FALSE,
             filename = "FIMO1000ChIPnocutexcl_noweights")

#### ANNOTATE NETWORKS FOR SUPP TABLES ####

allthree_equalweights <- read.table("output/GRNs/allthree_equalweights.txt",
                                    sep = "\t",
                                    header = TRUE)

allthree_equalweights[, "with_motif"] <- NA
allthree_equalweights[allthree_equalweights$source %in% FIMO_nohomo_1000$source, "with_motif"] <- FALSE
allthree_equalweights[paste0(allthree_equalweights$source, "-", allthree_equalweights$target) %in% paste0(FIMO_nohomo_1000$source, "-", FIMO_nohomo_1000$target), "with_motif"] <- TRUE

allthree_equalweights[, "in_ChIP"] <- NA
allthree_equalweights[allthree_equalweights$source %in% ChIP_HOTexcl_nocut$source, "in_ChIP"] <- FALSE
allthree_equalweights[paste0(allthree_equalweights$source, "-", allthree_equalweights$target) %in% paste0(ChIP_HOTexcl_nocut$source, "-", ChIP_HOTexcl_nocut$target), "in_ChIP"] <- TRUE

allthree_equalweights[, "in_eY1H"] <- NA
allthree_equalweights[allthree_equalweights$source %in% eY1H_net$source, "in_eY1H"] <- FALSE
allthree_equalweights[paste0(allthree_equalweights$source, "-", allthree_equalweights$target) %in% paste0(eY1H_net$source, "-", eY1H_net$target), "in_eY1H"] <- TRUE

allthree_equalweights[, "TF_WBgeneID"] <- fullsetTFs_BM[match(allthree_equalweights$source, fullsetTFs_BM$wormbase_gseq), "wormbase_gene"]

allthree_equalweights[, "TF_name"] <- fullsetTFs_BM[match(allthree_equalweights$source, fullsetTFs_BM$wormbase_gseq), "wormbase_locus"]
allthree_equalweights[allthree_equalweights$TF_name == "", "TF_name"] <- allthree_equalweights[allthree_equalweights$TF_name == "", "source"]

allthree_equalweights[, "target_WBgeneID"] <- wormRef::Cel_genes[match(allthree_equalweights$target, wormRef::Cel_genes$sequence_name), "wb_id"]

allthree_equalweights[, "target_name"] <- wormRef::Cel_genes[match(allthree_equalweights$target, wormRef::Cel_genes$sequence_name), "public_name"]

write.table(allthree_equalweights,
            "output/CelEsT_annotated.txt",
            sep = "\t",
            col.names = TRUE,
            row.names = FALSE)
