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

setwd("~/Cel_GRN_manuscript")

#### DEFINE FUNCTIONS ####

ecdf_fun <- function(x,perc) ecdf(x)(perc)

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


#### LOAD PACKAGES & FUNCTIONS ####

## First specify the packages of interest

packages <- c("stringr",
              "ggplot2",
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

STREMEtopmotif_orthology_probs <- readRDS("output/STREMEtopmotif_orthology_probs10k.rds")

allTFS_manualtoptargets_HOTexcl <- readRDS(file = "output/MODernENCODE_manualtoptargets_operon&HOTexcluded.rds")
allTFS_manualtoptargets <- readRDS(file = "output/MODernENCODE_manualtoptargets_operonexcluded.rds")

allTFs_labelused <- readRDS("output/allTFs_labelused.rds")

allmodERN_TFsONLY_BM <- readRDS("output/allmodERN_TFsONLY_BM.rds")

TFs <- names(allTFs_labelused)[match(allTFs_labelused[allTFs_labelused != "cfi-1"], names(allTFS_manualtoptargets_HOTexcl))]

bench_out_filelist <- list.files("output/benchmark_out/")

TF_orthology_probs <- readRDS("output/TF_orthology_probs10k.rds")

#### ORDER CHIP TARGETS BY CONSERVATION ####

neworder_targets_Hotexcl <- lapply(TFs, function(x){
print(x)
  temp <- allTFS_manualtoptargets_HOTexcl[[allTFs_labelused[x]]]
  
  if(!is.null(STREMEtopmotif_orthology_probs[[x]][match(temp, names(STREMEtopmotif_orthology_probs[[x]]))])){
  
  tempout <- data.frame("target" = temp, "probs" = 1) 
  tempout[, "probs"] <- STREMEtopmotif_orthology_probs[[x]][match(temp, names(STREMEtopmotif_orthology_probs[[x]]))]
  
  tempout[is.na(tempout$probs), "probs"] <- 1
  
  tempout <- tempout[order(tempout$probs), ]
  
  # order scores of 1 according to ChIP score. although it seems this order is retained, anyway
  
  tempout[tempout$probs == 1, "target"] <- temp[temp %in% tempout[tempout$probs == 1, "target"]]
  
  return(tempout$target)
  
} else {
    
  return(temp)

  }
  
})

names(neworder_targets_Hotexcl) <- TFs

neworder_targets <- lapply(TFs, function(x){
  print(x)
  temp <- allTFS_manualtoptargets[[allTFs_labelused[x]]]
  
  if(!is.null(STREMEtopmotif_orthology_probs[[x]][match(temp, names(STREMEtopmotif_orthology_probs[[x]]))])){
    
    tempout <- data.frame("target" = temp, "probs" = 1) 
    tempout[, "probs"] <- STREMEtopmotif_orthology_probs[[x]][match(temp, names(STREMEtopmotif_orthology_probs[[x]]))]
    
    tempout[is.na(tempout$probs), "probs"] <- 1
    
    tempout <- tempout[order(tempout$probs), ]
    
    # order scores of 1 according to ChIP score. although it seems this order is retained, anyway
    
    tempout[tempout$probs == 1, "target"] <- temp[temp %in% tempout[tempout$probs == 1, "target"]]
    
    return(tempout$target)
    
  } else {
    
    return(temp)
    
  }
  
})

names(neworder_targets) <- TFs

#### WRITE GRNS ####

cutoffs_vec <- c(100, 500, 1000, 1500, 2000, 2500, 3000, 5000, 10000)

for(i in 1:length(cutoffs_vec)){
  
  message(paste0("Now doing cutoff ", cutoffs_vec[i]))
  make.ChIP.GRN(cutoff = cutoffs_vec[i],
                targets = neworder_targets,
                prefix = "ChIPorth",
                suffix = "HOTincl",
                min_targets = 15)
  
}

for(i in 1:length(cutoffs_vec)){

  message(paste0("Now doing cutoff ", cutoffs_vec[i]))
  make.ChIP.GRN(cutoff = cutoffs_vec[i],
                targets = neworder_targets_Hotexcl,
                prefix = "ChIPorth",
                suffix = "HOTexcl",
                min_targets = 15)
  
}

## make unfiltered for HOT exclusion with cutoff 1000, since through magical foresight I know this is the one that will come out ahead

make.ChIP.GRN(cutoff = 1000,
              targets = neworder_targets_Hotexcl,
              prefix = "ChIPorth",
              suffix = "HOTexcl",
              min_targets = NULL)

#### TRY WITH CISBP PROBS INSTEAD OF STREME FOR THOSE ONES WITH CISBP MOTIFS ####

STREMECisBP_orthology_probs <- STREMEtopmotif_orthology_probs

shared_TFs <- intersect(names(TF_orthology_probs),
                        names(STREMECisBP_orthology_probs))

STREMECisBP_orthology_probs[shared_TFs] <- TF_orthology_probs[shared_TFs]

STREMECisBPorder_targets_Hotexcl <- lapply(TFs, function(x){
  print(x)
  temp <- allTFS_manualtoptargets_HOTexcl[[allTFs_labelused[x]]]
  
  if(!is.null(STREMECisBP_orthology_probs[[x]][match(temp, names(STREMECisBP_orthology_probs[[x]]))])){
    
    tempout <- data.frame("target" = temp, "probs" = 1) 
    tempout[, "probs"] <- STREMECisBP_orthology_probs[[x]][match(temp, names(STREMECisBP_orthology_probs[[x]]))]
    
    tempout[is.na(tempout$probs), "probs"] <- 1
    
    tempout <- tempout[order(tempout$probs), ]
    
    # order scores of 1 according to ChIP score. although it seems this order is retained, anyway
    
    tempout[tempout$probs == 1, "target"] <- temp[temp %in% tempout[tempout$probs == 1, "target"]]
    
    return(tempout$target)
    
  } else {
    
    return(temp)
    
  }
  
})

names(STREMECisBPorder_targets_Hotexcl) <- TFs

for(i in 1:length(cutoffs_vec)){
  
  message(paste0("Now doing cutoff ", cutoffs_vec[i]))
  make.ChIP.GRN(cutoff = cutoffs_vec[i],
                targets = STREMECisBPorder_targets_Hotexcl,
                prefix = "STREMECisBP",
                suffix = "HOTexcl",
                min_targets = 15)
  
}



STREMECisBPorder_targets_Hotincl <- lapply(TFs, function(x){
  print(x)
  temp <- allTFS_manualtoptargets[[allTFs_labelused[x]]]
  
  if(!is.null(STREMECisBP_orthology_probs[[x]][match(temp, names(STREMECisBP_orthology_probs[[x]]))])){
    
    tempout <- data.frame("target" = temp, "probs" = 1) 
    tempout[, "probs"] <- STREMECisBP_orthology_probs[[x]][match(temp, names(STREMECisBP_orthology_probs[[x]]))]
    
    tempout[is.na(tempout$probs), "probs"] <- 1
    
    tempout <- tempout[order(tempout$probs), ]
    
    # order scores of 1 according to ChIP score. although it seems this order is retained, anyway
    
    tempout[tempout$probs == 1, "target"] <- temp[temp %in% tempout[tempout$probs == 1, "target"]]
    
    return(tempout$target)
    
  } else {
    
    return(temp)
    
  }
  
})

names(STREMECisBPorder_targets_Hotincl) <- TFs

for(i in 1:length(cutoffs_vec)){
  
  message(paste0("Now doing cutoff ", cutoffs_vec[i]))
  make.ChIP.GRN(cutoff = cutoffs_vec[i],
                targets = STREMECisBPorder_targets_Hotincl,
                prefix = "STREMECisBP",
                suffix = "HOTincl",
                min_targets = 15)
  
}
