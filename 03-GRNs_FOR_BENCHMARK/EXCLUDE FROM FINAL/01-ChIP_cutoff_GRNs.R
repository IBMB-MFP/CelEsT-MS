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

packages <- c()

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

#### ChIP-based GRNs ####

cutoffs_vec <- c(100, 500, 1000, 1500, 2000, 2500, 3000, 5000, 10000)

for(i in 1:length(cutoffs_vec)){
  
  message(paste0("Now doing cutoff ", cutoffs_vec[i]))
  make.ChIP.GRN(cutoff = cutoffs_vec[i],
                targets = allTFS_manualtoptargets,
                prefix = "allChIP",
                suffix = "HOTincl",
                min_targets = 15)
  
}

for(i in 1:length(cutoffs_vec)){
  
  message(paste0("Now doing cutoff ", cutoffs_vec[i]))
  make.ChIP.GRN(cutoff = cutoffs_vec[i],
                targets = allTFS_manualtoptargets_HOTexcl,
                prefix = "allChIP",
                suffix = "HOTexcl",
                min_targets = 15)
  
}

for(i in 1:length(cutoffs_vec)){
  
  message(paste0("Now doing cutoff ", cutoffs_vec[i]))
  make.ChIP.GRN(cutoff = cutoffs_vec[i],
                targets = allTFS_manualtoptargets_HOTexcl,
                prefix = "allChIP",
                suffix = "HOTexcl",
                min_targets = NULL)
  
}

#### MAKE HOT REGION CUTOFF GRNS ####

manualtargets_list <- list.files("output")

hotcutfiles <- manualtargets_list[str_detect(manualtargets_list, "HOTcut")]

lapply(hotcutfiles, function(thisfile){
  
  temptargets <- readRDS(paste0("output/", thisfile))
  
  Hotcut <- str_extract(thisfile, "HOTcut[0-9]{2}")
  
  make.ChIP.GRN(cutoff = 10000,
                targets = temptargets,
                prefix = "allChIP",
                suffix = Hotcut,
                min_targets = 15)
  
  NULL
  
})
