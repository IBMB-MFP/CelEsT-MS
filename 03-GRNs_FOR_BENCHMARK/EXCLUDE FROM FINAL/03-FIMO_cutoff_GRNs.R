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

packages <- c("stringr")

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

separate_by_motif <- readRDS("output/separate_by_motif.rds")

CisBP_TFinfo_withmotif <- readRDS("output/CisBP_TFinfo_withmotif.rds")

observations <- read.table("output/benchmark_observations.txt",
            sep = "\t",
            header = TRUE)

#### DEFINE FUNCTIONS ####

calc.score.by.motif <- function(motif_list = separate_by_motif,
                                homotypic_binding = TRUE,
                                returns_parameter = 3) {
  
  
  outlist <- lapply(motif_list, function(x){
    
    sapply(unique(x$sequence_name), function(y){
      
      thisonedata <- x[x$sequence_name == y, ]
      
      n <- nrow(thisonedata)
      
      if(n == 1){
        
        return(thisonedata$score)
        
      } else {
        
        if(homotypic_binding == TRUE){
          
          thisonedata <- thisonedata[order(thisonedata$score, decreasing = TRUE), ]
          
          sumscore <- 0
          
          for(i in 1:n){
            
            thisscore <- thisonedata[i, 'score']
            sumscore <- sumscore + thisscore * (1 / returns_parameter^(i-1))
            
          }
          
          return(sumscore)
          
        } else {
          
          return(thisonedata$score[1])
          
        }
        
      }
      
    })
    
  })
  
  outlist_order <- lapply(outlist, function(x){
    
    x[order(x, decreasing = TRUE)]
    
  })
  
  return(outlist_order)
  
}

#### FIMO-based GRNs ####

cutoffs_vec <- c(100, 500, 1000, 1500, 2000, 2500, 3000, 5000, 10000)

FIMO_nohomo_calcs <- calc.score.by.motif(motif_list = separate_by_motif,
                                         homotypic_binding = FALSE)

FIMO_homo_list <- lapply(c(1, 3, 5, 7), function(x){calc.score.by.motif(returns_parameter = x)})
names(FIMO_homo_list) <- c(paste0("param", c(1, 3, 5, 7)))

# write GRNs for no homo

# identify duplicates (kept when both in modERN... )
# in fact they are better not kept; it impacts performance to have them in the final GRN

FIMOduplicated_TFs <- CisBP_TFinfo_withmotif[CisBP_TFinfo_withmotif$wormbase_gseq %in% names(separate_by_motif), ]$wormbase_gseq[(duplicated(CisBP_TFinfo_withmotif[CisBP_TFinfo_withmotif$wormbase_gseq %in% names(separate_by_motif), ]$Motif_ID)|duplicated(CisBP_TFinfo_withmotif[CisBP_TFinfo_withmotif$wormbase_gseq %in% names(separate_by_motif), ]$Motif_ID, fromLast = TRUE))]

# any duplicates in my benchmarking set? NO
# so eliminate them all for this

for(i in 1:length(cutoffs_vec)){
  
  message(paste0("Now doing cutoff ", cutoffs_vec[i]))

  temp_calcscore_GRN <- reshape2::melt(lapply(FIMO_nohomo_calcs, function(x){names(x)[1:min(cutoffs_vec[i], length(x))]}))
  colnames(temp_calcscore_GRN) <- c("target", "source")
  temp_calcscore_GRN[, "weight"] <- 1
  
  temp_calcscore_GRN <- temp_calcscore_GRN[!temp_calcscore_GRN$source %in% FIMOduplicated_TFs, ]
  
  colinearity_check <- check_corr(temp_calcscore_GRN,
             .source = "source",
             .target = "target",
             .mor = "weight")
  
  colinear_TFs <- colinearity_check[colinearity_check$correlation > 0.99, c("source", "source.2")]
  
  if(any(unlist(colinear_TFs) %in% observations$target_gseq)){
    
    exclude = apply(colinear_TFs, 1, function(z){
      
      if(sum(z %in% observations$target_gseq) == 1){
        
        exclude = z[!z %in% observations$target_gseq]
        
        return(exclude)
        
      }
      
      if(sum(z %in% observations$target_gseq) == 0){
      
      exclude = z[1]
      
      return(exclude)
      
      }
      
      if(sum(z %in% observations$target_gseq == 2)){
        
        message(paste0("Duplicated TFs both in benchmarkingset: ", z[1], " and ", z[2]))
        
      }
      
    })
    
  }

  write.table(temp_calcscore_GRN[!temp_calcscore_GRN$source %in% exclude, ],
              file = paste0("output/GRNs/FIMO_nohomo_", cutoffs_vec[i], ".txt"),
              sep = "\t",
              col.names = TRUE,
              row.names = FALSE)
  
}

# for cutoff 1500 (final choice) output unfiltered GRN (solely for combination with other)

unfilt1500_calcscore_GRN <- reshape2::melt(lapply(FIMO_nohomo_calcs, function(x){names(x)[1:min(1500, length(x))]}))
colnames(unfilt1500_calcscore_GRN) <- c("target", "source")
unfilt1500_calcscore_GRN[, "weight"] <- 1

write.table(unfilt1500_calcscore_GRN,
            file = "output/GRNs/FIMO_nohomo_1500_unfiltered.txt",
            sep = "\t",
            col.names = TRUE,
            row.names = FALSE)

# homotypic binding GRNs

lapply(names(FIMO_homo_list), function(thisname){
  
  for(i in 1:length(cutoffs_vec)){
  
  message(paste0("Now doing ", thisname, ", cutoff ", cutoffs_vec[i]))
  
  thesescores <- FIMO_homo_list[[thisname]]
  
  temp_calcscore_GRN <- reshape2::melt(lapply(thesescores, function(x){names(x)[1:min(cutoffs_vec[i], length(x))]}))
  colnames(temp_calcscore_GRN) <- c("target", "source")
  temp_calcscore_GRN[, "weight"] <- 1
  
  temp_calcscore_GRN <- temp_calcscore_GRN[!temp_calcscore_GRN$source %in% FIMOduplicated_TFs, ]
  
  colinearity_check <- check_corr(temp_calcscore_GRN,
                                  .source = "source",
                                  .target = "target",
                                  .mor = "weight")
  
  colinear_TFs <- colinearity_check[colinearity_check$correlation > 0.99, c("source", "source.2")]
  
  if(any(unlist(colinear_TFs) %in% observations$target_gseq)){
    
    exclude = apply(colinear_TFs, 1, function(z){
      
      if(sum(z %in% observations$target_gseq) == 1){
        
        exclude = z[!z %in% observations$target_gseq]
        
        return(exclude)
        
      }
      
      if(sum(z %in% observations$target_gseq) == 0){
        
        exclude = z[1]
        
        return(exclude)
        
      }
      
      if(sum(z %in% observations$target_gseq == 2)){
        
        message(paste0("Duplicated TFs both in benchmarkingset: ", z[1], " and ", z[2]))
        
      }
      
    })
    
  }
  
  write.table(temp_calcscore_GRN[!temp_calcscore_GRN$source %in% exclude, ],
              file = paste0("output/GRNs/FIMO_homo", thisname, "_", cutoffs_vec[i], ".txt"),
              sep = "\t",
              col.names = TRUE,
              row.names = FALSE)
  
  NULL
  
  }
  
})