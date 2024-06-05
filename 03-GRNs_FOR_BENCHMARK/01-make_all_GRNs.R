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

biocmanager_packages <- c("decoupleR") 

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
fullset_TFs_BM <- readRDS("output/fullset_TFs_BM.rds")

separate_by_motif <- readRDS("output/separate_by_motif.rds")

allTFS_manualtoptargets <- readRDS("output/MODernENCODE_manualtoptargets_operonexcluded.rds")
allTFS_manualtoptargets_HOTexcl <- readRDS("output/MODernENCODE_manualtoptargets_operon&HOTexcluded.rds")

walhoutEV1 <- openxlsx::read.xlsx("~/Downloads/msb167131-sup-0002-datasetev1.xlsx")

# limit to interactions described as high quality
walhoutEV1 <- walhoutEV1[walhoutEV1$`In.high-quality.dataset?` == 'yes', ]

parasite_mart <- useMart("parasite_mart", dataset = "wbps_gene", host = "https://parasite.wormbase.org", port = 443)

wTF3 <- read.xlsx("input/WTF3.xlsx")

# read in SJARACNe output file
SJARACNeage_out <- read.table("output/SJARACNe/age_corrected_non_zero/consensus_network_ncol_.txt",
                              header = TRUE)

CisBP_TFinfo_withmotif <- readRDS("output/CisBP_TFinfo_withmotif.rds")

observations <- read.table("output/benchmark_observations.txt",
                           sep = "\t",
                           header = TRUE)

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

# do unfiltered networks for combination

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

#### FIMO-based GRNs ####

cutoffs_vec <- c(100, 500, 1000, 1500, 2000, 2500, 3000, 5000, 10000)

FIMO_nohomo_calcs <- calc.score.by.motif(motif_list = separate_by_motif,
                                         homotypic_binding = FALSE)

saveRDS(FIMO_nohomo_calcs,
        "output/FIMO_nohomo_calcs.rds")

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

#### eY1H GRN ####

walhoutTFs <- unique(walhoutEV1$Prey.sequence.Name)

walhoutTFs_cutoff15  <- names(table(walhoutEV1$Prey.sequence.Name)[table(walhoutEV1$Prey.sequence.Name) >= 15])

# do we have any that are not in the wTF3.0 resource? NO of course 
any(!walhoutTFs %in% wTF3$Sequence.name)

walhoutGRN <- walhoutEV1[walhoutEV1$Prey.sequence.Name %in% walhoutTFs_cutoff15, c("Prey.sequence.Name", "Bait.sequence.name")]
colnames(walhoutGRN) <- c("source", "target")

walhoutGRN <- walhoutGRN[!duplicated(paste(walhoutGRN$source, walhoutGRN$target)), ]
walhoutGRN[, "weight"] <- 1

write.table(walhoutGRN,
            "output/GRNs/walhout_highqual_cutoff15.txt",
            col.names = TRUE,
            row.names = FALSE,
            sep = "\t")

#### SJARACNE GRNs ####

# filter by p < 1e-5

SJARACNEage_out_filt <- SJARACNeage_out[SJARACNeage_out$p.value < 0.00001, ]

SJARACNEage_out_filt_separated <- lapply(unique(SJARACNEage_out_filt$source), function(x){
  
  thisone_SJ <- SJARACNEage_out_filt[SJARACNEage_out_filt$source == x, ]
  
  thisone_SJ[order(thisone_SJ$MI, decreasing = TRUE), ]
  
})

names(SJARACNEage_out_filt_separated) <- unique(SJARACNEage_out_filt$source)

SJARACNEage_GRN_list <- lapply(SJARACNEage_out_filt_separated, function(x){
  
  temp_df <- x[, c("source", "target")]
  
  temp_df[, "weight"] <- sign(x$pearson)
  
  temp_df
  
})

SJARACNEage_GRN <- do.call(rbind, SJARACNEage_GRN_list)

write.table(SJARACNEage_GRN,
            "output/GRNs/SJARACNEage_GRN.txt",
            col.names = TRUE,
            row.names = FALSE,
            sep = "\t")



