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

return.top.chip.peaks <- function(tf_name,
                                  GR_list = notHOTpeaks_by_TF_GRlist,
                                  cutoff = 500){

  thisTF_peaks_GR <-  GR_list[[tf_name]]
  
  # first the case in which there is only a single replicate
  if(length(unique(thisTF_peaks_GR$experiment)) == 1){
    
    if(length(thisTF_peaks_GR) > cutoff){
      
      top_targets <- thisTF_peaks_GR[order(thisTF_peaks_GR$score, decreasing = TRUE), ][1:cutoff]
      
    } else {
      
      top_targets <- thisTF_peaks_GR[order(thisTF_peaks_GR$score, decreasing = TRUE), ]
      
    }
    
  }
  
  if(length(unique(thisTF_peaks_GR$experiment)) > 1){
    
    thisTF_stages <- unique(thisTF_peaks_GR$experiment)
    
    thisTF_peaks <- lapply(thisTF_stages, function(thisstage){
      
      thisTF_peaks_GR[str_detect(thisTF_peaks_GR$experiment, thisstage), ]
      
    })
    
    # find overlaps between stages
    # find a way to do all combinations of stages
    # then we collapse it down to find the ranges with overlapping peaks
    # these will take priority in the cutoff and will be ordered amongst themselves by the maximum score observed
    
    combined_ranges <- c()
    
    for(i in 1:length(thisTF_peaks)){
      
      j_vec <- c(1:length(thisTF_peaks))[-i]
      
      joverlaps_list <- lapply(j_vec, function(j){
        
        tempoverlaps <- findOverlaps(thisTF_peaks[[i]], thisTF_peaks[[j]])
        c(thisTF_peaks[[i]][queryHits(tempoverlaps)], thisTF_peaks[[j]][subjectHits(tempoverlaps)])
        
      })
      
      joverlaps_vec <- do.call(c, joverlaps_list)
      
      combined_ranges <- c(combined_ranges, joverlaps_vec)
      
    }
    
    combined_ranges <- do.call(c, combined_ranges)
    
    mcols(combined_ranges)[, "number_of_stages"] <- countOverlaps(combined_ranges, combined_ranges)
    
    reduced_ranges <- reduce(combined_ranges)
    
    # now need to find a way to retain the counts
    mcols(reduced_ranges)[, "number_of_stages"] <- mcols(combined_ranges)[subjectHits(findOverlaps(reduced_ranges, combined_ranges))[!duplicated(queryHits(findOverlaps(reduced_ranges, combined_ranges)))], "number_of_stages"]
    
    # find top scores for reduced_ranges 
    scores_vec <- c()
    
    if(length(reduced_ranges) > 0){
      
      for(k in 1:length(reduced_ranges)){
        
        scores_vec[k] <- max(sapply(thisTF_peaks, function(thisstage){
          
          max(mcols(thisstage[subjectHits(findOverlaps(reduced_ranges[k], thisstage)), ])$score)
          
        }))
        
      }
      
    }
    
    mcols(reduced_ranges)[, "max_score"] <- scores_vec
    
    if(!is.null(reduced_ranges$max_score)){
      reduced_ranges <- reduced_ranges[order(reduced_ranges$max_score, decreasing = TRUE)]
      reduced_ranges <- reduced_ranges[order(reduced_ranges$number_of_stages, decreasing = TRUE)]
    }
    
    # now its time to combine with other ranges which didn't overlap
    # these will take lower priority than the overlappping ranges and will be ordered amongst themselves by score
    
    if(length(reduced_ranges) > 0){
      
      not_overlaps_list <- lapply(thisTF_peaks, function(thisnooverlap){
        
        thisnooverlap[-queryHits(findOverlaps(thisnooverlap, reduced_ranges)), ]
        
      })
      
    } else {
      
      not_overlaps_list <- thisTF_peaks
      
    }
    
    not_overlapping_GR <- do.call(c, not_overlaps_list)
    not_overlapping_GR <- not_overlapping_GR[order(not_overlapping_GR$score, decreasing = TRUE), ]
    
    all_ordered_GR <- c(reduced_ranges, not_overlapping_GR)
    
    if(length(all_ordered_GR) > cutoff){
      
      top_targets <- all_ordered_GR[1:cutoff, ]
      
    } else {
      
      top_targets <- all_ordered_GR
      
    }
    
  }
  
  # here convert top targets back into a GRanges object
  return(top_targets)
  
}

make.TF.orth.GRN <- function(cutoff, fdrcut){
  # cutoff = 1000
  # fdrcut = 0
  tempGRN <- do.call(rbind, lapply(names(TF_orthology_probs), function(thisTF){
    # thisTF <- names(TF_orthology_probs)[5]
    x <- TF_orthology_probs[[thisTF]]
    
    fdrcutx <- names(x[!x > fdrcut])
    
    elegans_x <- caeno_target_array[, , "elegans"][, thisTF]
    
    elegans_x_order <- elegans_x[order(elegans_x, decreasing = TRUE)]
    
    thiscutoff <- min(cutoff, sum(elegans_x_order != 0))
    
    fdrcut_names <- names(elegans_x_order[1:thiscutoff])[names(elegans_x_order[1:thiscutoff]) %in% fdrcutx]
    
    data.frame("source" = rep(thisTF, times = length(fdrcut_names)),
               "target" = fdrcut_names,
               "weight" = rep(1, times = length(fdrcut_names)))
    
  }))
  
  # let's exclude fewer than 10
  tempGRN <- tempGRN[!tempGRN$source %in% names(table(tempGRN$source))[table(tempGRN$source) < 10], ]
  
  write.table(tempGRN,
              file = paste0("output/GRNs/TF_MODernSTREMEmotiforthprobs_cut", cutoff, "_fdr", fdrcut, ".txt"),
              sep = "\t",
              row.names = FALSE,
              col.names = TRUE)
  
}

split.GR.byTF.tolist2 <- function(GR_object){
  # GR_object <- MODern_agg_GR
  all_in_peaks_files_nameused_TFonly <- readRDS("~/Cel_GRN_orthology/output/all_in_peaks_files_nameused_TFonly.rds")
  all_in_peaks_files_BM <- readRDS("~/Cel_GRN_orthology/output/all_in_peaks_files_BM.rds")
  
  templist <- lapply(all_in_peaks_files_nameused_TFonly, function(thisTF){
    
    GR_object[str_detect(GR_object$experiment, thisTF), ]
    
  })
  
  # could use all_in_peaks_files_gseq_TFonly but just to be sure will match to BM, in case e.g. all peaks from some TF fall into HOT regions etc.
  tempnames <- all_in_peaks_files_nameused_TFonly
  tempnames[which(!tempnames %in% all_in_peaks_files_BM$wormbase_gseq)] <- all_in_peaks_files_BM[match(tempnames[!tempnames %in% all_in_peaks_files_BM$wormbase_gseq], all_in_peaks_files_BM$wormbase_locus), "wormbase_gseq"]
  names(templist) <- tempnames
  
  templist
  
}

make.TF.ChIP.orth.GRN <- function(cutoff,
                                  fdrcut,
                                  orthology_probs = STREMEtopmotif_orthology_probs,
                                  prefix = "TF_MODernSTREMEorthprobs_cut",
                                  target_list = allTFS_manualtoptargets_HOTexcl,
                                  min_target_cut = 100){
  
  # cutoff = 2500
  # fdrcut = 0.5
  tempGRN <- do.call(rbind, lapply(names(orthology_probs), function(thisTF){
    # thisTF <- Y70C5C.1
    x <- orthology_probs[[thisTF]]
    
    fdrcutx <- names(x[!x > fdrcut])
    
    elegans_x <- allTFS_manualtoptargets_HOTexcl[[allTFs_labelused[thisTF]]]
    
    if(length(elegans_x) < min_target_cut){
      return(data.frame("source" = rep(thisTF, times = length(elegans_x)),
                        "target" = elegans_x,
                        "weight" = rep(1, times = length(elegans_x))))
    }
    
    thiscutoff <- min(cutoff, length(elegans_x))
    
    fdrcut_names <- elegans_x[1:thiscutoff][elegans_x[1:thiscutoff] %in% fdrcutx]
    
    data.frame("source" = rep(thisTF, times = length(fdrcut_names)),
               "target" = fdrcut_names,
               "weight" = rep(1, times = length(fdrcut_names)))
    
  }))

  # get rid of duplications
  tempGRN <- tempGRN[!duplicated(paste(tempGRN$source, tempGRN$target)), ]
  
  # let's exclude fewer than 10
  tempGRN <- tempGRN[!tempGRN$source %in% names(table(tempGRN$source))[table(tempGRN$source) < 10], ]
  
  write.table(tempGRN,
              file = paste0("output/GRNs/", prefix, cutoff, "_fdr", fdrcut, ".txt"),
              sep = "\t",
              row.names = FALSE,
              col.names = TRUE)
  
}

make.control.TF.ChIP.orth.GRN <- function(cutoff,
                                          fdrcut,
                                          orthology_probs = TF_orthology_probs,
                                          prefix = "TF_MODernSTREMEorthprobs_controlcut"){
  
  # cutoff = 2500
  # fdrcut = 0
  tempGRN <- do.call(rbind, lapply(names(orthology_probs), function(thisTF){
    # thisTF <- names(TF_orthology_probs)[10]
    x <- orthology_probs[[thisTF]]
    
    fdrcutx <- names(x[!x > fdrcut])
    
    elegans_x_GR <- MODern_agg_top2500_list[[thisTF]]
    
    # promoters excluding HOT
    promoters_thisone_overlaps <- findOverlaps(elegans_x_GR, promoters)
    
    elegans_x_targets <- promoters[subjectHits(promoters_thisone_overlaps)]
    elegans_x <- promoters_BM[match(elegans_x_targets$gene_id, promoters_BM$entrezgene_id), "wormbase_gseq"]
    elegans_x <- elegans_x[!is.na(elegans_x)]
    
    thiscutoff <- min(cutoff, length(elegans_x))
    
    fdrcut_names <- elegans_x[1:thiscutoff][elegans_x[1:thiscutoff] %in% fdrcutx]
    
    data.frame("source" = rep(thisTF, times = length(fdrcut_names)),
               "target" = fdrcut_names,
               "weight" = rep(1, times = length(fdrcut_names)))
    
  }))
  
  # get rid of duplications
  tempGRN <- tempGRN[!duplicated(paste(tempGRN$source, tempGRN$target)), ]
  
  # let's exclude fewer than 10
  tempGRN <- tempGRN[!tempGRN$source %in% names(table(tempGRN$source))[table(tempGRN$source) < 10], ]
  
  write.table(tempGRN,
              file = paste0("output/GRNs/", prefix, cutoff, "_fdr", fdrcut, ".txt"),
              sep = "\t",
              row.names = FALSE,
              col.names = TRUE)
  
}

#### LOAD PACKAGES & FUNCTIONS ####

## First specify the packages of interest

packages <- c("stringr",
              "biomaRt",
              "TxDb.Celegans.UCSC.ce11.refGene"
              )

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

# set biomaRt mart
parasite_mart <- useMart("parasite_mart", dataset = "wbps_gene", host = "https://parasite.wormbase.org", port = 443)

FIMO_allspp_modERN_denovomotifs_arrays <- readRDS("output/FIMO_allspp_modERN_denovomotifs_arrays.rds")

potential_targets <- row.names(FIMO_allspp_modERN_denovomotifs_arrays[[2]][, , 1])

convert_names_to_entrez <- readRDS("output/convert_names_to_entrez2.rds")
potential_targets_entrez <- convert_names_to_entrez[match(potential_targets, convert_names_to_entrez$wormbase_gseq), "entrezgene_id"]

# potential_promoters <- promoters(genes(TxDb.Celegans.UCSC.ce11.refGene)[genes(TxDb.Celegans.UCSC.ce11.refGene)$gene_id %in% potential_targets_entrez], 
#                        upstream = 1000, 
#                        downstream = 200)

CisBP_TFinfo_withmotif <- readRDS("output/CisBP_TFinfo_withmotif.rds")

allTFS_manualtoptargets_HOTexcl <- readRDS(file = "output/MODernENCODE_manualtoptargets_operon&HOTexcluded.rds")

allTFs_labelused <- readRDS("output/allTFs_labelused.rds")
# 
# allChIP_agg_GR <- readRDS("output/allChIP_agg_GR.rds")
# allChIP_agg_GR_list <- lapply(allTFs_labelused, function(TF_name){
#   # print(TF_name)
#   allChIP_agg_GR[str_detect(allChIP_agg_GR$V4, TF_name)]
#   
# })
# 
# names(allChIP_agg_GR_list) <- names(allTFs_labelused)


#### BUILD ARRAY OF TOP MOTIF TARGETS ####

FIMO_allspp_modERN_denovomotifs_topmotif_list <- lapply(FIMO_allspp_modERN_denovomotifs_arrays, function(thisTF){
  
  # thisTF <- FIMO_allspp_modERN_denovomotifs_arrays[[2]]
  
  thisTF_topmotif <- thisTF[, 1, ]
  
  temp_frame <- data.frame(matrix(nrow = length(potential_targets), ncol = length(dimnames(thisTF)[[3]])))
  row.names(temp_frame) <- potential_targets
  colnames(temp_frame) <- dimnames(thisTF)[[3]]
  
  temp_frame[] <- lapply(colnames(temp_frame), function(this_sp){
    # this_sp <- colnames(temp_frame)[1]
    thisTF_topmotif[, this_sp][potential_targets]
  
  })
  
  temp_frame
  
})

names(FIMO_allspp_modERN_denovomotifs_topmotif_list) <- names(FIMO_allspp_modERN_denovomotifs_arrays)

FIMO_allspp_modERN_denovomotifs_topmotif_list[which(sapply(FIMO_allspp_modERN_denovomotifs_topmotif_list, ncol) == 0)] <- NULL

# convert to array
FIMO_allspp_modERN_denovomotifs_topmotif_array <- array(unlist(FIMO_allspp_modERN_denovomotifs_topmotif_list), dim = c(nrow(FIMO_allspp_modERN_denovomotifs_topmotif_list[[1]]), ncol(FIMO_allspp_modERN_denovomotifs_topmotif_list[[1]]), length(FIMO_allspp_modERN_denovomotifs_topmotif_list)))

dimnames(FIMO_allspp_modERN_denovomotifs_topmotif_array) <- list(row.names(FIMO_allspp_modERN_denovomotifs_topmotif_list[[1]]),
                                                                 colnames(FIMO_allspp_modERN_denovomotifs_topmotif_list[[1]]),
                                                                 names(FIMO_allspp_modERN_denovomotifs_topmotif_list))

# transpose to match elegans-only arrays from previous scripts i.e. third dimension is species

FIMO_allspp_modERN_denovomotifs_topmotif_array <- aperm(FIMO_allspp_modERN_denovomotifs_topmotif_array, c(1, 3, 2))

#### MAKE SECOND ARRAY WITH CISBP MOTIF FOR THOSE THAT HAVE IT #####

### NO in fact do this at the stage of the TF_orthology probabilities. Because it doesnt make sense to compute them again, they will be the same.

# another question is whether to do it with only direct motifs or indirect too...

TFs_with_ChIP_and_directmotif <- names(FIMO_allspp_modERN_denovomotifs_arrays)[names(FIMO_allspp_modERN_denovomotifs_arrays) %in% CisBP_TFinfo_withmotif[CisBP_TFinfo_withmotif$TF_Status == "D", "wormbase_gseq"]]

# 90 of 360 or 25% have a direct motif from CisBP

#### CALCULATE ORTHOLOGY PROBABILITIES ####
totalb4 <- Sys.time()
STREMEtopmotif_orthology_probs <- lapply(dimnames(FIMO_allspp_modERN_denovomotifs_topmotif_array)[[2]], function(thisTFname){
  # thisTFname <- dimnames(FIMO_allspp_modERN_denovomotifs_topmotif_array)[[2]][2]
  b4 <- Sys.time()
  temp_it <- FIMO_allspp_modERN_denovomotifs_topmotif_array[, thisTFname, ]
  
  message(paste0("Doing ", thisTFname, " which is number ", which(dimnames(FIMO_allspp_modERN_denovomotifs_topmotif_array)[[2]] == thisTFname), " of ", length(dimnames(FIMO_allspp_modERN_denovomotifs_topmotif_array)[[2]])))
  
  # note here checked the ranking is as I like it. ie. for ties, the min is true (note because negative vector)
  # i.e. from 100, if number 2 and 3 tie, they both get rank 3 i.e. top 0.03%
  
  percentiles <- apply(temp_it, 2, function(Y){
    
    temprank <- base::rank(-Y[!is.na(Y)], ties.method = "max")
    
    Y[!is.na(Y)] <- temprank / length(Y[!is.na(Y)])
    
    Y
    
  })
  
  target_products <- apply(percentiles, 1, prod, na.rm = TRUE)
  
  # compare this product to a distribution 
  
  shuffprods <- sapply(1:10000, function(i){
    
    shuffled <- apply(temp_it[, 1:10], 2, function(X){
      
      X[!is.na(X)] <- sample(X[!is.na(X)], size = sum(!is.na(X)), replace = FALSE)
      X
      
    })
    
    shuff_percentiles <- apply(shuffled, 2, function(Y){
      
      temprank <- base::rank(-Y[!is.na(Y)], ties.method = "max")
      
      Y[!is.na(Y)] <- temprank / length(Y[!is.na(Y)])
      
      Y
      
    })
    
    apply(shuff_percentiles, 1, prod, na.rm = TRUE)
    
  })
  
  product_vs_shuffle <- sapply(1:length(target_products), function(i){
    
    ecdf_fun(shuffprods[i, ], target_products[i])
    
  })
  
  names(product_vs_shuffle) <- names(target_products)
  
  adjusted_prod_vs_shuff <- p.adjust(product_vs_shuffle, method = "BH")
  names(adjusted_prod_vs_shuff) <- names(product_vs_shuffle)
  
  message("Time for this one:")
  print(Sys.time() - b4)

  if(exists("totalb4")){
  message("Total time up to now:")
  print(Sys.time() - totalb4)
  }
  
  return(adjusted_prod_vs_shuff)
  
})

names(STREMEtopmotif_orthology_probs) <- dimnames(FIMO_allspp_modERN_denovomotifs_topmotif_array)[[2]]
saveRDS(STREMEtopmotif_orthology_probs, "output/STREMEtopmotif_orthology_probs10k.rds")
# STREMEtopmotif_orthology_probs <- readRDS("output/STREMEtopmotif_orthology_probs10k.rds")

#### MAKE GRNs TO JUDGE BY ####

cutoffs_vec <- c(500, 1000, 1500, 2000, 2500, 3000, 5000)
fdrvec <- c(0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8)

for(i in 1:length(cutoffs_vec)){
  
  for(j in 1:length(fdrvec)){
    # i <- 1
    # j <- 3
    message(paste0("Now doing cutoff ", cutoffs_vec[i], " with FDR ", fdrvec[j]))
    make.TF.ChIP.orth.GRN(cutoff = cutoffs_vec[i],
                          fdrcut = fdrvec[j],
                          orthology_probs = STREMEtopmotif_orthology_probs,
                          prefix = "TF_MODernSTREMEorthprobs_cut",
                          target_list = allTFS_manualtoptargets_HOTexcl,
                          min_target_cut = 300)

  }
  
}

# sum(sapply(allTFS_manualtoptargets_HOTexcl, length) > 300)
# hist(sapply(allTFS_manualtoptargets_HOTexcl, length))
