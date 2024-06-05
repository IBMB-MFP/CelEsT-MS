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

packages <- c("biomaRt",
              "dplyr",
              "stringr",
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

biocmanager_packages <- c("GenomicRanges", # for intersecting
                          "TxDb.Celegans.UCSC.ce11.refGene", # for promoter sequences
                          "BSgenome.Celegans.UCSC.ce11"
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

github_packages <- c("r-lib/conflicted",
                     "snystrom/memes") 

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

# this is a function from the now-defunct RADami package
write.DNAStringSet <- function(x,
                               format= c('phylip', 'fasta'),
                               padding = 30,
                               filename = 'DNAStringSetOut.phy',
                               fastaPrefix = '>') {
  
  # writes a sequence matrix to phylip or fasta format
  
  x.width <- width(x)[1]
  
  x <- as.character(x)
  
  if(format[1] == 'phylip') {
    
    for(i in 1:length(x)) x[i] <- paste(names(x)[i], paste(rep(" ", (padding - nchar(names(x)[i]))), collapse = ''), x[i], sep = '')
    
    writeLines(c(paste(length(x), x.width), x), filename)
    
  }
  
  if(format[1] == 'fasta') {
    
    out <- as.character(matrix(c(paste(fastaPrefix, names(x), sep = ''), x), nrow = 2, byrow = T))
    
    writeLines(out, filename)
    
  }
  
  return(0)
  
}

split.GR.byTF.tolist <- function(GR_object){
  
  templist <- lapply(allTFs_labelused, function(thisTF){
    
    GR_object[str_detect(GR_object$experiment, thisTF), ]
    
  })
  
  # could use all_in_peaks_files_gseq_TFonly but just to be sure will match to BM, in case e.g. all peaks from some TF fall into HOT regions etc.
  tempnames <- allTFs_labelused
  tempnames[which(!tempnames %in% allmodERN_TFsONLY_BM$wormbase_gseq)] <- allmodERN_TFsONLY_BM[match(tempnames[!tempnames %in% allmodERN_TFsONLY_BM$wormbase_gseq], allmodERN_TFsONLY_BM$wormbase_locus), "wormbase_gseq"]
  names(templist) <- tempnames
  
  templist
  
}

get.sequences.write.byTF <- function(GR_list,
                                     outdir){
  
  dir.create("output/PeakSeqs",
             showWarnings = FALSE)
  dir.create(paste0("output/PeakSeqs/", outdir))
  
  tempnames <- names(GR_list)
  
  tempseq <-  lapply(tempnames, function(x){
    
    get_sequence(GR_list[[x]],
                 BSgenome.Celegans.UCSC.ce11)
    
  })
  
  names(tempseq) <- tempnames
  
  for(i in 1:length(tempseq)){
    
    write.DNAStringSet(tempseq[[i]],
                       format = "fasta",
                       filename = paste0("output/PeakSeqs/", outdir, "/", names(tempseq)[i], "_peakseq.fasta"))
    
  }
  
}

return.manual.targets <- function(tf_name,
                                  peaks_GRanges = chip_peaks_overlap,
                                  cutoff = 500){
  
  thisTF_peaks_GR <- chip_peaks_overlap[str_detect(chip_peaks_overlap$V4, paste0(tf_name, "_"))]
  
  # first the case in which there is only a single replicate
  if(length(unique(thisTF_peaks_GR$V4[str_detect(thisTF_peaks_GR$V4, paste0(tf_name, "_"))])) == 1){
    
    thisTF_peaks_df <- as.data.frame(thisTF_peaks_GR)
    
    thisTF_peaks_collapse <- thisTF_peaks_df %>% group_by(manual_gseq) %>% summarize(sum = sum(V7))
    
    if(nrow(thisTF_peaks_collapse) > cutoff){
      
      top_targets <- unname(unlist(thisTF_peaks_collapse[order(thisTF_peaks_collapse$sum, decreasing = TRUE), "manual_gseq"])[1:cutoff])
      
    } else {
      
      top_targets <- unname(unlist(thisTF_peaks_collapse[order(thisTF_peaks_collapse$sum, decreasing = TRUE), "manual_gseq"]))
      
    }
    
  }
  
  if(length(unique(thisTF_peaks_GR$V4[str_detect(thisTF_peaks_GR$V4, paste0(tf_name, "_"))])) > 1){
    
    thisTF_stages <- unique(thisTF_peaks_GR$V4[str_detect(thisTF_peaks_GR$V4, paste0(tf_name, "_"))])
    
    thisTF_peaks <- lapply(thisTF_stages, function(thisstage){
      
      thisTF_peaks_GR[str_detect(thisTF_peaks_GR$V4, thisstage), ]
      
    })
    
    allstagetargets <- unique(unlist(lapply(thisTF_peaks, function(zz){zz$manual_gseq})))
    
    number_of_appearances_df <- data.frame("appearances" = sapply(allstagetargets, function(thistarget){
      
      sum(sapply(thisTF_peaks, function(thispeakobject){
        
        thistarget %in% thispeakobject$manual_gseq
        
      }))
      
    }))
    
    # here we look at the targets in each stage, and add up the peaks targeting each target for each stage. 
    # then we return the mean of that sum for each stage overall to give an estimate of signal strength for the target.
    number_of_appearances_df[, "sum_signal"] <- sapply(row.names(number_of_appearances_df), function(thistarget){
      
      thisTF_peaks_subset <- thisTF_peaks[sapply(thisTF_peaks, function(wx){thistarget %in% wx$manual_gseq})]
      
      mean(sapply(thisTF_peaks_subset, function(wy){
        
        wy2 <- as.data.frame(wy[wy$manual_gseq == thistarget, ])
        
        sumtibble <-  wy2 %>% group_by(manual_gseq) %>% summarize(sum = sum(V7))
        
        return(sumtibble$sum)
        
      }))
      
    })
    
    # we order targets primarily by number of stages they are bound at, and then by the average strength of summed peaks across stages
    number_of_appearances_df[] <- number_of_appearances_df[with(number_of_appearances_df, order(appearances, sum_signal, decreasing = TRUE)), ]
    
    if(nrow(number_of_appearances_df) > cutoff){
      
      top_targets <- row.names(number_of_appearances_df)[1:cutoff]
      
    } else {
      
      top_targets <- row.names(number_of_appearances_df)
      
    }
    
  }
  
  return(top_targets)
  
}

run.STREME.on.GR <- function(GRlist,
                             minsequence_cutoff =  100,
                             streme_objfun = "de",
                             BSgenome_obj = BSgenome.Celegans.UCSC.ce11,
                             prefix = "mySTREME") {
  
  if(streme_objfun ==  "cd"){
    
    control_arg <- NA
    
  } else {
    
    control_arg <- "shuffle"
    
  }
  
  tryCatch(
    {
      
      templist <- lapply(names(GRlist), function(thisname){
        message(thisname)
        thisGR <- GRlist[[thisname]]

        seqlengths(thisGR) <- seqlengths(BSgenome_obj)
        thisGR <- trim(thisGR)
        
        tempseq <- get_sequence(thisGR,
                                BSgenome_obj)
        
        if(length(tempseq) < minsequence_cutoff){
          
          return(NA)
          
        }
        
        streme_out <- runStreme(
          tempseq,
          control = control_arg,
          outdir = "auto",
          objfun = streme_objfun,
          meme_path = "/opt/local/bin"
        )
        
        saveRDS(streme_out, paste0("output/STREME_out/", prefix, "_", thisname, "_minseqcut", minsequence_cutoff, ".txt"))
        
        return(streme_out)
        
      })
      
      names(templist) <- names(GRlist)
      
      return(templist)
      
    }, error = function(e){
      
      return("Returned error")
      
    }
  )
  
}


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

#### INPUT DATA ####

parasite_mart <- useMart("parasite_mart", dataset = "wbps_gene", host = "https://parasite.wormbase.org", port = 443)

MODern_allTFs_data <- read.xlsx("input/Kudron2024_biorXiv_Supp1.xlsx",
                                sheet = 6)

# MODern ChIP-seq peaks file downloaded from https://epic.gs.washington.edu/modERNresource/, tab = Worm / ChIPSeq by Gene / C. elegans + TFs in Kudron 2024 / Download Optimal Peaks For The Experiments
MODern_agg <- read.table("input/PublishedOptimalWormPeaks",
                         sep = "\t",
                         header = FALSE,
                         fill = TRUE)

# Extra ENCODE TFs data
# There are 6 genes listed in the Kudron et al 2024 file as having good ChIP but they are not included in the Peaks file (presumably by omission)
# 4 are present in the ENCODE platform

# These are:
# B0019.2
# W05B10.2 / ccch-3
# T26A8.4
# W05H7.4 / zfp-3

# additionally two more:

# W03F9.2
# Y116A8C.19

# are present in original MODern release but not in most recent file, although no sign they have been revoked in Kudron 2024 supplementary

ENCODE_notinKudron <- c("B0019.2",
                        "ccch-3",
                        "T26A8.4",
                        "zfp-3")

ENCODE_filelist <- list.files("~/Cel_GRN_manuscript/input/ENCODE_extraTFs/")

# Limit list (produced before Kudron 2024 published) to include only those additional TFs

# This is a file compiled manually from ENCODE
ENCODE_metadata <- read.xlsx("input/encode_extra_TFs.xlsx")

ENCODE_filelist_notinKudron <- ENCODE_filelist[sapply(ENCODE_filelist, function(X){any(str_detect(X, ENCODE_metadata[ENCODE_metadata$Target.of.assay %in% ENCODE_notinKudron, "bednarrowpeakfile"]))})]

ENCODE_datalist <- lapply(ENCODE_filelist_notinKudron, function(x){
  
  read.table(paste0("~/Cel_GRN_manuscript/input/ENCODE_extraTFs/", x),
             sep = '\t',
             header = FALSE,
             fill = TRUE)
  
})

names(ENCODE_datalist) <- ENCODE_filelist_notinKudron

# process to have gene and stage information in the 4th column as for modERN processed aggregated peaks

ENCODE_datalist <- lapply(ENCODE_filelist_notinKudron, function(thisone){
  
  thisdata <- ENCODE_datalist[[thisone]]
  thismetadata <- unlist(ENCODE_metadata[match(str_remove(thisone, "\\.bed\\.gz"), ENCODE_metadata$bednarrowpeakfile), ])
  
  thisdata[, 4] <- paste0(thismetadata["Target.gene.symbol"], "_", str_replace_all(unlist(thismetadata["Life.stage"]), pattern = " ", replacement = "_"))
  
  thisdata
  
})

ENCODE_aggregated <- do.call(rbind, ENCODE_datalist)

## Now to recover the remaining files that had been missing from this release:
# W03F9.2
# Y116A8C.19

MODern_missingfrom2024biorxivrelease <- c("W03F9.2",
                                          "Y116A8C.19"
)

MODern_old <- read.table("input/annotatedPeak.bed",
                         sep = "\t",
                         header = FALSE,
                         fill = TRUE)

MODern_old_missing <- MODern_old[sapply(MODern_old$V4, function(X){any(str_detect(X, MODern_missingfrom2024biorxivrelease))}), ]

## Extra TF ChIPs from literature
## Pubmed IDs:
# 31346165
# 36897776
# 30956009
# 34499028
# 28874466
# 25773600

lit_chips_filelist <- list.files("input/additional_literature_chips/")

## REPLACE with supplementary for manuscript
lit_chips_metadata <- read.table("input/additional_literature_chips_filesforpeakcaller.txt",
                                 sep = "\t",
                                 header = TRUE)

lit_chips_datalist <- lapply(lit_chips_filelist, function(thisone){
  
  thisdata <- read.table(paste0("input/additional_literature_chips/", thisone),
                         sep = '\t',
                         header = FALSE,
                         fill = TRUE)
  
  thisname <- lit_chips_metadata[match(str_remove(thisone, "\\.narrowPeak"), lit_chips_metadata$Accession), "target.TF"]
  
  thisdata[, 4] <- paste0(thisname, "_unknownstage")
  
  thisdata
  
})

lit_chips_aggregated <- do.call(rbind, lit_chips_datalist)

# must rename chromosomes to match the ENCODE files

lit_chips_aggregated[lit_chips_aggregated$V1 != "MtDNA", 1] <- paste0("chr", lit_chips_aggregated[lit_chips_aggregated$V1 != "MtDNA", 1])
lit_chips_aggregated[lit_chips_aggregated$V1 == "MtDNA", 1] <- "chrM"

# combine all

allChIP_agg <- rbind(MODern_agg[, 1:10],
                     ENCODE_aggregated,
                     MODern_old_missing[, 1:10],
                     lit_chips_aggregated
)

#### MAKE GRANGES FROM PEAKS FILE ####

allTFs_labelused <- readRDS("output/allTFs_labelused.rds")

allmodERN_TFsONLY_BM <- readRDS("output/allmodERN_TFsONLY_BM.rds")

allChIP_agg_for_GR <- allChIP_agg
colnames(allChIP_agg)[1:3] <-  c("chromosome", "start", "end")

allChIP_agg_GR <- makeGRangesFromDataFrame(allChIP_agg, 
                                           keep.extra.columns = TRUE)

allChIP_agg_peaksummits_GR <- makeGRangesFromDataFrame(data.frame("chrom" = allChIP_agg$chromosome,
                                     "start" = (allChIP_agg$start + allChIP_agg$V10) - 50,
                                     "end" = (allChIP_agg$start + allChIP_agg$V10) + 50,
                                     "experiment" = allChIP_agg$V4,
                                     "score" = allChIP_agg$V7),
                                     keep.extra.columns = TRUE)

# # also do wide version for STREME's centrality mode
# allChIP_agg_peaksummits_wideGR <- makeGRangesFromDataFrame(data.frame("chrom" = allChIP_agg$V1,
#                                                                  "start" = (allChIP_agg$start + allChIP_agg$V10) - 250,
#                                                                  "end" = (allChIP_agg$start + allChIP_agg$V10) + 250,
#                                                                  "experiment" = allChIP_agg$V4),
#                                                       keep.extra.columns = TRUE)

### MODERN_agg is ENCODE narrowPeak format. peak summit is in column 10. 0-based offset from start

## Will want to exclude HOT regions. This is Table S10 of Kudron et al.

# download.file("https://figshare.com/ndownloader/files/10099716",
#               destfile = "~/Cel_GRN_orthology/input/Kudron_et_al_2018_modERN_TableS10.txt")

# The file has different number of columns because it lists on each line all the experiments which show a peak in the HOT region in question
#  just want the first 4 lines; chrom, start, end, number of peaks, without worrying too much which ones they are
# solution code Kudos to https://stackoverflow.com/questions/68139359/how-to-read-first-four-columns-from-a-file-with-different-number-of-columns-on-e

# read the file as text lines
HOTtxt_lines <- readLines("~/Cel_GRN_orthology/input/Kudron_et_al_2018_modERN_TableS10.txt")
# split by one or more spaces
HOTtxt <- base::strsplit(HOTtxt_lines, " +")
# keep only the vector elements with more than 0 chars
HOTtxt <- lapply(HOTtxt, function(x) x[sapply(x, nchar) > 0])
# the last line may have a '\n' only, remove it
HOTtxt <- HOTtxt[lengths(HOTtxt) > 0]
# now extract the first 4 elements of each vector
HOTtxt <- lapply(HOTtxt, '[', 1:4)
# and rbind to data.frame
HOT_df <- do.call(rbind.data.frame, HOTtxt)
names(HOT_df) <- c("chrom", "start", "end", "peakno")

HOT_df$peakno <- as.numeric(HOT_df$peakno)

## I think I don't like their conception of it - because there are some with way more peaks than TFs assayed
# for which I assume that the nuumber of peaks is the nuber of experiments, not the number of different TFs
# and I would rather the different TFs as replicate / diff stage for same TF binding in same place should not be any kind of worry

HOT_exp <- lapply(HOTtxt_lines, function(x){

  x_split <- base::strsplit(x, " +")
  
  if(unlist(x_split)[4] == 1){
    
    return(NULL)
    
  }
  
  x_temp <- unlist(x_split)[5:length(unlist(x_split))]
  
  ## WRAP in str_remove for two freaky exceptions that otherwise dont make much sense; remove the fem--2 (presumably refers to him-8 background) and remove xtl1186_
  # x_TFs <- str_remove(str_remove(str_remove(x_temp, "ce_"), "[A-Z]{2,3}[0-9]+_"), "_[A-Za-z0-9]+_[A-Z]{2}:.*$")
  x_TFs <- str_remove(str_remove(str_remove(str_remove(str_remove(x_temp, "ce_"), "[A-Z]{2,3}[0-9]+_"), "fem-2_"), "xtl1186_"), "_[A-Za-z0-9]+_[A-Z]{2}:.*$")

  unique(x_TFs)
  
})

HOT_df[, "distinctpeak_no"] <- sapply(HOT_exp, length)

# needs mitochondrial DNA to be same as in the other one
HOT_df[HOT_df$chrom == "chrMtDNA", "chrom"] <- "chrM"

# for my purposes want to define HOT peaks as.... 
# original modENCODE Gerstien et al. 2010 Science 330:6012 defined HOT as >15 TFs. If I take this definition, it refers to 7110 peaks/regions, leaving 51706 not HOT. That seems ok.
sum(HOT_df$distinctpeak_no > 50)
# sum(HOT_df$distinctpeak_no < 15)

HOT_cluster_GR <- makeGRangesFromDataFrame(HOT_df[HOT_df$distinctpeak_no > 50, 1:3])
allChIP_agg_peaksummits_HOToverlaps <- findOverlaps(allChIP_agg_peaksummits_GR, HOT_cluster_GR)
allChIP_agg_peaksummits_notHOT_GR <- allChIP_agg_peaksummits_GR[-queryHits(allChIP_agg_peaksummits_HOToverlaps), ]
# allChIP_agg_peaksummits_notHOT_wideGR <- allChIP_agg_peaksummits_wideGR[-queryHits(allChIP_agg_peaksummits_HOToverlaps), ]

# # get promoter sequences 
# 
# promoters <- promoters(genes(TxDb.Celegans.UCSC.ce11.refGene), 
#                        upstream = 1000, 
#                        downstream = 200)
# 
# # promoters excluding HOT
# promoters_MODnotHOT_overlaps <- findOverlaps(allChIP_agg_peaksummits_notHOT_GR, promoters)
# allChIP_agg_peaksummits_notHOTpromoters_GR <- allChIP_agg_peaksummits_notHOT_GR[queryHits(promoters_MODnotHOT_overlaps), ]
# # allChIP_agg_peaksummits_notHOTpromoters_wideGR <- allChIP_agg_peaksummits_notHOT_wideGR[queryHits(promoters_MODnotHOT_overlaps), ]

# # promoters including HOT
# promoters_MODall_overlaps <- findOverlaps(allChIP_agg_peaksummits_GR, promoters)
# allChIP_agg_peaksummits_promoters_GR <- allChIP_agg_peaksummits_GR[queryHits(promoters_MODall_overlaps), ]
# allChIP_agg_peaksummits_promoters_wideGR <- allChIP_agg_peaksummits_wideGR[queryHits(promoters_MODall_overlaps), ]

# so I have four modalities - all peaks, all peaks not HOT, all peaks in promoters, all peaks in promoters not HOT
# initially can do this only for the ones which overlap CisBP as this is my test case, to see how well we perform in terms of motif identification. 
# This will then determine how we proceed with the lot (in terms of the four modalities)

# MODern_agg_peaksummits_GR
# MODern_agg_peaksummits_promoters_GR
# MODern_agg_peaksummits_notHOT_GR
# MODern_agg_peaksummits_notHOTpromoters_GR

allpeaks_by_TF_GRlist <- split.GR.byTF.tolist(allChIP_agg_peaksummits_GR)

notHOTpeaks_by_TF_GRlist <- split.GR.byTF.tolist(allChIP_agg_peaksummits_notHOT_GR)

# allpromoterpeaks_by_TF_GRlist <- split.GR.byTF.tolist(allChIP_agg_peaksummits_promoters_GR)
# notHOTpromoterpeaks_by_TF_GRlist <- split.GR.byTF.tolist(MODern_agg_peaksummits_notHOTpromoters_GR)
# 
# allpeaks_by_TF_wideGRlist <- split.GR.byTF.tolist(MODern_agg_peaksummits_wideGR)
# notHOTpeaks_by_TF_wideGRlist <- split.GR.byTF.tolist(MODern_agg_peaksummits_notHOT_wideGR)
# allpromoterpeaks_by_TF_wideGRlist <- split.GR.byTF.tolist(MODern_agg_peaksummits_promoters_wideGR)
# notHOTpromoterpeaks_by_TF_wideGRlist <- split.GR.byTF.tolist(MODern_agg_peaksummits_notHOTpromoters_wideGR)
# 
# # limit to TFs with direct motif in CisBP for now

CisBP_TFinfo_withmotif <- readRDS("output/CisBP_TFinfo_withmotif.rds")

# # here write sequences to FASTA for STREME's command line tool
# get.sequences.write.byTF(GR_list = allpeaks_by_TF_GRlist[names(allpeaks_by_TF_GRlist) %in% CisBP_TFinfo_withmotif[CisBP_TFinfo_withmotif$TF_Status == "D", "wormbase_gseq"]],
#                          outdir = "allpeaks")
# 
# get.sequences.write.byTF(GR_list = notHOTpeaks_by_TF_GRlist[names(notHOTpeaks_by_TF_GRlist) %in% CisBP_TFinfo_withmotif[CisBP_TFinfo_withmotif$TF_Status == "D" & CisBP_TFinfo_withmotif$Motif_Type == "PBM", "wormbase_gseq"]],
#                          outdir = "notHOTpeaks")
# 
# get.sequences.write.byTF(GR_list = allpromoterpeaks_by_TF_GRlist[names(allpromoterpeaks_by_TF_GRlist) %in% CisBP_TFinfo_withmotif[CisBP_TFinfo_withmotif$TF_Status == "D", "wormbase_gseq"]],
#                          outdir = "allpromoterpeaks")
# 
# get.sequences.write.byTF(GR_list = notHOTpromoterpeaks_by_TF_GRlist[names(notHOTpromoterpeaks_by_TF_GRlist) %in% CisBP_TFinfo_withmotif[CisBP_TFinfo_withmotif$TF_Status == "D", "wormbase_gseq"]],
#                          outdir = "notHOTpromoterpeaks")

#### TRY RUNNNING STREME IN R ####

# 
# STREME_notHOTpromoterpeaks <- run.STREME.on.GR(GRlist = notHOTpromoterpeaks_by_TF_GRlist[names(notHOTpromoterpeaks_by_TF_GRlist) %in% CisBP_TFinfo_withmotif[CisBP_TFinfo_withmotif$TF_Status == "D", "wormbase_gseq"]],
#                  streme_objfun = "de",
#                  minsequence_cutoff =  100,
#                  prefix = "STREME_notHOTpromoterpeaks")
# 
# STREME_notHOTpeaks <- run.STREME.on.GR(GRlist = notHOTpeaks_by_TF_GRlist[names(notHOTpeaks_by_TF_GRlist) %in% CisBP_TFinfo_withmotif[CisBP_TFinfo_withmotif$TF_Status == "D", "wormbase_gseq"]],
#                                                streme_objfun = "de",
#                                                minsequence_cutoff =  100,
#                                        prefix = "STREME_notHOTpeaks")
# # 
# # STREME_allpromoterpeaks <- run.STREME.on.GR(GRlist = allpromoterpeaks_by_TF_GRlist[names(allpromoterpeaks_by_TF_GRlist) %in% CisBP_TFinfo_withmotif[CisBP_TFinfo_withmotif$TF_Status == "D", "wormbase_gseq"]],
# #                                                streme_objfun = "de",
# #                                                minsequence_cutoff =  100,
# #                                             prefix = "STREME_allpromoterpeaks")
# # 
# STREME_allpeaks <- run.STREME.on.GR(GRlist = allpeaks_by_TF_GRlist[names(allpeaks_by_TF_GRlist) %in% CisBP_TFinfo_withmotif[CisBP_TFinfo_withmotif$TF_Status == "D", "wormbase_gseq"]],
#                                                streme_objfun = "de",
#                                                minsequence_cutoff =  100,
#                                     prefix = "STREME_allpeaks")
# 
# saveRDS(STREME_notHOTpromoterpeaks, "output/STREME_notHOTpromoterpeaks.rds")
# saveRDS(STREME_notHOTpeaks, "output/STREME_notHOTpeaks.rds")
# saveRDS(STREME_allpromoterpeaks , "output/STREME_allpromoterpeaks.rds")
# saveRDS(STREME_allpeaks, "output/STREME_allpeaks.rds")
# 
# STREME_notHOTpromoterpeaks_cd <- run.STREME.on.GR(GRlist = notHOTpromoterpeaks_by_TF_wideGRlist[names(notHOTpromoterpeaks_by_TF_wideGRlist) %in% CisBP_TFinfo_withmotif[CisBP_TFinfo_withmotif$TF_Status == "D", "wormbase_gseq"]],
#                                                streme_objfun = "cd",
#                                                minsequence_cutoff =  100)
# 
# # run wide GRs with central distance to check performance
# STREME_notHOTpeaks_cd <- run.STREME.on.GR(GRlist = notHOTpeaks_by_TF_wideGRlist[names(notHOTpeaks_by_TF_wideGRlist) %in% CisBP_TFinfo_withmotif[CisBP_TFinfo_withmotif$TF_Status == "D", "wormbase_gseq"]],
#                                        streme_objfun = "cd",
#                                        minsequence_cutoff =  100)
# 
# STREME_allpromoterpeaks_cd <- run.STREME.on.GR(GRlist = allpromoterpeaks_by_TF_wideGRlist[names(allpromoterpeaks_by_TF_wideGRlist) %in% CisBP_TFinfo_withmotif[CisBP_TFinfo_withmotif$TF_Status == "D", "wormbase_gseq"]],
#                                             streme_objfun = "cd",
#                                             minsequence_cutoff =  100)
# 
# STREME_allpeaks_cd <- run.STREME.on.GR(GRlist = allpeaks_by_TF_wideGRlist[names(allpeaks_by_TF_wideGRlist) %in% CisBP_TFinfo_withmotif[CisBP_TFinfo_withmotif$TF_Status == "D", "wormbase_gseq"]],
#                                     streme_objfun = "cd",
#                                     minsequence_cutoff =  100)      
# 
# saveRDS(STREME_notHOTpromoterpeaks_cd, "output/STREME_notHOTpromoterpeaks_cd.rds")
# saveRDS(STREME_notHOTpeaks_cd, "output/STREME_notHOTpeaks_cd.rds")
# saveRDS(STREME_allpromoterpeaks_cd , "output/STREME_allpromoterpeaks_cd.rds")
# saveRDS(STREME_allpeaks_cd, "output/STREME_allpeaks_cd.rds")

#### LIMIT SEQUENCES AND FIND GOOD WAY TO DO IT ####

notHOT_top_regions <- lapply(names(notHOTpeaks_by_TF_GRlist), function(X){
  print(X)
  return.top.chip.peaks(X,
                        cutoff = 300,
                        GR_list = notHOTpeaks_by_TF_GRlist)
})

names(notHOT_top_regions) <- names(notHOTpeaks_by_TF_GRlist)
saveRDS(notHOT_top_regions, "output/notHOT_top_regions.rds")
# notHOT_top_regions <- readRDS("output/notHOT_top_regions.rds")

#### OK then I will want to run STREME on these regions

STREME_notHOT50_topregions <- run.STREME.on.GR(GRlist = notHOT_top_regions,
                                       streme_objfun = "de",
                                       minsequence_cutoff =  100)

saveRDS(STREME_notHOT50_topregions, "output/STREME_notHOT50_topregions_ALLMODern.rds") 

notHOT_top_regions_nocutoff <- lapply(names(notHOTpeaks_by_TF_GRlist), function(X){
  print(X)
  return.top.chip.peaks(X,
                        cutoff = 10000,
                        GR_list = notHOTpeaks_by_TF_GRlist)
})

names(notHOT_top_regions_nocutoff) <- names(notHOTpeaks_by_TF_GRlist)
saveRDS(notHOT_top_regions_nocutoff, "output/notHOT_top_regions_nocutoff.rds")
# notHOT_top_regions_nocutoff <- readRDS("output/notHOT_top_regions_nocutoff.rds")

STREME_notHOT50nocutoff_topregions <- run.STREME.on.GR(GRlist = notHOT_top_regions_nocutoff,
                                               streme_objfun = "de",
                                               minsequence_cutoff =  100)

saveRDS(STREME_notHOT50nocutoff_topregions, "output/STREME_notHOT50nocutoff_topregions_ALLMODern.rds")
