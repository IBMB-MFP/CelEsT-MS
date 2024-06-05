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

packages <- c("eulerr",
              "grid")

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

walhout_unfiltered <- read.table("output/walhoutGRN_unfiltered.txt",
                                 header = TRUE)

walhout_unfiltered_Tfs <- unique(walhout_unfiltered$source)

allmodERN_TFsONLY_BM <- readRDS("output/allmodERN_TFsONLY_BM.rds")
ChIP_TFs <- allmodERN_TFsONLY_BM$wormbase_gseq

separate_by_motif <- readRDS("output/separate_by_motif.rds")
FIMO_Tfs <- names(separate_by_motif)

# how many?
unique(c(walhout_unfiltered_Tfs, ChIP_TFs, FIMO_Tfs))

combo_GRN <- read.table("output/GRNs/allthree_equalweights.txt",
                        header = TRUE)

sum(ChIP_TFs %in% unique(combo_GRN$source))
sum(FIMO_Tfs %in% unique(combo_GRN$source))
sum(walhout_unfiltered_Tfs %in% unique(combo_GRN$source))

# this is incorrect. should be as below
Tfs_alldata_fit <- euler(c("Motif" = sum(!FIMO_Tfs %in% c(walhout_unfiltered_Tfs, ChIP_TFs)),
                "ChIP" = sum(!ChIP_TFs %in% c(walhout_unfiltered_Tfs, FIMO_Tfs)),
                "eY1H" = sum(!walhout_unfiltered_Tfs %in% c(FIMO_Tfs, ChIP_TFs)),
                "Motif&ChIP" = length(base::intersect(FIMO_Tfs, 
                                                   ChIP_TFs)),
                "Motif&eY1H" = length(base::intersect(FIMO_Tfs,
                                                  walhout_unfiltered_Tfs)),
                "ChIP&eY1H" = length(base::intersect(walhout_unfiltered_Tfs, 
                                                   ChIP_TFs)),
                "Motif&ChIP&eY1H" = length(base::intersect(base::intersect(FIMO_Tfs,
                                           walhout_unfiltered_Tfs), ChIP_TFs))))

saveRDS(Tfs_alldata_fit,
        "plotdata/TFs_eulerr.rds")

pdf("graphics/TFs_eulerr.pdf",
    height = 2,
    width = 2)

plot(Tfs_alldata_fit,
     fill = c("dodgerblue", "orange", "purple"),
     labels = list(cex = 0.8),
     quantities = list(cex = 0.5))

dev.off()

# do eulerr plot for interactions

FIMO_nohomo_1500 <- read.table("output/GRNs/FIMO_nohomo_1500.txt",
                               header = TRUE)

ChIP_HOTexcl_nocut <- read.table("output/GRNs/allChIP_10000_HOTexcl_unfiltered.txt",
                                 header = TRUE)

eY1H_net <- read.table("output/GRNs/walhoutGRN_unfiltered.txt",
                       header = TRUE)

interactions_present <- sapply(list(FIMO_nohomo_1500,
                                    ChIP_HOTexcl_nocut,
            eY1H_net), function(x){     
              
              paste0(combo_GRN$source, combo_GRN$target) %in% paste0(x$source, x$target)
              
            })

count_interactions <- apply(interactions_present, 1, function(x){

  c("FIMO" <- x[1] == TRUE & x[2] == FALSE & x[3] == FALSE,
    "ChIP" <- x[1] == FALSE & x[2] == TRUE & x[3] == FALSE,
    "Wal" <- x[1] == FALSE & x[2] == FALSE & x[3] == TRUE,
    "FIMOandChIP" <- x[1] == TRUE & x[2] == TRUE & x[3] == FALSE,
  "FIMOandWal" <- x[1] == TRUE & x[2] == FALSE & x[3] == TRUE,
  "ChIPandWal" <- x[1] == FALSE & x[2] == TRUE & x[3] == TRUE,
  "FIMOandChIPandWal" <- !any(x == FALSE)
  )
  
})

allthree_edges_euler <- euler(c("Motif" = rowSums(count_interactions)[1],
                "ChIP" = rowSums(count_interactions)[2],
                "eY1H" = rowSums(count_interactions)[3],
                "Motif&ChIP" = rowSums(count_interactions)[4],
                "Motif&eY1H" = rowSums(count_interactions)[5],
                "ChIP&eY1H" = rowSums(count_interactions)[6],
                "Motif&ChIP&eY1H" = rowSums(count_interactions)[7]))

saveRDS(allthree_edges_euler,
        "plotdata/allthree_edges_eulerr.rds")

pdf("graphics/edges_eulerr.pdf",
    height = 2,
    width = 2)

plot(allthree_edges_euler,
     fill = c("dodgerblue", "orange", "purple"),
     quantities = list(cex = 0.5))

dev.off()

# make euler plot for final TF contributions in the finished network

wal_in_allthree <- walhout_unfiltered_Tfs[walhout_unfiltered_Tfs %in% combo_GRN$source]
chip_in_allthree <- ChIP_TFs[ChIP_TFs %in% combo_GRN$source]
motif_in_allthree <- FIMO_Tfs[FIMO_Tfs %in% combo_GRN$source]

Tfs_finaldata_fit <- euler(c("Motif" = sum(!motif_in_allthree %in% c(wal_in_allthree, chip_in_allthree)),
                "ChIP" = sum(!chip_in_allthree %in% c(wal_in_allthree, motif_in_allthree)),
                "eY1H" = sum(!wal_in_allthree %in% c(chip_in_allthree, motif_in_allthree)),
                "Motif&ChIP" = length(base::intersect(motif_in_allthree, 
                                                      chip_in_allthree)[!base::intersect(motif_in_allthree, 
                                                                                        chip_in_allthree) %in% base::intersect(base::intersect(motif_in_allthree,
                                                                                                                                               wal_in_allthree), chip_in_allthree)]),
                "Motif&eY1H" = length(base::intersect(motif_in_allthree,
                                                     wal_in_allthree)[!base::intersect(motif_in_allthree,
                                                                                      wal_in_allthree) %in% base::intersect(base::intersect(motif_in_allthree,
                                                                                                                                            wal_in_allthree), chip_in_allthree)]),
                "ChIP&eY1H" = length(base::intersect(wal_in_allthree, 
                                                    chip_in_allthree)[!base::intersect(wal_in_allthree, 
                                                                                       chip_in_allthree) %in% base::intersect(base::intersect(motif_in_allthree,
                                                                                                                                              wal_in_allthree), chip_in_allthree)]),
                "Motif&ChIP&eY1H" = length(base::intersect(base::intersect(motif_in_allthree,
                                                                          wal_in_allthree), chip_in_allthree))))


pdf("graphics/TFs_finalnetwork_eulerr.pdf",
    height = 2.222,
    width = 2.222)

p <- plot(Tfs_finaldata_fit,
          fill = c("dodgerblue", "orange", "purple"),
          labels = list(cex = 0.7),
          quantities = list(cex = 0.5))

p$vp$width <- unit(0.9, "npc")
p$vp$height <- unit(0.9, "npc")

p

dev.off()

# make euler plot for interaction sharing among common TFs

commonTFs <- base::intersect(base::intersect(walhout_unfiltered_Tfs, ChIP_TFs), FIMO_Tfs)

combo_GRN_common <- combo_GRN[combo_GRN$source %in% commonTFs, ]
nrow(combo_GRN_common)
unique(combo_GRN$source)

common_interactions_present <- sapply(list(FIMO_nohomo_1500,
                                    ChIP_HOTexcl_nocut,
                                    eY1H_net), function(x){     
                                      
                                      paste0(combo_GRN_common$source, combo_GRN_common$target) %in% paste0(x$source, x$target)
                                      
                                    })

common_count_interactions <- apply(common_interactions_present, 1, function(x){
  
  c("FIMO" <- x[1] == TRUE & x[2] == FALSE & x[3] == FALSE,
    "ChIP" <- x[1] == FALSE & x[2] == TRUE & x[3] == FALSE,
    "Wal" <- x[1] == FALSE & x[2] == FALSE & x[3] == TRUE,
    "FIMOandChIP" <- x[1] == TRUE & x[2] == TRUE & x[3] == FALSE,
    "FIMOandWal" <- x[1] == TRUE & x[2] == FALSE & x[3] == TRUE,
    "ChIPandWal" <- x[1] == FALSE & x[2] == TRUE & x[3] == TRUE,
    "FIMOandChIPandWal" <- !any(x == FALSE)
  )
  
})

allthreecommon_edges_euler <- euler(c("Motif" = rowSums(common_count_interactions)[1],
                                "ChIP" = rowSums(common_count_interactions)[2],
                                "eY1H" = rowSums(common_count_interactions)[3],
                                "Motif&ChIP" = rowSums(common_count_interactions)[4],
                                "Motif&eY1H" = rowSums(common_count_interactions)[5],
                                "ChIP&eY1H" = rowSums(common_count_interactions)[6],
                                "Motif&ChIP&eY1H" = rowSums(common_count_interactions)[7]))

saveRDS(allthreecommon_edges_euler,
        "plotdata/allthreecommon_edges_eulerr.rds")

pdf("graphics/common_edges_eulerr.pdf",
    height = 2,
    width = 2)

plot(allthreecommon_edges_euler,
     fill = c("dodgerblue", "orange", "purple"),
     quantities = list(cex = 0.5))

dev.off()
