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

packages <- c("CHNOSZ",
              "eulerr")

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

CisBP_TFinfo_withmotif <- readRDS("output/CisBP_TFinfo_withmotif.rds")

observations <- read.table("output/benchmark_observations.txt",
                           sep = "\t",
                           header = TRUE)

TF_orth1500fdr5 <- read.table("output/GRNs/TF_orthprobs_cut1500_fdr0.5unfiltered.txt",
                               sep = "\t",
                               header = TRUE)

FIMO1500 <- read.table("output/GRNs/FIMO_nohomo_1500_unfiltered.txt",
                       sep = "\t",
                       header = TRUE)

# load unfiltered GRNs and filter after combination

ChIP_HOTexcl_nocut <- read.table("output/GRNs/allChIP_10000_HOTexcl.txt",
                                sep = "\t",
                                header = TRUE)

# ensure to exclude cfi-1
ChIP_HOTexcl_nocut <- ChIP_HOTexcl_nocut[ChIP_HOTexcl_nocut$source != "T23D8.8", ]

ChIPorth500 <- read.table("output/GRNs/ChIPorth_500_HOTincl_unfiltered.txt",
                          header = TRUE)

ChIPorth1000 <- read.table("output/GRNs/ChIPorth_1000_HOTexcl_unfiltered.txt",
                           header = TRUE)

eY1H_net <- read.table("output/GRNs/walhoutGRN_unfiltered.txt",
                       sep = "\t",
                       header = TRUE)

#### COMBINE NETWORKS ####

allTFs_targets <- readRDS("output/MODernENCODE_manualtoptargets_operon&HOTexcluded.rds")

# combine to form simple network

do.all.GRN.combinations(GRNlist = list(TF_orth1500fdr5,
                                       ChIPorth1000,
                                       eY1H_net),
                        prefix = "orthCelEsT",
                        weight_vec = c(0.2, 0.5, 0.3))

# it's good!

# stats? 

orthCelEsT <- read.table("output/GRNs/orthCelEsT_equalweights.txt", header = TRUE)

allthree_equal <- read.table("output/GRNs/allthree_equalweights.txt", header = TRUE)

length(unique(orthCelEsT$source))
length(unique(allthree_equal$source))

mean(table(orthCelEsT$source))
mean(table(allthree_equal$source))

# do eulerr plot for interactions

interactions_present <- sapply(list(TF_orth1500fdr5,
                                    ChIPorth1000,
                                    eY1H_net), function(x){     
                                      
                                      paste0(orthCelEsT$source, orthCelEsT$target) %in% paste0(x$source, x$target)
                                      
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

orthCelEsT_edges_euler <- euler(c("Motif" = rowSums(count_interactions)[1],
                                "ChIP" = rowSums(count_interactions)[2],
                                "eY1H" = rowSums(count_interactions)[3],
                                "Motif&ChIP" = rowSums(count_interactions)[4],
                                "Motif&eY1H" = rowSums(count_interactions)[5],
                                "ChIP&eY1H" = rowSums(count_interactions)[6],
                                "Motif&ChIP&eY1H" = rowSums(count_interactions)[7]))

saveRDS(orthCelEsT_edges_euler,
        "plotdata/orthCelEsT_edges_eulerr.rds")

pdf("graphics/orthCelEsT_edges_eulerr.pdf",
    height = 2,
    width = 2)

plot(orthCelEsT_edges_euler,
     fill = c("dodgerblue", "orange", "purple"),
     quantities = list(cex = 0.5))

dev.off()

# make euler plot for interaction sharing among common TFs

commonTFs <- base::intersect(base::intersect(unique(eY1H_net$source), unique(ChIPorth1000$source)), unique(TF_orth1500fdr5$source))

orthCelEsT_common <- orthCelEsT[orthCelEsT$source %in% commonTFs, ]

common_interactions_present <- sapply(list(TF_orth1500fdr5,
                                           ChIPorth1000,
                                           eY1H_net), function(x){     
                                             
                                             paste0(orthCelEsT_common$source, orthCelEsT_common$target) %in% paste0(x$source, x$target)
                                             
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

orthCelEsTcommon_edges_euler <- euler(c("Motif" = rowSums(common_count_interactions)[1],
                                      "ChIP" = rowSums(common_count_interactions)[2],
                                      "eY1H" = rowSums(common_count_interactions)[3],
                                      "Motif&ChIP" = rowSums(common_count_interactions)[4],
                                      "Motif&eY1H" = rowSums(common_count_interactions)[5],
                                      "ChIP&eY1H" = rowSums(common_count_interactions)[6],
                                      "Motif&ChIP&eY1H" = rowSums(common_count_interactions)[7]))

saveRDS(orthCelEsTcommon_edges_euler,
        "plotdata/orthCelEsTcommon_edges_eulerr.rds")

pdf("graphics/orthCelEsT_common_edges_eulerr.pdf",
    height = 2,
    width = 2)

plot(orthCelEsTcommon_edges_euler,
     fill = c("dodgerblue", "orange", "purple"),
     quantities = list(cex = 0.5))

dev.off()

#### ANNOTATE NETWORKS FOR SUPP TABLES ####

orthCelEsT_equalweights <- read.table("output/GRNs/orthCelEsT_equalweights.txt",
                                    sep = "\t",
                                    header = TRUE)

orthCelEsT_equalweights[, "with_motif"] <- NA
orthCelEsT_equalweights[orthCelEsT_equalweights$source %in% TF_orth1500fdr5$source, "with_motif"] <- FALSE
orthCelEsT_equalweights[paste0(orthCelEsT_equalweights$source, "-", orthCelEsT_equalweights$target) %in% paste0(TF_orth1500fdr5$source, "-", TF_orth1500fdr5$target), "with_motif"] <- TRUE

orthCelEsT_equalweights[, "in_ChIP"] <- NA
orthCelEsT_equalweights[orthCelEsT_equalweights$source %in% ChIPorth1000$source, "in_ChIP"] <- FALSE
orthCelEsT_equalweights[paste0(orthCelEsT_equalweights$source, "-", orthCelEsT_equalweights$target) %in% paste0(ChIPorth1000$source, "-", ChIPorth1000$target), "in_ChIP"] <- TRUE

orthCelEsT_equalweights[, "in_eY1H"] <- NA
orthCelEsT_equalweights[orthCelEsT_equalweights$source %in% eY1H_net$source, "in_eY1H"] <- FALSE
orthCelEsT_equalweights[paste0(orthCelEsT_equalweights$source, "-", orthCelEsT_equalweights$target) %in% paste0(eY1H_net$source, "-", eY1H_net$target), "in_eY1H"] <- TRUE

write.table(orthCelEsT_equalweights,
            "output/orthCelEsT_annotated.txt",
            sep = "\t",
            col.names = TRUE,
            row.names = FALSE)

#### MAX_CelEsT ####


# combine to form simple network

do.all.GRN.combinations(GRNlist = list(FIMO1500,
                                       ChIPorth500,
                                       eY1H_net),
                        prefix = "orthCelEsTMAXcov",
                        weight_vec = c(0.2, 0.5, 0.3))

# it's good!

# stats? 

MAX_CelEsT <- read.table("output/GRNs/orthCelEsTMAXcov_equalweights.txt", header = TRUE)
length(unique(MAX_CelEsT$source))
allthree_equal <- read.table("output/GRNs/allthree_equalweights.txt", header = TRUE)

length(unique(MAX_CelEsT$source))
length(unique(allthree_equal$source))

mean(table(MAX_CelEsT$source))
mean(table(allthree_equal$source))

# do eulerr plot for interactions

MAX_interactions_present <- sapply(list(FIMO1500,
                                    ChIPorth500,
                                    eY1H_net), function(x){     
                                      
                                      paste0(orthCelEsT$source, orthCelEsT$target) %in% paste0(x$source, x$target)
                                      
                                    })

MAX_count_interactions <- apply(MAX_interactions_present, 1, function(x){
  
  c("FIMO" <- x[1] == TRUE & x[2] == FALSE & x[3] == FALSE,
    "ChIP" <- x[1] == FALSE & x[2] == TRUE & x[3] == FALSE,
    "Wal" <- x[1] == FALSE & x[2] == FALSE & x[3] == TRUE,
    "FIMOandChIP" <- x[1] == TRUE & x[2] == TRUE & x[3] == FALSE,
    "FIMOandWal" <- x[1] == TRUE & x[2] == FALSE & x[3] == TRUE,
    "ChIPandWal" <- x[1] == FALSE & x[2] == TRUE & x[3] == TRUE,
    "FIMOandChIPandWal" <- !any(x == FALSE)
  )
  
})

MAX_CelEsT_edges_euler <- euler(c("Motif" = rowSums(MAX_count_interactions)[1],
                                  "ChIP" = rowSums(MAX_count_interactions)[2],
                                  "eY1H" = rowSums(MAX_count_interactions)[3],
                                  "Motif&ChIP" = rowSums(MAX_count_interactions)[4],
                                  "Motif&eY1H" = rowSums(MAX_count_interactions)[5],
                                  "ChIP&eY1H" = rowSums(MAX_count_interactions)[6],
                                  "Motif&ChIP&eY1H" = rowSums(MAX_count_interactions)[7]))

saveRDS(MAX_CelEsT_edges_euler,
        "plotdata/orthCelEsT_edges_eulerr.rds")

pdf("graphics/MAX_CelEsT_edges_eulerr.pdf",
    height = 2,
    width = 2)

plot(MAX_CelEsT_edges_euler,
     fill = c("dodgerblue", "orange", "purple"),
     quantities = list(cex = 0.5))

dev.off()

# make euler plot for interaction sharing among common TFs

MAX_commonTFs <- base::intersect(base::intersect(unique(eY1H_net$source), unique(ChIPorth500$source)), unique(FIMO1500$source))

MAX_CelEsT_common <- MAX_CelEsT[MAX_CelEsT$source %in% MAX_commonTFs, ]

MAX_common_interactions_present <- sapply(list(FIMO1500,
                                           ChIPorth500,
                                           eY1H_net), function(x){     
                                             
                                             paste0(MAX_CelEsT_common$source, MAX_CelEsT_common$target) %in% paste0(x$source, x$target)
                                             
                                           })

MAX_common_count_interactions <- apply(MAX_common_interactions_present, 1, function(x){
  
  c("FIMO" <- x[1] == TRUE & x[2] == FALSE & x[3] == FALSE,
    "ChIP" <- x[1] == FALSE & x[2] == TRUE & x[3] == FALSE,
    "Wal" <- x[1] == FALSE & x[2] == FALSE & x[3] == TRUE,
    "FIMOandChIP" <- x[1] == TRUE & x[2] == TRUE & x[3] == FALSE,
    "FIMOandWal" <- x[1] == TRUE & x[2] == FALSE & x[3] == TRUE,
    "ChIPandWal" <- x[1] == FALSE & x[2] == TRUE & x[3] == TRUE,
    "FIMOandChIPandWal" <- !any(x == FALSE)
  )
  
})

MAX_CelEsTMAX_common_edges_euler <- euler(c("Motif" = rowSums(MAX_common_count_interactions)[1],
                                        "ChIP" = rowSums(MAX_common_count_interactions)[2],
                                        "eY1H" = rowSums(MAX_common_count_interactions)[3],
                                        "Motif&ChIP" = rowSums(MAX_common_count_interactions)[4],
                                        "Motif&eY1H" = rowSums(MAX_common_count_interactions)[5],
                                        "ChIP&eY1H" = rowSums(MAX_common_count_interactions)[6],
                                        "Motif&ChIP&eY1H" = rowSums(MAX_common_count_interactions)[7]))

saveRDS(MAX_CelEsTcommon_edges_euler,
        "plotdata/MAX_CelEsTcommon_edges_eulerr.rds")

pdf("graphics/MAX_CelEsT_common_edges_eulerr.pdf",
    height = 2,
    width = 2)

plot(MAX_CelEsTcommon_edges_euler,
     fill = c("dodgerblue", "orange", "purple"),
     quantities = list(cex = 0.5))

dev.off()

#### ANNOTATE MAXCELEST ####

MAX_CelEsT_equalweights <- read.table("output/GRNs/orthCelEsTMAXcov_equalweights.txt",
                                      sep = "\t",
                                      header = TRUE)

MAX_CelEsT_equalweights[, "with_motif"] <- NA
MAX_CelEsT_equalweights[MAX_CelEsT_equalweights$source %in% FIMO1500$source, "with_motif"] <- FALSE
MAX_CelEsT_equalweights[paste0(MAX_CelEsT_equalweights$source, "-", MAX_CelEsT_equalweights$target) %in% paste0(FIMO1500$source, "-", FIMO1500$target), "with_motif"] <- TRUE

MAX_CelEsT_equalweights[, "in_ChIP"] <- NA
MAX_CelEsT_equalweights[MAX_CelEsT_equalweights$source %in% ChIPorth500$source, "in_ChIP"] <- FALSE
MAX_CelEsT_equalweights[paste0(MAX_CelEsT_equalweights$source, "-", MAX_CelEsT_equalweights$target) %in% paste0(ChIPorth500$source, "-", ChIPorth500$target), "in_ChIP"] <- TRUE

MAX_CelEsT_equalweights[, "in_eY1H"] <- NA
MAX_CelEsT_equalweights[MAX_CelEsT_equalweights$source %in% eY1H_net$source, "in_eY1H"] <- FALSE
MAX_CelEsT_equalweights[paste0(MAX_CelEsT_equalweights$source, "-", MAX_CelEsT_equalweights$target) %in% paste0(eY1H_net$source, "-", eY1H_net$target), "in_eY1H"] <- TRUE

write.table(MAX_CelEsT_equalweights,
            "output/MAX_CelEsT_annotated.txt",
            sep = "\t",
            col.names = TRUE,
            row.names = FALSE)



