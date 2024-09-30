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

#### LOAD PACKAGES & FUNCTIONS ####

## First specify the packages of interest

packages <- c("epitools")

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

TForth_cut1500_fdr0.8_GRN <- read.table("~/Cel_GRN_manuscript/output/GRNs/TF_orthprobs_cut1500_fdr0.8.txt",
                                        header = TRUE,
                                        sep = "\t")

FIMO1000 <- read.table("~/Cel_GRN_manuscript/output/GRNs/FIMO_nohomo_1000.txt",
                       header = TRUE,
                       sep = "\t")

FIMO_nohomo_calcs <- readRDS("~/Cel_GRN_manuscript/output/FIMO_nohomo_calcs.rds")

CelEsT <- read.table("output/GRNs/allthree_equalweights.txt",
                     sep = '\t',
                     header = TRUE)

# make control network from FIMO1000 but with top targets

controlGRN <- do.call(rbind, lapply(unique(TForth_cut1500_fdr0.8_GRN$source)[unique(TForth_cut1500_fdr0.8_GRN$source) %in% names(FIMO_nohomo_calcs)], function(x){

  tempnumber <- sum(TForth_cut1500_fdr0.8_GRN$source == x)
  
  if(tempnumber == 0){return(NULL)}
  
  temptargets <- names(FIMO_nohomo_calcs[[x]][1:tempnumber])
  
  temp_GRN <- data.frame(source = rep(x, times = length(temptargets)),
                         target = temptargets)
  
  temp_GRN[, "weight"] <- 1
  
  temp_GRN
  
}))

write.table(controlGRN,
            "output/GRNs/TF_orthprobs_cut1500_fdr0.8_controlnet.txt",
            row.names = FALSE,
            col.names = TRUE,
            sep = "\t")

# CelEsT control

CelEsT_control <- CelEsT[CelEsT$source %in% unique(TForth_cut1500_fdr0.8_GRN$source), ]

write.table(CelEsT_control,
            "output/GRNs/TF_orthprobs_cut1500_fdr0.8_CelEsTcontrolnet.txt",
            row.names = FALSE,
            col.names = TRUE,
            sep = "\t")

# limit full FIMO to same number of TFs

FIMO1000_limit <- FIMO1000[FIMO1000$source %in% unique(TForth_cut1500_fdr0.8_GRN$source), ]

write.table(FIMO1000_limit,
            "output/GRNs/TF_orthprobs_cut1500_fdr0.8_fullFIMOcontrolnet.txt",
            row.names = FALSE,
            col.names = TRUE,
            sep = "\t")

#### more likely to be in the ChIP?

ChIPGRN <- read.table("~/Cel_GRN_manuscript/output/GRNs/allChIP_10000_HOTexcl_unfiltered.txt",
                      header = TRUE,
                      sep = "\t")

sharedwithChIP <-sapply(unique(TForth_cut1500_fdr0.8_GRN$source), function(x){
  x %in% ChIPGRN$source
})

sharedwithChIP <- names(sharedwithChIP)[sharedwithChIP]

inChIP <-t(sapply(sharedwithChIP, function(x){

  thisnet_targ <- TForth_cut1500_fdr0.8_GRN[TForth_cut1500_fdr0.8_GRN$source == x, "target"]
  
  control_targ <- controlGRN[controlGRN$source == x, "target"]
  
  ChIP_targ <- ChIPGRN[ChIPGRN$source == x, "target"]
  
  FullFIMO_targ <- FIMO1000_limit[FIMO1000_limit$source == x, "target"]
  
  
  c("TForth_in_ChIP" = sum(thisnet_targ %in% ChIP_targ),
             "control_in_ChIP" = sum(control_targ %in% ChIP_targ),
    "FullFIMO_in_ChIP" = sum(FullFIMO_targ %in% ChIP_targ))
  
  
}))

notinChIP <- t(sapply(sharedwithChIP, function(x){
  
  thisnet_targ <- TForth_cut1500_fdr0.8_GRN[TForth_cut1500_fdr0.8_GRN$source == x, "target"]
  
  control_targ <- controlGRN[controlGRN$source == x, "target"]
  
  ChIP_targ <- ChIPGRN[ChIPGRN$source == x, "target"]
  
  FullFIMO_targ <- FIMO1000_limit[FIMO1000_limit$source == x, "target"]
  
  c("TForth_notin_ChIP" = sum(!thisnet_targ %in% ChIP_targ),
             "control_notin_ChIP" = sum(!control_targ %in% ChIP_targ),
    "FullFIMO_notin_ChIP" = sum(!FullFIMO_targ %in% ChIP_targ))
  
  
}))

comparewithallFIMO <- rbind(colSums(inChIP)[c(1,3)],
                            colSums(notinChIP)[c(1,3)])

test <- chisq.test(comparewithallFIMO)

options(scipen = 1)

test$p.value

oddsratio.wald(comparewithallFIMO)
# significant and odds ratio about 1.7
