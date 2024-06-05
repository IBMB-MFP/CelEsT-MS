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

packages <- c("openxlsx",
              "stringr",
              "biomaRt")

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

fullset_TFs_BM <- readRDS("output/fullset_TFs_BM.rds")

walhoutEV1 <- openxlsx::read.xlsx("~/Downloads/msb167131-sup-0002-datasetev1.xlsx")

# limit to interactions described as high quality
walhoutEV1 <- walhoutEV1[walhoutEV1$`In.high-quality.dataset?` == 'yes', ]

parasite_mart <- useMart("parasite_mart", dataset = "wbps_gene", host = "https://parasite.wormbase.org", port = 443)

wTF3 <- read.xlsx("input/WTF3.xlsx")

#### PROCESS INTO GRN ####

walhoutTFs <- unique(walhoutEV1$Prey.sequence.Name)

sum(walhoutTFs %in% fullset_TFs_BM$wormbase_gseq)
sum(!walhoutTFs %in% fullset_TFs_BM$wormbase_gseq)

walhoutTFs_cutoff15  <- names(table(walhoutEV1$Prey.sequence.Name)[table(walhoutEV1$Prey.sequence.Name) >= 15])

# do we have any that are not in the wTF3.0 resource? NO
any(!walhoutTFs_cutoff15 %in% wTF3$Sequence.name)
                              
walhout_BM <- fullset_TFs_BM[match(walhoutTFs_cutoff15, fullset_TFs_BM$wormbase_gseq), ]

walhoutGRN <- walhoutEV1[walhoutEV1$Prey.sequence.Name %in% walhoutTFs_cutoff15, c("Prey.sequence.Name", "Bait.sequence.name")]
colnames(walhoutGRN) <- c("source", "target")

walhoutGRN <- walhoutGRN[!duplicated(paste(walhoutGRN$source, walhoutGRN$target)), ]
walhoutGRN[, "weight"] <- 1

write.table(walhoutGRN,
            "output/GRNs/walhout_highqual_cutoff15.txt",
            col.names = TRUE,
            row.names = FALSE,
            sep = "\t")





