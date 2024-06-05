#start with empty workspace

rm(list = ls(all = TRUE))

# clear all loaded packages
invisible(lapply(paste0("package:", names(sessionInfo()$otherPkgs)),
                 detach,
                 character.only = TRUE, unload = TRUE))

# turn off scientific notation for plots

options(scipen=10000)

#### set working directory ####

# here create new folder and set working directory within it

setwd("~/Cel_GRN_manuscript/")

#### DEFINE FUNCTIONS ####

#### LOAD PACKAGES & FUNCTIONS ####

## First specify the packages of interest

packages <- c("openxlsx",
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

parasite_mart <- useMart("parasite_mart", dataset = "wbps_gene", host = "https://parasite.wormbase.org", port = 443)

walhoutEV1 <- openxlsx::read.xlsx("~/Downloads/msb167131-sup-0002-datasetev1.xlsx")

highqualwalhoutEV1 <- walhoutEV1[walhoutEV1$`In.high-quality.dataset?` == "yes", ]

#### DEFINE TARGETS ####

walhoutTFs <- names(table(highqualwalhoutEV1$Prey.sequence.Name))

walhoutBM <- getBM(values = walhoutTFs,
                   attributes = c("wormbase_gseq",
                                  "wormbase_locus",
                                  "wormbase_gene"),
                   filters = "wormbase_gseqname",
                   mart = parasite_mart)

saveRDS(walhoutBM, "output/walhoutBM.rds")

walhoutGRN <- highqualwalhoutEV1[, c("Prey.sequence.Name", "Bait.sequence.name")]
colnames(walhoutGRN) <- c("source", "target")

walhoutGRN <- walhoutGRN[!duplicated(paste(walhoutGRN$source, walhoutGRN$target)), ]

walhoutGRN[, "weight"] <- 1

write.table(walhoutGRN,
            "output/GRNs/walhoutGRN_unfiltered.txt",
            col.names = TRUE,
            row.names = FALSE,
            sep = "\t")

