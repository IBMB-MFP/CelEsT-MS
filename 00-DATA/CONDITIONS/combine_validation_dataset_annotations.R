#start with empty workspace

rm(list = ls(all = TRUE))                                               

# # clear all loaded packages
# invisible(lapply(paste0("package:", names(sessionInfo()$otherPkgs)),
#                  detach,
#                  character.only = TRUE, unload = TRUE))

# turn off scientific notation for plots

options(scipen = 10000)

#### set working directory ####

# here create new folder and set working directory within it

setwd("~/Cel_GRN_manuscript/")

dir.create("input")
dir.create("output")

#### DEFINE FUNCTIONS ####

#### LOAD PACKAGES & FUNCTIONS ####

## First specify the packages of interest

packages <- c("openxlsx")

## Now load or install&load all
package.check <- lapply(packages, function(y){
  
  if (!require(y, character.only = TRUE)) {
    
    install.packages(y, dependencies = TRUE)
    
    library(y, character.only = TRUE)
    
  }
  
})

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


HS_SRA <- read.table("input/HS_SRA.txt",
                     header = TRUE,
                     sep = "\t")


PA14_SRA <- read.table("input/PA14_SRA.txt",
                     header = TRUE,
                     sep = "\t")

daf16_SRA <- read.xlsx("input/IIS_SRA.xlsx",
                       sheet = 1)

daf18_SRA <- read.xlsx("input/IIS_SRA.xlsx",
                       sheet = 2)

daf2_SRA <- read.xlsx("input/IIS_SRA.xlsx",
                      sheet = 3)

daf2daf16_SRA <- read.xlsx("input/IIS_SRA.xlsx",
                           sheet = 4)

daf2daf18_SRA <- read.xlsx("input/IIS_SRA.xlsx",
                           sheet = 5)

all_IIS_SRA <- do.call(rbind, lapply(c("daf16_SRA",
                                   "daf18_SRA",
                                   "daf2_SRA",
                                   "daf2daf16_SRA",
                                   "daf2daf18_SRA"), function(thisgene){
                                     
                                     thisdata <- get(thisgene)
                                     
                                     thisdata[, "genotype"] <- str_remove(thisgene, "_SRA")
                                     
                                     thisdata
                                     
                                   }))

MALE_SRA <- read.table("input/male_SRA.txt",
                     header = TRUE,
                     sep = "\t")

write.xlsx(list("IIS_studies" = all_IIS_SRA,
                "HS_studies" = HS_SRA,
                "PA14_studies" = PA14_SRA,
                "MALE_studies" = MALE_SRA),
           "output/validation_datasets_for_TS6.xlsx")
