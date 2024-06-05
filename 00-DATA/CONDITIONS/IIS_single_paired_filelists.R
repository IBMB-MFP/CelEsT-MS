#start with empty workspace

rm(list = ls(all = TRUE))                                               

# turn off scientific notation for plots

options(scipen = 10000)

#### set working directory ####

# here create new folder and set working directory within it

setwd("~/Cel_GRN_manuscript/")

#### LOAD PACKAGES & FUNCTIONS ####

## First specify the packages of interest

packages <- c("stringr",
              "openxlsx")

## Now load or install&load all
package.check <- lapply(packages, function(y){
  
  if (!require(y, character.only = TRUE)) {
    
    install.packages(y, dependencies = TRUE)
    
    library(y, character.only = TRUE)
    
  }
  
})

#### INPUT DATA ####

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

all_SRA <- rbind(daf16_SRA,
                 daf18_SRA,
                 daf2_SRA,
                 daf2daf16_SRA,
                 daf2daf18_SRA)

#### SPLIT SAMPLES BY SINGLE OR PAIRED ####

write.table(unlist(c(str_split(all_SRA[all_SRA$PairedOrSingle == "SINGLE", "treatment"], ", "), 
                     str_split(all_SRA[all_SRA$PairedOrSingle == "SINGLE", "control"], ", "))),
            "input/IISsinglesamples.txt",
            quote = FALSE,
            col.names = FALSE,
            row.names = FALSE)

write.table(unlist(c(str_split(all_SRA[all_SRA$PairedOrSingle == "PAIRED", "treatment"], ", "), 
                     str_split(all_SRA[all_SRA$PairedOrSingle == "PAIRED", "control"], ", "))),
            "input/IISpairedsamples.txt",
            quote = FALSE,
            col.names = FALSE,
            row.names = FALSE)



