#start with empty workspace

rm(list = ls(all = TRUE))                                               

# turn off scientific notation for plots

options(scipen = 10000)

#### set working directory ####

# here create new folder and set working directory within it

setwd("~/Cel_GRN_manuscript/")

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

#### INPUT DATA ####

ENCODE_extraTFs <- read.xlsx("input/encode_extra_TFs.xlsx")

#### write file ####

write.table(ENCODE_extraTFs$bednarrowpeakfile,
            "input/ENCODE_extraTFs_files.txt",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)
