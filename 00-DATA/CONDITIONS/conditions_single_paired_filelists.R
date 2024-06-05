#start with empty workspace

rm(list = ls(all = TRUE))                                               

# turn off scientific notation for plots

options(scipen = 10000)

#### set working directory ####

# here create new folder and set working directory within it

setwd("~/Cel_GRN_manuscript/")

#### LOAD PACKAGES & FUNCTIONS ####

## First specify the packages of interest

packages <- c("stringr")

## Now load or install&load all
package.check <- lapply(packages, function(y){
  
  if (!require(y, character.only = TRUE)) {
    
    install.packages(y, dependencies = TRUE)
    
    library(y, character.only = TRUE)
    
  }
  
})

#### INPUT DATA ####
getwd()
males_SRA <- read.table("input/male_SRA.txt",
           sep = "\t",
           header = TRUE,
           fill = TRUE)

PA14_SRA <- read.table("input/PA14_SRA.txt",
                        sep = "\t",
                        header = TRUE,
                        fill = TRUE)

HS_SRA <- read.table("input/HS_SRA.txt",
                       sep = "\t",
                       header = TRUE,
                       fill = TRUE)

all_SRA <- rbind(males_SRA,
                 PA14_SRA,
                 HS_SRA)

all_SRA <- all_SRA[!all_SRA$control == "", ]

#### SPLIT SAMPLES BY SINGLE OR PAIRED ####

write.table(unlist(c(str_split(all_SRA[all_SRA$PairedOrSingle == "SINGLE", "treatment"], ", "), 
                     str_split(all_SRA[all_SRA$PairedOrSingle == "SINGLE", "control"], ", "))),
            "input/conditionssinglesamples.txt",
            quote = FALSE,
            col.names = FALSE,
            row.names = FALSE)

write.table(unlist(c(str_split(all_SRA[all_SRA$PairedOrSingle == "PAIRED", "treatment"], ", "), 
                     str_split(all_SRA[all_SRA$PairedOrSingle == "PAIRED", "control"], ", "))),
            "input/conditionspairedsamples.txt",
            quote = FALSE,
            col.names = FALSE,
            row.names = FALSE)



