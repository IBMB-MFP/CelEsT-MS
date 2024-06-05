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

benchmark_SRA_all <- read.table("input/benchmark_SRA_all.txt",
           sep = "\t",
           header = TRUE,
           fill = TRUE)

#### SPLIT SAMPLES BY SINGLE OR PAIRED ####

write.table(unlist(c(str_split(benchmark_SRA_all[benchmark_SRA_all$single_or_paired == "SINGLE", "treatment_samples"], ", "), 
                     str_split(benchmark_SRA_all[benchmark_SRA_all$single_or_paired == "SINGLE", "control_samples"], ", "))),
            "input/benchmarksinglesamples.txt",
            quote = FALSE,
            col.names = FALSE,
            row.names = FALSE)

write.table(unlist(c(str_split(benchmark_SRA_all[benchmark_SRA_all$single_or_paired == "PAIRED", "treatment_samples"], ", "), 
                     str_split(benchmark_SRA_all[benchmark_SRA_all$single_or_paired == "PAIRED", "control_samples"], ", "))),
            "input/benchmarkpairedsamples.txt",
            quote = FALSE,
            col.names = FALSE,
            row.names = FALSE)



