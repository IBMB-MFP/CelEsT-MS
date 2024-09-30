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

packages <- c("openxlsx",
              "dplyr",
              "stringr")

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

github_packages <- c()

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

benchmark_DEstats_RAPToR <- read.table("output/benchmark_DEstats_RAPToR.txt",
                                       header = TRUE)

benchmark_DEstats_RAW <- read.table("output/benchmark_DEstats_RAW.txt",
                                       header = TRUE)

benchmark_obs <- read.table("output/benchmark_observations.txt",
                            sep = "\t",
                            header = TRUE)

if(!file.exists("input/WTF3.xlsx")){
  
  download.file("https://www.embopress.org/doi/full/10.15252/msb.20167131#msb167131-sup-0004",
                destfile = "input/WTF3.xlsx")
  
  
}

wTF3 <- read.xlsx("input/WTF3.xlsx")

#### TF FAMILIES ####

table(wTF3$DBD)[order(table(wTF3$DBD), decreasing = TRUE)]

wTF3[, "family"] <- str_remove(wTF3$DBD, " - [0-9]{1,2} finger.*")
wTF3[, "family"] <- str_remove(wTF3$family, " - [0-9]{1,2} domain.*")

wTF3[, "family"] <- str_remove(wTF3$family, " x[0-9]{1}$")

wTF3[, "family"] <- str_remove(wTF3$family, " $")

wTF3[str_detect(wTF3$DBD, "^WH"), "family"] <- "WH"

wTF3[str_detect(wTF3$DBD, "^AT Hook"), "family"] <- "AT Hook"

wTF3[str_detect(wTF3$DBD, "^HD"), "family"] <- "HD"

wTF3[str_detect(wTF3$DBD, "^ZF C2H2"), "family"] <- "ZF - C2H2"

main_families <- names(table(wTF3$family)[order(table(wTF3$family), decreasing = TRUE)])[1:9]

wTF3[!wTF3$family %in% main_families, "family"] <- "Other"


table(wTF3$family)[order(table(wTF3$family), decreasing = TRUE)]

# totals in Fuxman bass add up to 963 / 941; but its not double count in these families because the only underrepresented family in my count is 'other'...

wTF3[, "bench_no"] <- sapply(wTF3$Sequence.name, function(x){sum(benchmark_obs$target_gseq %in% x)})

saveRDS(wTF3, "output/wTF3_modified.rds")

bench_by_fam <- wTF3 %>% group_by(family) %>% summarise(totalbench = sum(bench_no),
                                                        uniquebench = sum(bench_no != 0),
                                                        count = n())

bench_by_fam[, "rep_frac"] <- (bench_by_fam$uniquebench / bench_by_fam$count) * 100

# ZF CCCH families and T-box absolutely no bench mark knockdowns
# AT Hook has only one; exclude. 
# Remaining 6 families can be assessed. Should annotate with exp, unique and family totals. 

lapply(main_families[!main_families %in% c("ZF - CCCH", "T-box", "AT Hook")] , function(x){
  
temp_obs <- benchmark_obs[benchmark_obs$target_gseq %in% wTF3[wTF3$family == x, "Sequence.name"], ]
temp_bench <- benchmark_DEstats_RAPToR[row.names(temp_obs), ]

write.table(temp_obs,
            paste0("output/", x, "_family_obs.txt"),
            sep = "\t",
            row.names = TRUE,
            col.names = TRUE)

write.table(temp_bench,
            paste0("output/", x, "_family_benchRAPToR.txt"),
            sep = "\t",
            row.names = TRUE,
            col.names = TRUE)

x

})

# and for RAW

lapply(main_families[!main_families %in% c("ZF - CCCH", "T-box", "AT Hook")] , function(x){
  
  temp_obs <- benchmark_obs[benchmark_obs$target_gseq %in% wTF3[wTF3$family == x, "Sequence.name"], ]
  temp_bench <- benchmark_DEstats_RAW[row.names(temp_obs), ]
  
  write.table(temp_obs,
              paste0("output/", x, "_family_obs.txt"),
              sep = "\t",
              row.names = TRUE,
              col.names = TRUE)
  
  write.table(temp_bench,
              paste0("output/", x, "_family_benchRAW.txt"),
              sep = "\t",
              row.names = TRUE,
              col.names = TRUE)
  
  x
  
})



