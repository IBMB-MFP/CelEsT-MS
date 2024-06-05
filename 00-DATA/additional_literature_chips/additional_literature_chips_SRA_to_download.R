#start with empty workspace

rm(list = ls(all = TRUE))

setwd("~/Cel_GRN_manuscript/")

#### LOAD PACKAGES & FUNCTIONS ####

## First specify the packages of interest

packages <- c("stringr")

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

#### make file of samples to download ####

simple <- read.table("input/additional_literature_chips_simple.txt",
                     sep = "\t",
                     header = TRUE)

complex <- read.table("input/additional_literature_chips_complex.txt",
                      sep = "\t",
                      header = TRUE)

single_sample_vec <- unlist(c(str_split(simple$target_samples, pattern = ", "),
                              str_split(simple$input_samples, pattern = ", "),
                              str_split(complex[complex$paired_or_single == "SINGLE", "target_sample1_files"], pattern = ", "),
                              str_split(complex[complex$paired_or_single == "SINGLE", "input_sample1_files"], pattern = ", "),
                              str_split(complex[complex$paired_or_single == "SINGLE", "target_sample2_files"], pattern = ", "),
                              str_split(complex[complex$paired_or_single == "SINGLE", "input_sample2_files"], pattern = ", ")))

paired_sample_vec <- unlist(c(str_split(complex[complex$paired_or_single == "PAIRED", "target_sample1_files"], pattern = ", "),
                       str_split(complex[complex$paired_or_single == "PAIRED", "input_sample1_files"], pattern = ", "),
                       str_split(complex[complex$paired_or_single == "PAIRED", "target_sample2_files"], pattern = ", "),
                       str_split(complex[complex$paired_or_single == "PAIRED", "input_sample2_files"], pattern = ", ")))

paired_sample_vec <- paired_sample_vec[paired_sample_vec != ""]                     

write.table(single_sample_vec,
            "output/additional_literature_chips_SRA_download_single.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)      

write.table(paired_sample_vec,
            "output/additional_literature_chips_SRA_download_paired.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)  
                     
                     