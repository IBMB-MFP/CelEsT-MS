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

packages <- c("stringr",
              "UniprotR")

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

biocmanager_packages <- c("reshape2") 

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

# Wormbase descriptions

download.file(url = "https://downloads.wormbase.org/releases/WS288/species/c_elegans/PRJNA13758/annotation/c_elegans.PRJNA13758.WS288.functional_descriptions.txt.gz",
              "input/c_elegans.PJRNA13758.WS288.functional_descriptions.txt.gz")

Wormbase_descriptions <- readLines("input/c_elegans.PJRNA13758.WS288.functional_descriptions.txt.gz")

Wormbase_descriptions_list <- lapply(unlist(str_split(paste(Wormbase_descriptions, collapse = " "), " = ")), function(x){
  
  thisvec <-  unlist(str_split(x, "\t"))
  
  deconstructed_thisvec <- unlist(str_split(unlist(str_split(unlist(str_split(thisvec, " Concise description: ")), " Automated description: ")), " Gene class description: "))
  
  names(deconstructed_thisvec) <- c("WBGeneID",
                                    "gene_name",
                                    "gseq",
                                    "concise_description",
                                    "automated_description",
                                    "gene_class")
  
  return(deconstructed_thisvec)
  
})

Wormbase_descriptions_df <- data.frame(do.call(rbind, Wormbase_descriptions_list))
Wormbase_descriptions_df[1, 1] <- str_extract(Wormbase_descriptions_df[1, 1], "WBGene.*$")

Wormbase_descriptions_df$gseq <- str_remove(Wormbase_descriptions_df$gseq, "[a-z]+$")

saveRDS(Wormbase_descriptions_df,
        "output/Wormbase_descriptions_df.rds")


Wormbase_descriptions_ourTFS <- Wormbase_descriptions_df[str_remove(Wormbase_descriptions_df$gseq, "[a-z]{1}$") %in% fullset_TFs_BM$wormbase_gseq, ]

Wormbase_negative_regulators <- Wormbase_descriptions_ourTFS[str_detect(Wormbase_descriptions_ourTFS$automated_description, "negative regulation of transcription") & !str_detect(Wormbase_descriptions_ourTFS$automated_description, "positive regulation of transcription"), ]
Wormbase_positive_regulators <- Wormbase_descriptions_ourTFS[str_detect(Wormbase_descriptions_ourTFS$automated_description, "positive regulation of transcription") & !str_detect(Wormbase_descriptions_ourTFS$automated_description, "negative regulation of transcription"), ]

Wormbase_bidirectional_regulators <- Wormbase_descriptions_ourTFS[str_detect(Wormbase_descriptions_ourTFS$automated_description, "negative regulation of transcription") & str_detect(Wormbase_descriptions_ourTFS$automated_description, "positive regulation of transcription"), ]

## Add Uniprot description directionality
getOption('timeout')
options(timeout = 200)
PF <- GetProteinFunction(fullset_TFs_BM$uniprot_swissprot_accession[fullset_TFs_BM$uniprot_swissprot_accession != ""])

write.xlsx(PF,
           "output/Uniprot_proteinfunc_raw.xlsx",
           rowNames = TRUE)

PF_reviewed <- read.xlsx("output/Uniprot_proteinfunc_REVIEWED.xlsx")

Uniprot_negative_regulators <- fullset_TFs_BM[match(PF_reviewed[PF_reviewed$manual_activity == -1, "UNIPROT"], fullset_TFs_BM$uniprot_swissprot_accession), "wormbase_gseq"]
Uniprot_negative_regulators <- Uniprot_negative_regulators[!is.na(Uniprot_negative_regulators)]

Uniprot_positive_regulators <- fullset_TFs_BM[match(PF_reviewed[PF_reviewed$manual_activity == 1, "UNIPROT"], fullset_TFs_BM$uniprot_swissprot_accession), "wormbase_gseq"]
Uniprot_positive_regulators <- Uniprot_positive_regulators[!is.na(Uniprot_positive_regulators)]

# check for any default conflicts. None found
base::intersect(Uniprot_negative_regulators, Wormbase_positive_regulators$gseq)
base::intersect(Uniprot_positive_regulators, Wormbase_negative_regulators$gseq)

fullset_TFs_BM[, "Wormbase_mode_of_regulation"] <- NA
fullset_TFs_BM[fullset_TFs_BM$wormbase_gseq %in% Wormbase_positive_regulators$gseq, "Wormbase_mode_of_regulation"] <- 1
fullset_TFs_BM[fullset_TFs_BM$wormbase_gseq %in% Wormbase_negative_regulators$gseq, "Wormbase_mode_of_regulation"] <- -1
fullset_TFs_BM[fullset_TFs_BM$wormbase_gseq %in% Wormbase_bidirectional_regulators$gseq, "Wormbase_mode_of_regulation"] <- 0

fullset_TFs_BM[, 'Uniprot_mode_of_regulation'] <- NA
fullset_TFs_BM[fullset_TFs_BM$uniprot_swissprot_accession %in% PF_reviewed[PF_reviewed$manual_activity == 1, "UNIPROT"], 'Uniprot_mode_of_regulation'] <- 1
fullset_TFs_BM[fullset_TFs_BM$uniprot_swissprot_accession %in% PF_reviewed[PF_reviewed$manual_activity == 0, "UNIPROT"], 'Uniprot_mode_of_regulation'] <- 0
fullset_TFs_BM[fullset_TFs_BM$uniprot_swissprot_accession %in% PF_reviewed[PF_reviewed$manual_activity == -1, "UNIPROT"], 'Uniprot_mode_of_regulation'] <- -1

write.xlsx(fullset_TFs_BM,
           "output/fullset_TFs_BM_with_MOR_for_TS2.xlsx")

