#start with empty workspace

rm(list = ls(all = TRUE))

# turn off scientific notation for plots

options(scipen = 10000)

#### set working directory ####

# here create new folder and set working directory within it

setwd("~/Cel_GRN_manuscript")

#### DEFINE FUNCTIONS ####

#### LOAD PACKAGES & FUNCTIONS ####

## First specify the packages of interest

packages <- c("biomaRt",
              "readr",
              "openxlsx")

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

andersen_ae2 <- readRDS("output/andersen_ageestimates_CelYA2.rds")

replicate_ages <- andersen_ae2$age.estimates[, "age.estimate"]

# C. elegans CeNDR data (supplementary files to GEO accession number GSE186719)

CeNDR_raw_counts <- read.table(file = "~/NNMT_manuscript/input/GSE186719_Celegans_208strains_609samples_rawCounts.tsv",
                               sep = "\t",
                               fill = TRUE,
                               quote = "",
                               header = TRUE)

CeNDR_normalised_counts_collapse <- readRDS("output/Cel_CeNDR_normalised_counts_gene_level.rds")

CeNDR_normalised_counts_collapse_mostlyexpressed <- readRDS("output/CeNDR_normalised_counts_collapse_mostlyexpressed.rds")

#### PROCESS COUNT DATA AND CORRECT FOR AGE ####

# for age correction and so on we will exclude all genes with any 0 counts
# as fitting splines with missing data points will likely skew the spines excessively.

Cendrcounts_nonzero <- CeNDR_normalised_counts_collapse[!apply(CeNDR_normalised_counts_collapse,
                                                               1, function(x){any(x == 0)}), ]

#### do variance stabilising transformation ####

# convert raw counts into MRN pseudocounts 

# DESeq2 wants a colData object. Not actually used for the normalisation. Here we can use the sample IDs with tissue.
CeNDR_col_data <- str_extract(colnames(CeNDR_raw_counts)[2:ncol(CeNDR_raw_counts)], "^[0-9A-Z]+")
CeNDR_col_data <- as.matrix(CeNDR_col_data)

rownames(CeNDR_col_data) <- colnames(CeNDR_raw_counts)[2:ncol(CeNDR_raw_counts)]
colnames(CeNDR_col_data) <- c("Strain")

# convert to matrix of integers

CeNDR_raw_counts_mat <- as.matrix(CeNDR_raw_counts[, 2:ncol(CeNDR_raw_counts)])
row.names(CeNDR_raw_counts_mat) <- CeNDR_raw_counts[, 1]
colnames(CeNDR_raw_counts_mat) <- colnames(CeNDR_raw_counts)[2:ncol(CeNDR_raw_counts)]

CeNDR_raw_counts_mat <- data.frame(CeNDR_raw_counts_mat)
CeNDR_raw_counts_mat[, "gene"] <- str_extract(row.names(CeNDR_raw_counts_mat), "^[0-9A-Z]+\\.[0-9]{1,2}")

# collapse to gene level by summing transcript counts
CeNDR_raw_counts_mat_collapse <- stats::aggregate(. ~ gene, data = CeNDR_raw_counts_mat, FUN = base::sum)

row.names(CeNDR_raw_counts_mat_collapse) <- CeNDR_raw_counts_mat_collapse[, 1]
CeNDR_raw_counts_mat_collapse <- CeNDR_raw_counts_mat_collapse[, 2:ncol(CeNDR_raw_counts_mat_collapse)]

# some raw counts are not integers. need to be rounded
CeNDR_raw_counts_mat_collapse <- round(CeNDR_raw_counts_mat_collapse)

# need to create a DESeq2 object. Design set to ~1 allows for use of estimateSizeFactors
tempdds <- DESeqDataSetFromMatrix(countData = CeNDR_raw_counts_mat_collapse, colData = CeNDR_col_data, design = ~ 1)

vst <- vst(tempdds, blind = TRUE)
vst_assay <- assay(vst)

write.table(t(vst_assay),
            file = "output/vst_cendr.txt",
            sep = "\t",
            row.names = TRUE,
            col.names = TRUE)

vst_age <- t(vst_assay)

vst_age <- vst_age[, colnames(vst_age) %in% row.names(Cendrcounts_nonzero)]
vst_age <- as.data.frame(vst_age)
vst_age[, "age"] <- andersen_ae2$age.estimates[match(row.names(vst_age), row.names(andersen_ae2$age.estimates)), "age.estimate"]

age_corrected_nonzero_counts <- apply(vst_age[1:(ncol(vst_age)-1)], 2, function(thisgene){
  
  logspline <- smooth.spline(x = vst_age$age, y = thisgene, df = 6)
  # predict(thisgene_spline_df6, x = replicate_ages)$y
  # plot(x = replicate_ages, y = log10(thisgenevec[names(replicate_ages)]), pch = 20)
  # points(x = replicate_ages, y = predict(thisgene_spline_df6, x = replicate_ages)$y, col = "red")
  
  return(residuals(logspline))
  
})

row.names(age_corrected_nonzero_counts) <- row.names(vst_age)
age_corrected_nonzero_counts_df <- as.data.frame(t(age_corrected_nonzero_counts))

# in our case we will fill both columns with the gene and ensure that they have the correct column headers
age_corrected_vst_ARACNe <- cbind(row.names(age_corrected_nonzero_counts_df),
                                  row.names(age_corrected_nonzero_counts_df),
                                  age_corrected_nonzero_counts_df)

# the columns have to be headed exactly like this NB capital letters
colnames(age_corrected_vst_ARACNe)[1:2] <- c("isoformId",
                                             "geneSymbol")

write.table(age_corrected_vst_ARACNe,
            file = "output/CeNDR_age_corrected_vst_for_ARACNe.txt",
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE,
            quote = FALSE)

# filter TFs by ones presnt in non-zero counts; 123 TFs drop out at this stage
TFs_for_ARACNe <- fullset_TFs_BM$wormbase_gseq[fullset_TFs_BM$wormbase_gseq %in% row.names(age_corrected_vst_ARACNe)]

write.table(TFs_for_ARACNe,
            "output/TFs_for_agecorr_vst_ARACNe.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)

# create output directories for running both the age-corrected counts and the raw pseudocounts 

dir.create("output/SJARACNe")
dir.create("output/SJARACNe/age_corrected_non_zero")
