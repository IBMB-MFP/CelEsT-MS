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

packages <- c("stringr",
              "biomaRt",
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

# set biomaRt mart
parasite_mart <- useMart("parasite_mart", dataset = "wbps_gene", host = "https://parasite.wormbase.org", port = 443)

CisBP_TFinfo_withmotif <- readRDS("output/CisBP_TFinfo_withmotif.rds")

# load up output from FIMO on command line. 
# This output file has been placed in the output folder prior to running this script

CisBP_FIMO_full <- read.table("output/fimo_out/fimo.tsv",
                              sep = "\t",
                              header = TRUE)

allTFS_manualtoptargets <- readRDS(file = "output/MODernENCODE_manualtoptargets_operonexcluded.rds")

allmodERN_TFsONLY_BM <- readRDS("output/allmodERN_TFsONLY_BM.rds")

#### ELIMINATE DUPLICATED MOTIFS ####

# we have some duplicated motifs; we got rid of duplicated TFs but some TFs share the same motif
# this occurs because inferred motifs by similarity are often taken from another worm TF. We don't want that
# here we are inferring duplication by the not-so-likely scenario of independent motifs finding the exact same number of targets
# in a minority of cases this will be just coincidence; we will deal with this slightly later
table(CisBP_FIMO_full$motif_id)

# Initially I thought that where we have ChIP-seq data directly for a TF,
# we want to keep the duplicated motif, since it will serve to back up the other line of evidence 
# However I see some problematic cases where oscillatory dynamics seem to bleed through into TFs without oscillatory expression
# due to indirect motif influence. so I want to ELIMINATE all indirect motifs that come from other worm TFs
# how many would I have if I got rid of all indirect motifs?
# get rid of duplicated TFs (motifs based on other worm TFs) that don't have ChIP-seq data

duplication_table <- table(CisBP_FIMO_full$motif_id)[duplicated(table(CisBP_FIMO_full$motif_id))|duplicated(table(CisBP_FIMO_full$motif_id), fromLast = TRUE)]

CisBP_full_separate_by_motif <- lapply(unique(CisBP_FIMO_full$motif_id), function(x){
  
  CisBP_FIMO_full[CisBP_FIMO_full$motif_id == x,  ]
  
})

names(CisBP_full_separate_by_motif) <- unique(CisBP_FIMO_full$motif_id)

# some apparent duplications (same number of targets) are coincidental.
# need a way to ensure I am not sorting out legitimate ones.

duplicated_targets <- unlist(lapply(unique(duplication_table), function(thismany){

  gene_vec <- names(duplication_table[duplication_table == thismany])
  
  gene_target_lists <- lapply(gene_vec, function(thisgene){
    
    thisgene_start <- CisBP_full_separate_by_motif[[thisgene]][, 3:7]
    row.names(thisgene_start) <- NULL
    return(head(thisgene_start, n = 20))
    
  })
  
  output_vec <- duplicated(gene_target_lists) | duplicated(gene_target_lists, fromLast = TRUE)
  names(output_vec) <- gene_vec  
  
  return(output_vec)
  
}))

# figure out the duplicated suites
# this will make a list of all of the sets of genes sharing a motif

dup_suites <- lapply(unique(duplication_table), function(x){
  
  dupsuite_names <- names(duplication_table[duplication_table == x])
  
  dupsuite_names[dupsuite_names %in% names(duplicated_targets)[duplicated_targets]]
  
})

# any suites with spurious duplications? they would show up as having 2 or more separate motifs. 
# there are none. 
sapply(lapply(dup_suites, function(x){
  
  CisBP_TFinfo_withmotif[CisBP_TFinfo_withmotif$wormbase_gseq %in% x, "Motif_ID"]
  
}), function(y){length(unique(y))})

# here we identify duplicated suites where at least one member has a directly determined motif
dup_suites_with_direct <- dup_suites[sapply(dup_suites, function(x){
  
  any(CisBP_TFinfo_withmotif[match(x, CisBP_TFinfo_withmotif$wormbase_gseq), "TF_Status"] == "D")
  
})]

# definitely want to exclude ones with indirect motifs duplicated from direct worm TFs
# on revisit, additionally all of these, whether in MODern or not
# this will produce a vector of ones to exclude, in addition to more to come

EXCLUDE_indirect_in_direct_suite <- unlist(lapply(dup_suites_with_direct, function(x){
  
  theseonesinfo <- CisBP_TFinfo_withmotif[match(x, CisBP_TFinfo_withmotif$wormbase_gseq), ]
  theseones_indirect_names <- theseonesinfo[theseonesinfo$TF_Status == "I", "wormbase_gseq"]
  
  return(theseones_indirect_names)
  
}))

# from the suites with *only* indirect members: 
# pick the one with best similarity regression score and exclude the others. also keep in mind expression level etc.

dup_suites_with_onlyindirect <- dup_suites[sapply(dup_suites, function(x){
  
  all(CisBP_TFinfo_withmotif[match(x, CisBP_TFinfo_withmotif$wormbase_gseq), "TF_Status"] == "I")
  
})]

# lastly pick the best scores for the ones that are not in MODern. 
# This will need to be done manually from the Web Browser because the actual SR scores don't appear in the info file, annoyingly.

# based on similarity regression scores from manual search of CisBP website keep these three from those suites
KEEP_dup_onlyindirect <- c("C33D12.1", # ceh-30, same SR score, arbitrary choice
                           "W05B5.3", #nhr-85 ver nhr-2 based on SR score
                           "ZC247.3", # lin-11 over mec-3 based on SR score alone, check expression etc.
                           "F48G7.11", # nhr-190 over nhrs 133, 184, 179 based on SR score in that order. check expression etc.
                           "F44C8.11",   # nhr-96
                           "T26H2.9",    # nhr-79
                           "F40E10.2")   # sox-3

EXCLUDE_dup_onlyindirect <- unlist(dup_suites_with_onlyindirect)[!unlist(dup_suites_with_onlyindirect) %in% KEEP_dup_onlyindirect]

all_to_exclude <- c(EXCLUDE_dup_onlyindirect,
                    EXCLUDE_indirect_in_direct_suite)

# so 81 TFs will be discarded because they have only non-unique indirect motifs available
# that gives us 260 to work with
# although perhaps not...? 252 seem to come out of the file.

CisBP_FIMO_censored <- CisBP_FIMO_full[!CisBP_FIMO_full$motif_id %in% all_to_exclude, ]
FIMO_TFs <- unique(CisBP_FIMO_censored$motif_id)
FIMO_TFs <- FIMO_TFs[order(FIMO_TFs)]

# add Walhout TFs with >=15 targets
walhoutEV1 <- openxlsx::read.xlsx("~/Downloads/msb167131-sup-0002-datasetev1.xlsx")
# limit to interactions described as high quality
walhoutEV1 <- walhoutEV1[walhoutEV1$`In.high-quality.dataset?` == 'yes', ]
walhoutTFs_cutoff15  <- names(table(walhoutEV1$Prey.sequence.Name)[table(walhoutEV1$Prey.sequence.Name) >= 15])

fullset_TFs <- unique(c(allmodERN_TFsONLY_BM$wormbase_gseq,
                        FIMO_TFs,
                        walhoutEV1$Prey.sequence.Name))

fullset_TFs_BM <- getBM(mart = parasite_mart, 
                              filters = "wormbase_gseqname",
                              value = fullset_TFs,
                              attributes = c("wormbase_gene",
                                             "wormbase_locus",
                                             "wormbase_gseq",
                                             "entrezgene_id",
                                             "description",
                                             "uniprot_swissprot_accession"))

# some duplications of certain TFs, one entry with Uniprot and the other without
# eliminate duplicated ones without Uniprot

fullset_TFs_BM <- fullset_TFs_BM[!((duplicated(fullset_TFs_BM$wormbase_gseq)|
                 duplicated(fullset_TFs_BM$wormbase_gseq, fromLast = TRUE)) &
                 fullset_TFs_BM$uniprot_swissprot_accession == ""), ]

# total set is 596 unique TFs

saveRDS(fullset_TFs_BM,
        "output/fullset_TFs_BM.rds")

# also save as txt file which will be used downstream with SJARACNe as the list of hub genes

write.table(fullset_TFs,
            file = "output/fullsetTFs.txt",
            sep = "\t",
            quote = FALSE,
            row.names = FALSE,
            col.names = FALSE)

#### CALCULATE COMPOSITE SCORE TO ORDER TARGETS ####

#### HERE introduce tweaking / parmeter sweep originally done in benchmarking folder

# FIMO gives an enrichment score
# However also want to take into account homotypic binding
# This is the phenomenon where in many organisms multiple TFBSs ar present in regulated genes
# so a promoter with 5 significant TFBSs is likely regulated by that gene even if non of them is a perfect perfect match
# taking only the best scores as computed by FIMO would discard that richness of information
# so what I will do is compute a score whereby a target gains extra points by virtue of having additional TFBSs in the promoters
# but to stop multiple more spurious sites from overwhelming the genuine perfect matches, compute a score with rapidly diminishing returns for additional binding sites
# I add the scores for additional TFBSs after the best; but divided by 5^(X-1) where X is rank
# i.e. the score for the second TFBSs is added but divided by 5
# the score for the third TFBSs is added but divided by 5 etc.
# parameter of 5 is chosen based on future benchmarking experiments

separate_by_motif <- lapply(unique(CisBP_FIMO_censored$motif_id), function(x){
  
  CisBP_FIMO_censored[CisBP_FIMO_censored$motif_id == x,  ]
  
})

names(separate_by_motif) <- unique(CisBP_FIMO_censored$motif_id)

saveRDS(separate_by_motif,
        "output/separate_by_motif.rds")

sequence_counts_by_motif <- lapply(separate_by_motif, function(x){table(x$sequence_name)})
names(sequence_counts_by_motif) <- unique(CisBP_FIMO_censored$motif_id)

# here looking to see any bias by TF in terms of extent of binding
sites_by_TF <- sapply(sequence_counts_by_motif, table)

percent_not_homotypic <- sapply(sites_by_TF, function(x){
  
  x[1] / sum(x)
  
})

percent_not_homotypic[order(percent_not_homotypic)]

# Note a few with clearly palindromic motifs have unusually high rates of binding - because almost by definition, the appearance of one instance means another on the other strand
write.table(data.frame(percent_not_homotypic, row.names = str_remove(names(percent_not_homotypic), "\\.1$")),
            file = "output/percent_not_homotypic.txt",
            col.names = FALSE)


