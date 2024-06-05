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

setwd("~/Cel_GRN_orthology")

#### DEFINE FUNCTIONS ####

#### LOAD PACKAGES & FUNCTIONS ####

## First specify the packages of interest

packages <- c("stringr",
              "biomaRt")

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

CisBP_TFinfo_withmotif <- readRDS("~/Cel_GRN_assembly/output/CisBP_TFinfo_withmotif.rds")

# load up output from FIMO on command line. 
# This output file has been placed in the output folder prior to running this script

#### test on opinata ####

CisBP_FIMO_full_inopinata <- read.table("output/FIMO_output/inopinata/fimo.tsv",
                              sep = "\t",
                              header = TRUE)

sum(CisBP_FIMO_full_inopinata$q.value < 0.1) / length(unique(CisBP_FIMO_full_inopinata$motif_id))

library(dplyr)

by_motif <- lapply(unique(CisBP_FIMO_full_inopinata$motif_id), function(x){
  
  CisBP_FIMO_full_inopinata[CisBP_FIMO_full_inopinata$motif_id == x, ]
  
})

names(by_motif) <- unique(CisBP_FIMO_full_inopinata$motif_id)

qvalpointtwo <- sapply(by_motif, function(x){
  
  sum(x$q.value < 0.2)
  
})

# average of 350 significant with q-value (/FDR) 0.2
# yeah so for the most part, it's none at all. but some 13000!!! 
# I have a feeling the q values are not computed by motif. I could try that, couldnt I.
# they do seem to be done by motif, but I suppose the information content is not good enough. 
# so lets do the best hits. 

#### ELIMINATE DUPLICATED MOTIFS ####

# we have some duplicated motifs; we got rid of duplicated TFs but some TFs share the same motif
# this occurs because inferred motifs by similarity are often taken from another worm TF. We don't want that
# here we are inferring duplication by the not-so-likely scenario of independent motifs finding the exact same number of targets
# in a minority of cases this will be just coincidence; we will deal with this slightly later
table(CisBP_FIMO_full_inopinata$motif_id)

# Initially I thought that where we have ChIP-seq data directly for a TF,
# we want to keep the duplicated motif, since it will serve to back up the other line of evidence 
# However I see some problematic cases where oscillatory dynamics seem to bleed through into TFs without oscillatory expression
# due to indirect motif influence. so I want to ELIMINATE all indirect motifs that come from other worm TFs
# how many would I have if I got rid of all indirect motifs?
# get rid of duplicated TFs (motifs based on other worm TFs) that don't have ChIP-seq data

duplication_table <- table(CisBP_FIMO_full_inopinata$motif_id)[duplicated(table(CisBP_FIMO_full_inopinata$motif_id))|duplicated(table(CisBP_FIMO_full_inopinata$motif_id), fromLast = TRUE)]
duplication_table[order(duplication_table)]

CisBP_full_separate_by_motif <- lapply(unique(CisBP_FIMO_full_inopinata$motif_id), function(x){
  
  CisBP_FIMO_full_inopinata[CisBP_FIMO_full_inopinata$motif_id == x,  ]
  
})

names(CisBP_full_separate_by_motif) <- unique(CisBP_FIMO_full_inopinata$motif_id)

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

# from the suites with *only* indirect members, keep if in MODern; 
# otherwise pick the one with best similarity regression score and exclude the others
# this seems ok on revisiting; because with only indirect members
dup_suites_with_onlyindirect <- dup_suites[sapply(dup_suites, function(x){
  
  all(CisBP_TFinfo_withmotif[match(x, CisBP_TFinfo_withmotif$wormbase_gseq), "TF_Status"] == "I")
  
})]

# lastly pick the best scores for the ones that are not in MODern. 
# This will need to be done manually from the Web Browser because the actual SR scores don't appear in the info file, annoyingly.

# based on similarity regression scores from manual search of CisBP website keep these three from those suites
KEEP_dup_onlyindirect <- c("C33D12.1", # ceh-31; since ceh-30 is possibly repressed in hermaphrodites 
  "ZC247.3", # slightly better score
  "F44C8.8", # nhr-133; nhr-190 has a slightly better score here but it much more neuronal specific
  "K08A8.2", # sox-2; sox-3 has slightly better score but ox-2 is expressed in more structures and is more highly expressed
  "F44C8.11",   # nhr-96
"T26H2.9"    # nhr-79
)   

EXCLUDE_dup_onlyindirect <- unlist(dup_suites_with_onlyindirect)[!unlist(dup_suites_with_onlyindirect) %in% KEEP_dup_onlyindirect]

# put the exclusions together. To sum up;
# exclude ones that are duplicated with other worm TFs, where all have indirect motifs but others with that motif are in MODern;
# exclude ones that are duplicated with other worm TFs, where all have indirect motifs, none in MODern but these don't have the best SR scores;
# exclude ones that are duplicated with other worm TFs that have direct motifs (UNLESS they are in MODern)

all_to_exclude <- c(EXCLUDE_dup_onlyindirect,
                    EXCLUDE_indirect_in_direct_suite)

# so 82! TFs will be discarded because they have only non-unique indirect motifs available
# that gives us 263 to work with
# although perhaps not...? 255 seem to come out of the file.

CisBP_FIMO_censored <- CisBP_FIMO_full_inopinata[!CisBP_FIMO_full_inopinata$motif_id %in% all_to_exclude, ]
# FIMO_TFs <- unique(CisBP_FIMO_censored$motif_id)
# FIMO_TFs <- FIMO_TFs[order(FIMO_TFs)]

#### let's try to look at overlaps with elegans ####

CisBP_FIMO_full_elegans <- read.table("~/Cel_GRN_assembly/output/fimo_out/fimo.tsv",
                              sep = "\t",
                              header = TRUE)

elegans_by_motif <- lapply(unique(CisBP_FIMO_full_elegans$motif_id), function(x){
  
  CisBP_FIMO_full_elegans[CisBP_FIMO_full_elegans$motif_id == x,  ]
  
})

names(elegans_by_motif) <- unique(CisBP_FIMO_full_elegans$motif_id)

# pick top 2000; although it looks like homotypic binding is not so great to do,
# we could use extra binding sites here to help with the cutoff.

# get one 2 one homology to restrict this.

Cinop_one2one <- readRDS("output/Cinopinata_one2one_promoterseq.rds")

elegans_by_motif_one2oneinop <- lapply(elegans_by_motif, function(x){

  x[x$sequence_name %in% Cinop_one2one$wormbase_gseq, ]

})

jacc_motif <- sapply(intersect(names(by_motif), names(elegans_by_motif)), function(thismotif){
thismotif <- intersect(names(by_motif), names(elegans_by_motif))[1]
  
  inop_ex <- by_motif[[thismotif]][order(by_motif[[thismotif]]$p.value, decreasing = FALSE), ][1:2000, ]
  
  elegans_ex <- elegans_by_motif_one2oneinop[[thismotif]][order(elegans_by_motif_one2oneinop[[thismotif]]$p.value, decreasing = FALSE), ][1:2000, ]
  
  jaccard(Cinop_one2one[match(inop_ex$sequence_name, Cinop_one2one$cainopprjdb5687_gene), "wormbase_gseq"],
          elegans_ex$sequence_name)
  
})

plot(jacc_motif)

random_distro <- sapply(intersect(names(by_motif), names(elegans_by_motif)), function(thismotif){
  
inop_no <- min(nrow(by_motif[[thismotif]]), 2000)
elegans_no <- min(nrow(elegans_by_motif[[thismotif]]), 2000)
  
jaccard(sample(Cinop_one2one$wormbase_gseq, size = inop_no),
sample(Cinop_one2one$wormbase_gseq, size = elegans_no))
  
})



# there is clearly less overlap than random...

jacc_for_plot <- rbind(data.frame(jacc = jacc_motif[names(jacc_motif) %in% Cinop_one2one$wormbase_gseq], 
                 group = "one2one_orthologues"),
      data.frame(jacc = jacc_motif[!names(jacc_motif) %in% Cinop_one2one$wormbase_gseq], 
                 group = "notone2one_orthologues"),
      data.frame(jacc = random_distro, 
                 group = "random"))
library(ggplot2)
ggplot(jacc_for_plot, aes(x = group, y = jacc)) + 
  geom_boxplot() + 
  geom_jitter(width = 0.1) + 
  theme_classic()

boxplot(jacc_motif[names(jacc_motif) %in% Cinop_one2one$wormbase_gseq],
     jacc_motif[!names(jacc_motif) %in% Cinop_one2one$wormbase_gseq],
     random_distro)

jaccard <- function(a, b) {
  intersection = length(intersect(a, b))
  union = length(a) + length(b) - intersection
  return (intersection/union)
}



# will want to test for conservation of targets. 
# at the level of species and at the level of TFs.
# would want to do with all species data. For TFs, 


#### quick look at another species ####
