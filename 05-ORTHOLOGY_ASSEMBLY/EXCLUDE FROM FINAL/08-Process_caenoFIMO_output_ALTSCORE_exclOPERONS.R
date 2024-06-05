#start with empty workspace

rm(list = ls(all = TRUE))

# clear all loaded packages
# invisible(lapply(paste0("package:", names(sessionInfo()$otherPkgs)),
#                  detach,
#                  character.only = TRUE, unload = TRUE))

# turn off scientific notation for plots

options(scipen=10000)

#### set working directory ####

# here create new folder and set working directory within it

setwd("~/Cel_GRN_orthology")

#### DEFINE FUNCTIONS ####

ecdf_fun <- function(x,perc) ecdf(x)(perc)

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

github_packages <- c("r-lib/conflicted",
                     "etam4260/kneedle") # kneedle to automatially find elbow, if that works

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

elegans_FIMO_full <- read.table("~/Cel_GRN_assembly/output/fimo_out/fimo.tsv",
                              sep = "\t",
                              header = TRUE)

caeno_FIMO_sp <- list.files("output/FIMO_output/")

caeno_FIMO_list <- lapply(caeno_FIMO_sp, function(x){

  message(paste0("Reading FIMO output for ", x))
  
  read.table(paste0("output/FIMO_output/", x, "/fimo.tsv"),
             sep = "\t",
             header = TRUE)
  
})

names(caeno_FIMO_list) <- caeno_FIMO_sp

#### DEFINE OPERON LIST FOR EXCLUSION ####

# elegans predicted operons from supplementary 3 from blumenthal paper [Allen et al. Genome Research 2011 10.1101/gr.113811.110]
blumenthal_s3_sh2exactmatch <- read_excel("~/Downloads/Supplemental_File3.xls",
                                          sheet = 2,
                                          col_names = FALSE)

blumenthal_s3_sh3partmatch <- read_excel("~/Downloads/Supplemental_File3.xls",
                                         sheet = 3,
                                         col_names = TRUE)

blumenthal_s3_sh5notmatch <- read_excel("~/Downloads/Supplemental_File3.xls",
                                        sheet = 5,
                                        col_names = FALSE)

exactmatch_operon_list <- str_split(str_remove_all(str_remove_all(str_remove_all(blumenthal_s3_sh2exactmatch$...2, "\\["), "\\]"), "\\'"), ", ")
partmatch_operon_list <- str_split(str_remove_all(str_remove_all(str_remove_all(  blumenthal_s3_sh3partmatch$`predicted genes`, "\\["), "\\]"), "\\'"), ", ")
nomatch_operon_list <- str_split(str_remove_all(str_remove_all(str_remove_all(blumenthal_s3_sh5notmatch$...1, "\\["), "\\]"), "\\'"), ", ")

operon_list <- c(exactmatch_operon_list,
                 partmatch_operon_list,
                 nomatch_operon_list)

downstream_operon_genes <- unlist(lapply(operon_list, function(x){
  
  x[-1]
  
}))

# operons in briggsae

# Table S7 from Jhaveri et al. 2023 (G3, 10.1093/g3journal/jkac101)
briggsae_confirmed_operons <- read.xlsx("~/Cel_GRN_orthology/input/Jhaveri_et_al_2022_briggsaeoperons_TableS7.xlsx",
                                        sheet = 1)

colnames(briggsae_confirmed_operons) <- briggsae_confirmed_operons[6, ]
briggsae_confirmed_operons <- briggsae_confirmed_operons[7:nrow(briggsae_confirmed_operons), ]

# remove hash which indicates additional experimental verification for a single gene
briggsae_confirmed_operons$`gene 1` <- str_remove(briggsae_confirmed_operons$`gene 1`, "#")

# convert to operon list

briggsae_conf_oper_list <- apply(briggsae_confirmed_operons, 1, function(x){
  
  x[str_detect(names(x), "^gene")]
  
})

briggsae_first_genes <- briggsae_conf_oper_list[1, ]
briggsae_downstream_genes <- unlist(briggsae_conf_oper_list)[!unlist(briggsae_conf_oper_list) %in% briggsae_first_genes]
briggsae_downstream_genes <- briggsae_downstream_genes[!is.na(briggsae_downstream_genes)]

# get WB Gene IDs for these. They seem to all be returned (plus some extras?)
# some few duplications occur in the table but they are due to two different entrezgene IDs for a few genes
# wormbase identifiers map one 2 one

briggsae_op_downstr_BM <- getBM(mart = parasite_mart,
                                values = briggsae_downstream_genes,
                                filters = "wormbase_gseqname",
                                attributes = c("wormbase_gseq",
                                               "entrezgene_id",
                                               "wormbase_gene",
                                               "wormbase_locus",
                                               "caelegprjna13758_gene",
                                               "caelegprjna13758_orthology_type"
                                ))

briggsae_op_downstr_elegWBID <- briggsae_op_downstr_BM[briggsae_op_downstr_BM$caelegprjna13758_orthology_type == "ortholog_one2one", "caelegprjna13758_gene"]

briggsae_op_downstr_eleg_gseq <- getBM(mart = parasite_mart,
                                       values = briggsae_op_downstr_elegWBID,
                                       filters = "wbps_gene_id",
                                       attributes = c("wormbase_gseq",
                                                      "wormbase_gene"))

briggsae_op_downstr_eleg_gseq <- unique(briggsae_op_downstr_eleg_gseq$wormbase_gseq)

all_downstream_operon_genes <- unique(c(briggsae_op_downstr_eleg_gseq, downstream_operon_genes))

#### ELIMINATE DUPLICATED MOTIFS ####

# we have some duplicated motifs; we got rid of duplicated TFs but some TFs share the same motif
# this occurs because inferred motifs by similarity are often taken from another worm TF. We don't want that
# here we are inferring duplication by the not-so-likely scenario of independent motifs finding the exact same number of targets
# in a minority of cases this will be just coincidence; we will deal with this slightly later
table(elegans_FIMO_full$motif_id)

# Initially I thought that where we have ChIP-seq data directly for a TF,
# we want to keep the duplicated motif, since it will serve to back up the other line of evidence 
# However I see some problematic cases where oscillatory dynamics seem to bleed through into TFs without oscillatory expression
# due to indirect motif influence. so I want to ELIMINATE all indirect motifs that come from other worm TFs
# how many would I have if I got rid of all indirect motifs?
# get rid of duplicated TFs (motifs based on other worm TFs) that don't have ChIP-seq data

duplication_table <- table(elegans_FIMO_full$motif_id)[duplicated(table(elegans_FIMO_full$motif_id))|duplicated(table(elegans_FIMO_full$motif_id), fromLast = TRUE)]
duplication_table[order(duplication_table)]

CisBP_full_separate_by_motif <- lapply(unique(elegans_FIMO_full$motif_id), function(x){
  
  elegans_FIMO_full[elegans_FIMO_full$motif_id == x,  ]
  
})

names(CisBP_full_separate_by_motif) <- unique(elegans_FIMO_full$motif_id)

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
# this will make a list of all of the sets of genes likely sharing a motif

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

# from the suites with *only* indirect members: 
# pick the one with best similarity regression score and exclude the others. also keep in mind expression level etc.

dup_suites_with_onlyindirect <- dup_suites[sapply(dup_suites, function(x){
  
  all(CisBP_TFinfo_withmotif[match(x, CisBP_TFinfo_withmotif$wormbase_gseq), "TF_Status"] == "I")
  
})]
# 
# EXCLUDE_indirect_onlyindirectsuite_otherMODern <- unlist(lapply(dup_suites_with_onlyindirect, function(x){
#   
#   if(any(x %in% names(allTFS_manualtoptargets))){
#     
#     return(x[!x %in% names(allTFS_manualtoptargets)])
#     
#   } else { 
#     
#     return(NULL)
#     
#   }
#   
# }))

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

# 81 exclusions. 252 here.

elegans_FIMO_censored <- elegans_FIMO_full[!elegans_FIMO_full$motif_id %in% all_to_exclude, ]
elegans_FIMO_TFs <- unique(elegans_FIMO_censored$motif_id)
elegans_FIMO_TFs <- elegans_FIMO_TFs[order(elegans_FIMO_TFs)]

#### PROCESS TARGETS ####

separate_by_motif <- lapply(unique(elegans_FIMO_censored$motif_id), function(x){
  
  elegans_FIMO_censored[elegans_FIMO_censored$motif_id == x,  ]
  
})

names(separate_by_motif) <- unique(elegans_FIMO_censored$motif_id)

## Separate by motifs for different species

caeno_sp_by_motif <- lapply(caeno_FIMO_list, function(thissp_FIMO){
  
  temp_by_motif <- lapply(unique(elegans_FIMO_censored$motif_id), function(x){
    
    thissp_FIMO[thissp_FIMO$motif_id == x,  ]
    
  })
  
  names(temp_by_motif) <- unique(elegans_FIMO_censored$motif_id)
  
  temp_by_motif
  
})

names(caeno_sp_by_motif) <- names(caeno_FIMO_list)

#### calc score ####

# to start off we will not use any homotypic binding parameters

calc_score_by_motif <- lapply(separate_by_motif, function(x){

  sapply(unique(x$sequence_name), function(y){

    max(x[x$sequence_name == y, "score"])
    
  })
  
})

names(calc_score_by_motif) <- names(separate_by_motif)

calc_score_by_motif_order <- lapply(calc_score_by_motif, function(x){
  
  x[order(x, decreasing = TRUE)]
  
})

names(calc_score_by_motif_order) <- names(calc_score_by_motif)

caeno_calc_score_list <- lapply(caeno_sp_by_motif, function(thisFIMO_list){
  
  temp_calc_score_by_motif <- lapply(thisFIMO_list, function(x){
    
    sapply(unique(x$sequence_name), function(y){
      
      max(x[x$sequence_name == y, "score"])
      
    })
    
  })
  
  names(temp_calc_score_by_motif) <- names(thisFIMO_list)
  
  temp_calc_score_by_motif_order <- lapply(temp_calc_score_by_motif, function(x){
    
    x[order(x, decreasing = TRUE)]
    
  })
  
  names(temp_calc_score_by_motif_order) <- names(temp_calc_score_by_motif)
  
  temp_calc_score_by_motif_order
  
})

names(caeno_calc_score_list) <- names(caeno_sp_by_motif)

#### BUILD ARRAY OF TARGETS ####

elegans_targets <- unique(elegans_FIMO_censored$sequence_name)

parasites_BM_list <- readRDS("output/parasites_BM_list.rds")

elegans_target_table <- as.data.frame(lapply(calc_score_by_motif, function(thismotif_hits){
  
  tempvec <- vector(length = length(elegans_targets))
  names(tempvec) <- elegans_targets
  tempvec[] <- 0
  
  tempvec[match(names(thismotif_hits), names(tempvec))] <- thismotif_hits
  
  tempvec
  
}))

caeno_target_tables <- lapply(names(caeno_calc_score_list), function(thisspecies){

  thisspecies_scores <- caeno_calc_score_list[[thisspecies]]
  
  thissp_one2one_object <- readRDS(paste0("output/C", thisspecies, "_one2one_promoterseq.rds"))
  
  thisspecies_scores_ortholog <- as.data.frame(lapply(thisspecies_scores,  function(thisspthisTF){
    
    tempvec <- vector(length = length(elegans_targets))
    names(tempvec) <- elegans_targets
    tempvec[] <- 0
    
    thissp_thisTF_gseqnames <- thissp_one2one_object[match(names(thisspthisTF), thissp_one2one_object[, 5]), "wormbase_gseq"]
    
    thissp_thisTF_scores <- thisspthisTF[thissp_thisTF_gseqnames %in% elegans_targets]
    names(thissp_thisTF_scores) <- thissp_thisTF_gseqnames[thissp_thisTF_gseqnames %in% elegans_targets]

    tempvec[names(thissp_thisTF_scores)] <- thissp_thisTF_scores
    
    # make NA for ones that are not in the one2one ortholog BM and so were not considered for inclusion in the study
    
    tempvec[!names(tempvec) %in% thissp_one2one_object$wormbase_gseq] <- NA

    tempvec
    
  }))
  
  # filter by TFs with one to one
  thisspecies_scores_ortholog[, !colnames(thisspecies_scores_ortholog) %in% parasites_BM_list[[thisspecies]][parasites_BM_list[[thisspecies]][, 6] == "ortholog_one2one", "wormbase_gseq"]] <- NA
  
  thisspecies_scores_ortholog
    
  })

names(caeno_target_tables) <- caeno_FIMO_sp

caeno_target_array <- array(unlist(caeno_target_tables), dim = c(nrow(caeno_target_tables[[1]]), ncol(caeno_target_tables[[1]]), 10 ) )

dimnames(caeno_target_array) <- list(elegans_targets,
                                     unique(elegans_FIMO_censored$motif_id),
                                     caeno_FIMO_sp)

saveRDS(caeno_target_array,
        "output/caeno_target_array.RDS")

caeno_target_array <- readRDS("output/caeno_target_array.RDS")

TF_orthology_probs <- lapply(dimnames(caeno_target_array)[[2]], function(thisTFname){

  temp_it <- caeno_target_array[, thisTFname, ]
  
  message(paste0("Doing ", thisTFname, " which is number ", which(dimnames(caeno_target_array)[[2]] == thisTFname), " of ", length(dimnames(caeno_target_array)[[2]])))
  
  # note here checked the ranking is as I like it. ie. for ties, the min is true (note because negative vector)
  # i.e. from 100, if number 2 and 3 tie, they both get rank 3 i.e. top 0.03%
  
  percentiles <- apply(temp_it, 2, function(Y){

    temprank <- base::rank(-Y[!is.na(Y)], ties.method = "max")

    Y[!is.na(Y)] <- temprank / length(Y[!is.na(Y)])
    
    Y
    
  })
  
  target_products <- apply(percentiles, 1, prod, na.rm = TRUE)
  
  # compare this product to a distribution 
  
  shuffprods <- sapply(1:1000, function(i){
    
    shuffled <- apply(temp_it[, 1:10], 2, function(X){
      
      X[!is.na(X)] <- sample(X[!is.na(X)], size = sum(!is.na(X)), replace = FALSE)
      X
      
    })
    
    shuff_percentiles <- apply(shuffled, 2, function(Y){
      
      temprank <- base::rank(-Y[!is.na(Y)], ties.method = "max")
      
      Y[!is.na(Y)] <- temprank / length(Y[!is.na(Y)])
      
      Y
      
    })
    
    apply(shuff_percentiles, 1, prod, na.rm = TRUE)
    
  })
  
  product_vs_shuffle <- sapply(1:length(target_products), function(i){
    
    target_products[i]
    shuffprods[i, ]
    
    ecdf_fun(shuffprods[i, ], target_products[i])
    
  })
  
  names(product_vs_shuffle) <- names(target_products)

  adjusted_prod_vs_shuff <- p.adjust(product_vs_shuffle, method = "BH")
  names(adjusted_prod_vs_shuff) <- names(product_vs_shuffle)
  
  return(adjusted_prod_vs_shuff)
  
})

names(TF_orthology_probs) <- dimnames(caeno_target_array)[[2]]
saveRDS(TF_orthology_probs, "output/TF_orthology_probs.rds")
TF_orthology_probs <- readRDS("output/TF_orthology_probs.rds")

make.TF.orth.GRN <- function(cutoff,
                             fdrcut,
                             TF_orthology_probs = TF_orthology_probs,
                             elegans_target_table = elegans_target_table,
                             suffix = NULL){
# cutoff = 1000
# fdrcut = 0
tempGRN <- do.call(rbind, lapply(names(TF_orthology_probs), function(thisTF){
  # thisTF <- names(TF_orthology_probs)[5]
  x <- TF_orthology_probs[[thisTF]]
  
  fdrcutx <- names(x[!x > fdrcut])
  
  elegans_x <- elegans_target_table[, thisTF]
  names(elegans_x) <- row.names(elegans_target_table)
  
  elegans_x_order <- elegans_x[order(elegans_x, decreasing = TRUE)]
  
  thiscutoff <- min(cutoff, sum(elegans_x_order != 0))
  
  fdrcut_names <- names(elegans_x_order[1:thiscutoff])[names(elegans_x_order[1:thiscutoff]) %in% fdrcutx]
  
  data.frame("source" = rep(thisTF, times = length(fdrcut_names)),
               "target" = fdrcut_names,
               "weight" = rep(1, times = length(fdrcut_names)))
  
}))

# let's exclude fewer than 10
tempGRN <- tempGRN[!tempGRN$source %in% names(table(tempGRN$source))[table(tempGRN$source) < 10], ]

write.table(tempGRN,
            file = paste0("output/GRNs/TF_orthprobs_cut", cutoff, "_fdr", fdrcut, suffix,".txt"),
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE)

}

cutoffs_vec <- c(500, 1000, 1500, 2000, 2500, 3000, 5000, 10000)
fdrvec <- c(0, 0.05, 0.1)

for(i in 1:length(cutoffs_vec)){
  
  for(j in 1:length(fdrvec)){
    
    message(paste0("Now doing cutoff ", cutoffs_vec[i], " with FDR ", fdrvec[j]))
    make.TF.orth.GRN(cutoff = cutoffs_vec[i],
                     fdrcut = fdrvec[j])
    
    
  }
  
}

# from here, compare the best cutoffs etc.

top2500_fdr0 <- read.table("output/GRNs/TF_orthprobs_cut2500_fdr0.txt",
                           header = TRUE)

top2500_fdr0_noorthcontrol <- do.call(rbind, lapply(unique(top2500_fdr0$source), function(thisTF){
  # thisTF <- unique(top2500_fdr0$source)[1]
  number_of_targets <- table(top2500_fdr0$source)[thisTF]
  
  toptargets <- row.names(elegans_target_table[order(elegans_target_table[, thisTF], decreasing = TRUE), ])[1:number_of_targets]
  
  data.frame("source" = rep(thisTF, times = number_of_targets),
             "target" = toptargets,
             "weight" = rep(1, times = number_of_targets))
  
}))

write.table(top2500_fdr0_noorthcontrol,
            file = paste0("output/GRNs/top2500_fdr0_noorthcontrol.txt"),
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE)

## FIMOD 2000 control for comparison

FIMOD2000 <- read.table("~/Cel_GRN_benchmark/output/FIMOD_2000_nosign_equalweights.txt",
                        header = TRUE)


FIMOD2000_censored <- FIMOD2000[FIMOD2000$source %in% unique(top2500_fdr0$source), ]

write.table(FIMOD2000_censored,
            "output/GRNs/FIMOD2000_orth2500fdr0_control.txt",
            sep = "\t",
            col.names = TRUE,
            row.names = FALSE)


#### let's try to exclude orthologues of operon genes

# can do once for all of the targets

down_op_excl_list <- lapply(dimnames(caeno_target_array)[[3]], function(thissp){

  # thissp <- dimnames(caeno_target_array)[[3]][1]
  thisslice <- caeno_target_array[, , thissp]
  
  thisslice[row.names(thisslice) %in% all_downstream_operon_genes, ] <- NA

  return(thisslice)
  
})

down_op_excl_array <- array(unlist(down_op_excl_list), dim = c(nrow(down_op_excl_list[[1]]), ncol(down_op_excl_list[[1]]), length(down_op_excl_list)))

dimnames(down_op_excl_array) <- list(dimnames(caeno_target_array)[[1]],
                                     dimnames(caeno_target_array)[[2]],
                                     dimnames(caeno_target_array)[[3]])

TF_orthology_opexcl_probs <- lapply(dimnames(down_op_excl_array)[[2]], function(thisTFname){
  # thisTFname <- dimnames(caeno_target_array)[[2]][1]
  temp_it <- down_op_excl_array[, thisTFname, ]
  
  message(paste0("Doing ", thisTFname, " which is number ", which(dimnames(down_op_excl_array)[[2]] == thisTFname), " of ", length(dimnames(down_op_excl_array)[[2]])))
  
  # note here checked the ranking is as I like it. ie. for ties, the min is true (note because negative vector)
  # i.e. from 100, if number 2 and 3 tie, they both get rank 3 i.e. top 0.03%
  
  percentiles <- apply(temp_it, 2, function(Y){
    
    temprank <- base::rank(-Y[!is.na(Y)], ties.method = "max")
    
    Y[!is.na(Y)] <- temprank / length(Y[!is.na(Y)])
    
    Y
    
  })
  
  target_products <- apply(percentiles, 1, prod, na.rm = TRUE)
  
  # compare this product to a distribution 
  
  shuffprods <- sapply(1:1000, function(i){
    
    shuffled <- apply(temp_it[, 1:10], 2, function(X){
      
      X[!is.na(X)] <- sample(X[!is.na(X)], size = sum(!is.na(X)), replace = FALSE)
      X
      
    })
    
    shuff_percentiles <- apply(shuffled, 2, function(Y){
      
      temprank <- base::rank(-Y[!is.na(Y)], ties.method = "max")
      
      Y[!is.na(Y)] <- temprank / length(Y[!is.na(Y)])
      
      Y
      
    })
    
    apply(shuff_percentiles, 1, prod, na.rm = TRUE)
    
  })
  
  product_vs_shuffle <- sapply(1:length(target_products), function(i){
    
    target_products[i]
    shuffprods[i, ]
    
    ecdf_fun(shuffprods[i, ], target_products[i])
    
  })
  
  names(product_vs_shuffle) <- names(target_products)
  
  adjusted_prod_vs_shuff <- p.adjust(product_vs_shuffle, method = "BH")
  names(adjusted_prod_vs_shuff) <- names(product_vs_shuffle)
  
  return(adjusted_prod_vs_shuff)
  
})

names(TF_orthology_opexcl_probs) <- dimnames(down_op_excl_array)[[2]]
saveRDS(TF_orthology_opexcl_probs, "output/TF_orthology_opexcl_probs.rds")
# TF_orthology_opexcl_probs <- readRDS("output/TF_orthology_opexcl_probs.rds")

# to try, use best cutoffs from no exclusion for first look, which was 2500 and FDR = 0
make.TF.orth.GRN(cutoff = 2500,
                 fdrcut = 0,
                 suffix = "_operonEXCL",
                 elegans_target_table = elegans_target_table[!row.names(elegans_target_table) %in% all_downstream_operon_genes, ],
                 TF_orthology_probs = TF_orthology_opexcl_probs)


# all together - dont use operons








checkmeans <- apply(caeno_target_array, c(1, 2), mean, na.rm = TRUE)
checkmeans[is.nan(checkmeans)] <- NA

# how many species do we have this interaction for?
numbermeans <- apply(caeno_target_array, c(1, 2), function(x){sum(!is.na(x))})

# censor out ones that are not present in at least 4 other species
censormeans <- numbermeans
censormeans[censormeans < 4] <- NA

# take away ones that are not there in elegans
censormeans[elegans_target_table == 0] <- NA

# how many species do we have this interaction for?
significant_no <- apply(caeno_target_array, c(1, 2), function(x){sum(!is.na(x) & x != 0)})

# what proportion, considering only potential interactions possible in 4 species?
significant_prop <- significant_no / censormeans

conserved <- apply(significant_prop, 2, function(x){sum(x > 0.6, na.rm = TRUE)})
conserved[order(conserved, decreasing = TRUE)]

# how many interactions do we have?
table(numbermeans)

# will put a cutoff at 5 interactions.

checkmeans_morethan5 <- checkmeans
checkmeans_morethan5[numbermeans < 5] <- NA

# let's also take away the ones that are not interactions in C. elegans.
checkmeans_morethan5[elegans_target_table == 0] <- NA

# take top 100 interactions for each TF

top100_nohomo_GRN_noweight <- do.call(rbind, lapply(colnames(checkmeans_morethan5), function(thisTF){

  x <- checkmeans_morethan5[, thisTF]
  
  if(!any(!is.na(x))){
    
    return(NULL)
    
  } else {
    
  no_pos_not_na <- sum(!is.na(x) & x != 0)
  
  select_no <- min(no_pos_not_na, 100)

  targets <- names(x[order(x, decreasing = TRUE)][1:select_no])
  
  data.frame("source" = rep(thisTF, times = select_no),
                         "target" = targets,
                         "weight" = rep(1, times = select_no))
  
  }
  
}))

top100_nohomo_GRN_orthweight <- do.call(rbind, lapply(colnames(checkmeans_morethan5), function(thisTF){
  
  x <- checkmeans_morethan5[, thisTF]
  
  if(!any(!is.na(x))){
    
    return(NULL)
    
  } else {
    
    no_pos_not_na <- sum(!is.na(x) & x != 0)
    
    select_no <- min(no_pos_not_na, 100)
    
  targets <- x[order(x, decreasing = TRUE)][1:select_no]
  
  # weight between 0.5 for minimum and 1 for maximum
  weight <- 0.5 + ((targets - min(targets)) / (max(targets) - min(targets)))/2
  
  data.frame("source" = rep(thisTF, times = select_no),
             "target" = names(targets),
             "weight" = weight)
  
  }
  
}))

top100_nohomo_GRN_elegweight <- do.call(rbind, lapply(colnames(checkmeans_morethan5), function(thisTF){

  x <- checkmeans_morethan5[, thisTF]
  
  if(!any(!is.na(x))){
    
    return(NULL)
    
  } else {
  
    no_pos_not_na <- sum(!is.na(x) & x != 0)
    
    select_no <- min(no_pos_not_na, 100)
    
  targets <- x[order(x, decreasing = TRUE)][1:select_no]
  
  eleg_target_scores <- elegans_target_table[names(targets), thisTF]
  
  # weight between 0.5 for minimum and 1 for maximum
  weight <- 0.5 + ((eleg_target_scores - min(eleg_target_scores)) / (max(eleg_target_scores) - min(eleg_target_scores)))/2
  
  data.frame("source" = rep(thisTF, times = select_no),
             "target" = names(targets),
             "weight" = weight)
  }
  
}))

top100_for_R <- top100_nohomo_GRN_elegweight
colnames(top100_for_R)[3] <- "mor"

decoupleR::check_corr(top100_for_R)

dir.create("output/GRNs")

write.table(top100_nohomo_GRN_noweight,
            file = "output/GRNs/top100_nohomo_GRN_noweight.txt",
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE)

write.table(top100_nohomo_GRN_orthweight,
            file = "output/GRNs/top100_nohomo_GRN_orthweight.txt",
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE)

write.table(top100_nohomo_GRN_elegweight,
            file = "output/GRNs/top100_nohomo_GRN_elegweight.txt",
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE)




elegans_mean_cor <- sapply(1:ncol(elegans_target_table), function(i){

if(any(!is.na(checkmeans_morethan5[, i]))){
  
  return(cor.test(elegans_target_table[, i], checkmeans_morethan5[, i], method = "spearman", na.rm = TRUE)$estimate)
  
} else {
  
  return (NULL)
  
}
  
})

# let's try adding the elegans scores to the orthology scores

top100_nohomo_GRN_elegandorth <- do.call(rbind, lapply(colnames(checkmeans_morethan5), function(thisTF){

  x <- checkmeans_morethan5[, thisTF]
  elegans_x <- elegans_target_table[, thisTF]
  names(elegans_x) <- row.names(elegans_target_table)

  elegandorth <- elegans_x + x/2

  if(!any(!is.na(elegandorth))){
    
    return(NULL)
    
  } else {
    
    no_pos_not_na <- sum(!is.na(elegandorth) & elegandorth != 0)
    
    select_no <- min(no_pos_not_na, 100)
    
    targets <- elegandorth[order(elegandorth, decreasing = TRUE)][1:select_no]
    
    data.frame("source" = rep(thisTF, times = select_no),
               "target" = names(targets),
               "weight" = rep(1, times = select_no))
  }
  
}))

write.table(top100_nohomo_GRN_elegandorth,
            file = "output/GRNs/top100_nohomo_GRN_elegandorth.txt",
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE)

top500_nohomo_GRN_elegandorth <- do.call(rbind, lapply(colnames(checkmeans_morethan5), function(thisTF){
  
  x <- checkmeans_morethan5[, thisTF]
  elegans_x <- elegans_target_table[, thisTF]
  names(elegans_x) <- row.names(elegans_target_table)
  
  elegandorth <- elegans_x + x/2
  
  if(!any(!is.na(elegandorth))){
    
    return(NULL)
    
  } else {
    
    no_pos_not_na <- sum(!is.na(elegandorth) & elegandorth != 0)
    
    select_no <- min(no_pos_not_na, 500)
    
    targets <- elegandorth[order(elegandorth, decreasing = TRUE)][1:select_no]
    
    data.frame("source" = rep(thisTF, times = select_no),
               "target" = names(targets),
               "weight" = rep(1, times = select_no))
  }
  
}))

write.table(top500_nohomo_GRN_elegandorth,
            file = "output/GRNs/top500_nohomo_GRN_elegandorth.txt",
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE)

top1000_nohomo_GRN_elegandorth <- do.call(rbind, lapply(colnames(checkmeans_morethan5), function(thisTF){
  
  x <- checkmeans_morethan5[, thisTF]
  elegans_x <- elegans_target_table[, thisTF]
  names(elegans_x) <- row.names(elegans_target_table)
  
  elegandorth <- elegans_x + x/2
  
  if(!any(!is.na(elegandorth))){
    
    return(NULL)
    
  } else {
    
    no_pos_not_na <- sum(!is.na(elegandorth) & elegandorth != 0)
    
    select_no <- min(no_pos_not_na, 1000)
    
    targets <- elegandorth[order(elegandorth, decreasing = TRUE)][1:select_no]
    
    data.frame("source" = rep(thisTF, times = select_no),
               "target" = names(targets),
               "weight" = rep(1, times = select_no))
  }
  
}))

write.table(top1000_nohomo_GRN_elegandorth,
            file = "output/GRNs/top1000_nohomo_GRN_elegandorth.txt",
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE)

top500_nohomo_GRN_elegandorth_equal <- do.call(rbind, lapply(colnames(checkmeans_morethan5), function(thisTF){
  
  x <- checkmeans_morethan5[, thisTF]
  elegans_x <- elegans_target_table[, thisTF]
  names(elegans_x) <- row.names(elegans_target_table)
  
  elegandorth <- elegans_x + x
  
  if(!any(!is.na(elegandorth))){
    
    return(NULL)
    
  } else {
    
    no_pos_not_na <- sum(!is.na(elegandorth) & elegandorth != 0)
    
    select_no <- min(no_pos_not_na, 500)
    
    targets <- elegandorth[order(elegandorth, decreasing = TRUE)][1:select_no]
    
    data.frame("source" = rep(thisTF, times = select_no),
               "target" = names(targets),
               "weight" = rep(1, times = select_no))
  }
  
}))

write.table(top500_nohomo_GRN_elegandorth_equal,
            file = "output/GRNs/top500_nohomo_GRN_elegandorth_equal.txt",
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE)

top500_nohomo_GRN_elegandorth_equal_plusweight <- do.call(rbind, lapply(colnames(checkmeans_morethan5), function(thisTF){

  x <- checkmeans_morethan5[, thisTF]
  elegans_x <- elegans_target_table[, thisTF]
  names(elegans_x) <- row.names(elegans_target_table)
  
  elegandorth <- elegans_x + x
  
  if(!any(!is.na(elegandorth))){
    
    return(NULL)
    
  } else {
    
    no_pos_not_na <- sum(!is.na(elegandorth) & elegandorth != 0)
    
    select_no <- min(no_pos_not_na, 500)
    
    targets <- elegandorth[order(elegandorth, decreasing = TRUE)][1:select_no]
    
    weight <- 0.5 + ((targets - min(targets)) / (max(targets) - min(targets)))/2
    
    data.frame("source" = rep(thisTF, times = select_no),
               "target" = names(targets),
               "weight" = weight)
  }
  
}))

write.table(top500_nohomo_GRN_elegandorth_equal_plusweight,
            file = "output/GRNs/top500_nohomo_GRN_elegandorth_equal_plusweight.txt",
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE)

top500_nohomo_GRN_elegansonly <- do.call(rbind, lapply(colnames(checkmeans_morethan5), function(thisTF){

  elegans_x <- elegans_target_table[, thisTF]
  names(elegans_x) <- row.names(elegans_target_table)
  
  if(!any(!is.na(elegans_x))){
    
    return(NULL)
    
  } else {
    
    no_pos_not_na <- sum(!is.na(elegans_x) & elegans_x != 0)
    
    select_no <- min(no_pos_not_na, 500)
    
    targets <- elegans_x[order(elegans_x, decreasing = TRUE)][1:select_no]
  
    data.frame("source" = rep(thisTF, times = select_no),
               "target" = names(targets),
               "weight" = rep(1, times = select_no))
  }
  
}))

# throw ot ones that are not in the other sets
top500_nohomo_GRN_elegansonly <- top500_nohomo_GRN_elegansonly[top500_nohomo_GRN_elegansonly$source %in% top500_nohomo_GRN_elegandorth_equal$source, ]

colnames(top500_nohomo_GRN_elegansonly)[3] <- "mor"
elegansonly_corrcheck <- decoupleR::check_corr(top500_nohomo_GRN_elegansonly)

top500_nohomo_GRN_elegansonly <- top500_nohomo_GRN_elegansonly[!top500_nohomo_GRN_elegansonly$source %in% unlist(elegansonly_corrcheck[elegansonly_corrcheck$correlation > 0.99, c("source", "source.2")]), ]
colnames(top500_nohomo_GRN_elegansonly)[3] <- "weight"

write.table(top500_nohomo_GRN_elegansonly,
            file = "output/GRNs/top500_nohomo_GRN_elegansonly.txt",
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE)

# try one here where we have a yes/no; use the top elegans one but filter out any that are not conserved in at least 4 other species

top500_nohomo_GRN_elegans_caenopresent <- do.call(rbind, lapply(colnames(checkmeans_morethan5), function(thisTF){
  
  x <- checkmeans_morethan5[, thisTF]
  elegans_x <- elegans_target_table[, thisTF]
  names(elegans_x) <- row.names(elegans_target_table)
  
  elegandorth <- elegans_x + x
  
  if(!any(!is.na(elegans_x))){
    
    return(NULL)
    
  } else {
    
    no_pos_not_na <- sum(!is.na(elegans_x) & elegans_x != 0)
    
    select_no <- min(no_pos_not_na, 500)
    
    targets <- elegans_x[order(elegans_x, decreasing = TRUE)][1:select_no]
    
    data.frame("source" = rep(thisTF, times = select_no),
               "target" = names(targets),
               "weight" = rep(1, times = select_no))
  }
  
}))

write.table(top500_nohomo_GRN_elegansonly,
            file = "output/GRNs/top500_nohomo_GRN_elegansonly.txt",
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE)

#### look at distributions of overlaps for particular TF across species ####

all_caeno_list <- c(list(elegans = elegans_target_table), caeno_target_tables)
all_caeno_target_array <- array(unlist(all_caeno_list), dim = c(nrow(all_caeno_list[[1]]), ncol(all_caeno_list[[1]]), length(all_caeno_list)))

dimnames(all_caeno_target_array) <- list(elegans_targets,
                                     unique(elegans_FIMO_censored$motif_id),
                                     names(all_caeno_list))

# let's do this test on TFs only present in at least 5 species
distributions <- apply(all_caeno_target_array, 2, function(temp_it){
# temp_it <- all_caeno_target_array[, 5, ]

if(sum(apply(temp_it, 2, function(x){any(!is.na(x))})) < 5){
  return(NULL)
}

# these are the actual targets i.e. where a significant motif occurs in the promoter
temp_it_actualtargets <- apply(temp_it, 2, function(x){
  
  row.names(temp_it)[!is.na(x) & x > 0]
  
})

# these are the potential targets i.e. where the one2one orthology pair exists
targets_tempit <- apply(temp_it, 2, function(x){
  
  row.names(temp_it)[!is.na(x)]
  
})

number_of_targets <- apply(temp_it, 2, function(x){sum(!is.na(x) & x > 0)})

tempsamples <- lapply(1:length(targets_tempit), function(i){sample(targets_tempit[[i]], size = number_of_targets[[i]])})

wilcox.test(table(unlist(temp_it_actualtargets)),
            table(unlist(tempsamples)))$p.value
            
})

distributions_adj <- p.adjust(unlist(distributions), method = "BH")

plot(distributions_adj)
sum(distributions_adj < 0.1, na.rm = TRUE)

# 150/210 TFs present in at least 5 of the species have significant enrichment of overlapping targets.

# Here we calculate for each transcription factor a target cutoff of significant proportion beyond which 95% of the cutoffs 
signif_prop_distributions <- apply(all_caeno_target_array, 2, function(temp_it){
  
  ## This line for troubleshooting
  # temp_it <- all_caeno_target_array[, "W03C9.4", ]
  
  # if the transcription factor is present in fewer than 5 species, return NULL
  if(sum(apply(temp_it, 2, function(x){any(!is.na(x))})) < 5){
    return(NULL)
  }
  
  # these are the actual targets i.e. where a significant motif occurs in the promoter
  temp_it_actualtargets <- apply(temp_it, 2, function(x){
    
    row.names(temp_it)[!is.na(x) & x > 0]
    
  })
  
  # these are the potential targets i.e. where the one2one orthology pair exists
  targets_tempit <- apply(temp_it, 2, function(x){
    
    row.names(temp_it)[!is.na(x)]
    
  })
  
  number_of_targets <- apply(temp_it, 2, function(x){sum(!is.na(x) & x > 0)})
  
  iterate100distro <- lapply(1:100, function(arbitrary){
  
  tempsamples <- lapply(1:length(targets_tempit), function(i){

    sample(targets_tempit[[i]], size = number_of_targets[[i]])
    
    })
  
  names(tempsamples) <- names(targets_tempit)
  
  sample_array <- data.frame(matrix(nrow = nrow(temp_it),
                         ncol = ncol(temp_it)))
  
  row.names(sample_array) <- row.names(temp_it)
  colnames(sample_array) <- colnames(temp_it)
  
  sample_array[] <- FALSE
  
  sample_array[is.na(temp_it)] <- NA
  
  sample_array[] <- lapply(colnames(sample_array), function(thisspp){

    thisTFdata <- sample_array[, thisspp]
    thisTFdata[row.names(sample_array) %in% tempsamples[[thisspp]]] <- TRUE
    
    thisTFdata

  })
  
  sample_array_morethan5 <- sample_array[, 2:11]
  sample_array_morethan5[apply(sample_array_morethan5, 1, function(x){sum(!is.na(x))}) < 5, ] <- NA
  
  # let's also take away the ones that are not interactions in C. elegans.
  sample_array_morethan5[sample_array[, 1] == 0, ] <- NA

  prop_distro <- apply(sample_array_morethan5, 1, function(row){

    sum(row == TRUE, na.rm = TRUE) / sum(!is.na(row))
    
  })
    
  })
  
  alldistro <- unlist(iterate100distro)

  alldistro
  
})

signif_prop_distributions <- lapply(signif_prop_distributions, function(x){
  
  x[!is.nan(x)]
  
})

saveRDS(signif_prop_distributions,
        "~/Cel_GRN_orthology/signif_prop_distributions.rds")

signif_quantiles <- sapply(signif_prop_distributions, quantile, 0.95)

# also try just filtering worm ones on whats conserved

top500_nohomo_GRN_elegans_conserved <- do.call(rbind, lapply(colnames(checkmeans_morethan5), function(thisTF){

  elegans_x <- elegans_target_table[, thisTF]
  names(elegans_x) <- row.names(elegans_target_table)
  
  x <- significant_prop[, thisTF]
  
  significant_thisone <- x[x>0.6]
  significant_thisone <- significant_thisone[!is.na(significant_thisone)]
  
  if(!any(!is.na(elegans_x))){
    
    return(NULL)
    
  } else {
    
    no_pos_not_na <- sum(!is.na(elegans_x) & elegans_x != 0)
    
    select_no <- min(no_pos_not_na, 500)
    
    targets <- elegans_x[order(elegans_x, decreasing = TRUE)][1:select_no]
    
    conserved_targets <- targets[names(targets) %in% names(significant_thisone)]
    
    data.frame("source" = rep(thisTF, times = length(conserved_targets)),
               "target" = names(conserved_targets),
               "weight" = rep(1, times = length(conserved_targets)))
  }
  
}))

table(top500_nohomo_GRN_elegans_conserved$source)


top500_nohomo_GRN_elegans_conservedcontrol <- do.call(rbind, lapply(colnames(checkmeans_morethan5), function(thisTF){

  elegans_x <- elegans_target_table[, thisTF]
  names(elegans_x) <- row.names(elegans_target_table)
  
  x <- significant_prop[, thisTF]
  
  significant_thisone <- x[x>0.6]
  significant_thisone <- significant_thisone[!is.na(significant_thisone)]
  
  if(!any(!is.na(elegans_x))){
    
    return(NULL)
    
  } else {
    
    no_pos_not_na <- sum(!is.na(elegans_x) & elegans_x != 0)
    
    select_no <- min(no_pos_not_na, 500)
    
    targets <- elegans_x[order(elegans_x, decreasing = TRUE)][1:select_no]
    
    conserved_targets <- targets[names(targets) %in% names(significant_thisone)]
    
    data.frame("source" = rep(thisTF, times = length(conserved_targets)),
               "target" = names(targets)[0:length(conserved_targets)],
               "weight" = rep(1, times = length(conserved_targets)))
  }
  
}))

#### Here with empirically determined significant proportion cutoff ####

top500_nohomo_GRN_elegans_conservedempcutoff <- do.call(rbind, lapply(colnames(checkmeans_morethan5), function(thisTF){

  elegans_x <- elegans_target_table[, thisTF]
  names(elegans_x) <- row.names(elegans_target_table)
  
  x <- significant_prop[, thisTF]
  
  cutoff <- max(quantile(signif_prop_distributions[[thisTF]], 0.95), 0.6)
  
  significant_thisone <- x[x > cutoff]
  significant_thisone <- significant_thisone[!is.na(significant_thisone)]
  
  if(!any(!is.na(elegans_x))){
    
    return(NULL)
    
  } else {
    
    no_pos_not_na <- sum(!is.na(elegans_x) & elegans_x != 0)
    
    select_no <- min(no_pos_not_na, 500)
    
    targets <- elegans_x[order(elegans_x, decreasing = TRUE)][1:select_no]
    
    conserved_targets <- targets[names(targets) %in% names(significant_thisone)]
    
    data.frame("source" = rep(thisTF, times = length(conserved_targets)),
               "target" = names(conserved_targets),
               "weight" = rep(1, times = length(conserved_targets)))
  }
  
}))

table(top500_nohomo_GRN_elegans_conservedempcutoff$source)

top500_nohomo_GRN_elegans_conservedempcutoffcontrol <- do.call(rbind, lapply(colnames(checkmeans_morethan5), function(thisTF){
  
  elegans_x <- elegans_target_table[, thisTF]
  names(elegans_x) <- row.names(elegans_target_table)
  
  x <- significant_prop[, thisTF]
  
  cutoff <- max(quantile(signif_prop_distributions[[thisTF]], 0.95), 0.6)
  
  significant_thisone <- x[x > cutoff]
  significant_thisone <- significant_thisone[!is.na(significant_thisone)]
  
  if(!any(!is.na(elegans_x))){
    
    return(NULL)
    
  } else {
    
    no_pos_not_na <- sum(!is.na(elegans_x) & elegans_x != 0)
    
    select_no <- min(no_pos_not_na, 500)
    
    targets <- elegans_x[order(elegans_x, decreasing = TRUE)][1:select_no]
    
    conserved_targets <- targets[names(targets) %in% names(significant_thisone)]
    
    data.frame("source" = rep(thisTF, times = length(conserved_targets)),
               "target" = names(targets)[0:length(conserved_targets)],
               "weight" = rep(1, times = length(conserved_targets)))
  }
  
}))

write.table(top500_nohomo_GRN_elegans_conservedempcutoff,
            file = "output/GRNs/top500_nohomo_GRN_elegans_conservedempcutoff.txt",
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE)

write.table(top500_nohomo_GRN_elegans_conservedempcutoffcontrol,
            file = "output/GRNs/top500_nohomo_GRN_elegans_conservedempcutoffcontrol.txt",
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE)

# using an empirically determined cutoff doesn't work as well as just drawing the line at 0.6
# what about other arbitrary cutoffs then

# 0.5 much worse, 0.7 similaer but slightly worse, 0.8 too few targets
# I just hit the arbitrary sweet spot

#### CHANGE ARBITRARY CUTOFF ####
top500_nohomo_GRN_elegans_conservedpoint7 <- do.call(rbind, lapply(colnames(checkmeans_morethan5), function(thisTF){
  
  elegans_x <- elegans_target_table[, thisTF]
  names(elegans_x) <- row.names(elegans_target_table)
  
  x <- significant_prop[, thisTF]
  
  significant_thisone <- x[x>0.7]
  significant_thisone <- significant_thisone[!is.na(significant_thisone)]
  
  if(!any(!is.na(elegans_x))){
    
    return(NULL)
    
  } else {
    
    no_pos_not_na <- sum(!is.na(elegans_x) & elegans_x != 0)
    
    select_no <- min(no_pos_not_na, 500)
    
    targets <- elegans_x[order(elegans_x, decreasing = TRUE)][1:select_no]
    
    conserved_targets <- targets[names(targets) %in% names(significant_thisone)]
    
    data.frame("source" = rep(thisTF, times = length(conserved_targets)),
               "target" = names(conserved_targets),
               "weight" = rep(1, times = length(conserved_targets)))
  }
  
}))

table(top500_nohomo_GRN_elegans_conservedpoint6$source)


top500_nohomo_GRN_elegans_conservedcontrolpoint7 <- do.call(rbind, lapply(colnames(checkmeans_morethan5), function(thisTF){
  
  elegans_x <- elegans_target_table[, thisTF]
  names(elegans_x) <- row.names(elegans_target_table)
  
  x <- significant_prop[, thisTF]
  
  significant_thisone <- x[x>0.7]
  significant_thisone <- significant_thisone[!is.na(significant_thisone)]
  
  if(!any(!is.na(elegans_x))){
    
    return(NULL)
    
  } else {
    
    no_pos_not_na <- sum(!is.na(elegans_x) & elegans_x != 0)
    
    select_no <- min(no_pos_not_na, 500)
    
    targets <- elegans_x[order(elegans_x, decreasing = TRUE)][1:select_no]
    
    conserved_targets <- targets[names(targets) %in% names(significant_thisone)]
    
    data.frame("source" = rep(thisTF, times = length(conserved_targets)),
               "target" = names(targets)[0:length(conserved_targets)],
               "weight" = rep(1, times = length(conserved_targets)))
  }
  
}))

write.table(top500_nohomo_GRN_elegans_conservedpoint7,
            file = "output/GRNs/top500_nohomo_GRN_elegans_conserved_0.7.txt",
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE)

write.table(top500_nohomo_GRN_elegans_conservedcontrolpoint7,
            file = "output/GRNs/top500_nohomo_GRN_elegans_conservedcontrol_0.7.txt",
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE)

#### HERE CUTOFF 1000 ####

top1000_nohomo_GRN_elegans_conserved <- do.call(rbind, lapply(colnames(checkmeans_morethan5), function(thisTF){
  
  elegans_x <- elegans_target_table[, thisTF]
  names(elegans_x) <- row.names(elegans_target_table)
  
  x <- significant_prop[, thisTF]
  
  significant_thisone <- x[x>0.6]
  significant_thisone <- significant_thisone[!is.na(significant_thisone)]
  
  if(!any(!is.na(elegans_x))){
    
    return(NULL)
    
  } else {
    
    no_pos_not_na <- sum(!is.na(elegans_x) & elegans_x != 0)
    
    select_no <- min(no_pos_not_na, 1000)
    
    targets <- elegans_x[order(elegans_x, decreasing = TRUE)][1:select_no]
    
    conserved_targets <- targets[names(targets) %in% names(significant_thisone)]
    
    data.frame("source" = rep(thisTF, times = length(conserved_targets)),
               "target" = names(conserved_targets),
               "weight" = rep(1, times = length(conserved_targets)))
  }
  
}))

table(top1000_nohomo_GRN_elegans_conserved$source)


top1000_nohomo_GRN_elegans_conservedcontrol <- do.call(rbind, lapply(colnames(checkmeans_morethan5), function(thisTF){
  
  elegans_x <- elegans_target_table[, thisTF]
  names(elegans_x) <- row.names(elegans_target_table)
  
  x <- significant_prop[, thisTF]
  
  significant_thisone <- x[x>0.6]
  significant_thisone <- significant_thisone[!is.na(significant_thisone)]
  
  if(!any(!is.na(elegans_x))){
    
    return(NULL)
    
  } else {
    
    no_pos_not_na <- sum(!is.na(elegans_x) & elegans_x != 0)
    
    select_no <- min(no_pos_not_na, 1000)
    
    targets <- elegans_x[order(elegans_x, decreasing = TRUE)][1:select_no]
    
    conserved_targets <- targets[names(targets) %in% names(significant_thisone)]
    
    data.frame("source" = rep(thisTF, times = length(conserved_targets)),
               "target" = names(targets)[0:length(conserved_targets)],
               "weight" = rep(1, times = length(conserved_targets)))
  }
  
}))

top500_nohomo_GRN_elegans_conserved2 <- top500_nohomo_GRN_elegans_conservedcontrol <- do.call(rbind, lapply(colnames(checkmeans_morethan5), function(thisTF){

  elegans_x <- elegans_target_table[, thisTF]
  names(elegans_x) <- row.names(elegans_target_table)
  
  x <- significant_prop[, thisTF]
  
  significant_thisone <- x[x>0.6]
  significant_thisone <- significant_thisone[!is.na(significant_thisone)]
  
  if(!any(!is.na(elegans_x))){
    
    return(NULL)
    
  } else {
    
    # no_pos_not_na <- sum(!is.na(elegans_x) & elegans_x != 0)
    
    targets <- elegans_x[order(elegans_x, decreasing = TRUE)]
    
    conserved_targets <- targets[names(targets) %in% names(significant_thisone)]
    
    select_no <- min(length(conserved_targets), 500)
    
  data.frame("source" = rep(thisTF, times = length(conserved_targets[0:select_no])),
               "target" = names(conserved_targets)[0:select_no],
               "weight" = rep(1, times = select_no))
  }
  
}))

top500_nohomo_GRN_elegans_conservedcontrol2 <- do.call(rbind, lapply(colnames(checkmeans_morethan5), function(thisTF){
  
  elegans_x <- elegans_target_table[, thisTF]
  names(elegans_x) <- row.names(elegans_target_table)
  
  x <- significant_prop[, thisTF]
  
  significant_thisone <- x[x>0.6]
  significant_thisone <- significant_thisone[!is.na(significant_thisone)]
  
  if(!any(!is.na(elegans_x))){
    
    return(NULL)
    
  } else {
    
    # no_pos_not_na <- sum(!is.na(elegans_x) & elegans_x != 0)
    # 
    # select_no <- min(no_pos_not_na, 500)
    
    targets <- elegans_x[order(elegans_x, decreasing = TRUE)]
    
    conserved_targets <- targets[names(targets) %in% names(significant_thisone)]
    
    select_no <- min(length(conserved_targets), 500)
    
    data.frame("source" = rep(thisTF, times = length(conserved_targets[0:select_no])),
               "target" = names(targets)[0:select_no],
               "weight" = rep(1, times = select_no))

  }
  
}))

write.table(top500_nohomo_GRN_elegans_conserved,
            file = "output/GRNs/top500_nohomo_GRN_elegans_conserved.txt",
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE)

write.table(top500_nohomo_GRN_elegans_conservedcontrol,
            file = "output/GRNs/top500_nohomo_GRN_elegans_conservedcontrol.txt",
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE)

write.table(top1000_nohomo_GRN_elegans_conserved,
            file = "output/GRNs/top1000_nohomo_GRN_elegans_conserved.txt",
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE)

write.table(top1000_nohomo_GRN_elegans_conservedcontrol,
            file = "output/GRNs/top1000_nohomo_GRN_elegans_conservedcontrol.txt",
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE)

write.table(top500_nohomo_GRN_elegans_conserved2,
            file = "output/GRNs/top500_nohomo_GRN_elegans_conserved2.txt",
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE)

write.table(top500_nohomo_GRN_elegans_conservedcontrol2,
            file = "output/GRNs/top500_nohomo_GRN_elegans_conservedcontrol2.txt",
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE)

### get these GRNs right. conserved 1 and its control are to limit to 500 ans then compare. conserved2 is the other way around (compare then limit to 0500)

# so conserved1 works better than conserved2, by some margin, with fewer targets (mean of 58)
mean(table(top500_nohomo_GRN_elegans_conserved$source))
mean(table(top500_nohomo_GRN_elegans_conserved2$source))

library(openxlsx)
benchmark_SRA <- read.xlsx("~/Cel_GRN_benchmark/benchmark_SRA.xlsx",
                           sheet = 1)

benchmark_SRA_BM <- getBM(mart = parasite_mart,
                          values = benchmark_SRA$target_name,
                          filters = "wormbase_locusid",
                          attributes = c("wormbase_locus", "wormbase_gseq"))

# mean 67.5 targets for benchmarking set. benchmarking set is representative.
mean(table(top500_nohomo_GRN_elegans_conserved$source)[benchmark_SRA_BM[benchmark_SRA_BM$wormbase_gseq %in% names(table(top500_nohomo_GRN_elegans_conserved$source)), "wormbase_gseq"]])

## two things to try
# homotypic binding score for ordering elegans targets
# define conserved targets in a TF specific manner based on distribution
# this I do hybrid; 0.6, or if 95% percentile of distribution is higher, then that.

calc.score.by.motif <- function(motif_list = CisBP_full_separate_by_motif,
                                homotypic_binding = TRUE,
                                returns_parameter = 3) {
  
  
  outlist <- lapply(motif_list, function(x){

    sapply(unique(x$sequence_name), function(y){
      
      thisonedata <- x[x$sequence_name == y, ]
      
      n <- nrow(thisonedata)
      
      if(n == 1){
        
        return(thisonedata$score)
        
      } else {
        
        if(homotypic_binding == TRUE){
          
          thisonedata <- thisonedata[order(thisonedata$score, decreasing = TRUE), ]
          
          sumscore <- 0
          
          for(i in 1:n){
            
            thisscore <- thisonedata[i, 'score']
            sumscore <- sumscore + thisscore * (1 / returns_parameter^(i-1))
            
          }
          
          return(sumscore)
          
        } else {
          
          return(thisonedata$score[1])
          
        }
        
      }
      
    })
    
  })
  
  outlist_order <- lapply(outlist, function(x){
    
    x[order(x, decreasing = TRUE)]
    
  })
  
  return(outlist_order)
  
}

elegans_homotypic_score_3 <- calc.score.by.motif()

# with elegans homotypic score, replace all_caeno_array elegans scores with homotypic binding scores and repeat derivation of GRNs

elegans_homotypic_score_3 <- lapply(elegans_homotypic_score_3, function(x){
  
  tempout <- x[row.names(all_caeno_list[[1]])]
  tempout[is.na(tempout)] <- 0
  
  names(tempout) <- row.names(all_caeno_list[[1]])
  
  tempout
  
})

elegans_homotypic_score_3 <- elegans_homotypic_score_3[colnames(all_caeno_list[[1]])]

elegans_homotypic_score_3_df <- do.call(cbind, elegans_homotypic_score_3)

all_caeno_list_eleghomo <- all_caeno_list
all_caeno_list_eleghomo[["elegans"]] <- data.frame(elegans_homotypic_score_3_df)

top500_nohomo_GRN_elegans_conservedempcutoff_plushomo <- do.call(rbind, lapply(colnames(checkmeans_morethan5), function(thisTF){
  
  elegans_x <- elegans_homotypic_score_3_df[, thisTF]
  names(elegans_x) <- row.names(elegans_homotypic_score_3_df)
  
  x <- significant_prop[, thisTF]
  
  cutoff <- max(quantile(signif_prop_distributions[[thisTF]], 0.95), 0.6)
  
  significant_thisone <- x[x > cutoff]
  significant_thisone <- significant_thisone[!is.na(significant_thisone)]
  
  if(!any(!is.na(elegans_x))){
    
    return(NULL)
    
  } else {
    
    no_pos_not_na <- sum(!is.na(elegans_x) & elegans_x != 0)
    
    select_no <- min(no_pos_not_na, 500)
    
    targets <- elegans_x[order(elegans_x, decreasing = TRUE)][1:select_no]
    
    conserved_targets <- targets[names(targets) %in% names(significant_thisone)]
    
    data.frame("source" = rep(thisTF, times = length(conserved_targets)),
               "target" = names(conserved_targets),
               "weight" = rep(1, times = length(conserved_targets)))
  }
  
}))


top500_nohomo_GRN_elegans_conservedempcutoffcontrol_plushomo <- do.call(rbind, lapply(colnames(checkmeans_morethan5), function(thisTF){
  
  elegans_x <- elegans_homotypic_score_3_df[, thisTF]
  names(elegans_x) <- row.names(elegans_homotypic_score_3_df)
  
  x <- significant_prop[, thisTF]
  
  cutoff <- max(quantile(signif_prop_distributions[[thisTF]], 0.95), 0.6)
  
  significant_thisone <- x[x > cutoff]
  significant_thisone <- significant_thisone[!is.na(significant_thisone)]
  
  if(!any(!is.na(elegans_x))){
    
    return(NULL)
    
  } else {
    
    no_pos_not_na <- sum(!is.na(elegans_x) & elegans_x != 0)
    
    select_no <- min(no_pos_not_na, 500)
    
    targets <- elegans_x[order(elegans_x, decreasing = TRUE)][1:select_no]
    
    conserved_targets <- targets[names(targets) %in% names(significant_thisone)]
    
    data.frame("source" = rep(thisTF, times = length(conserved_targets)),
               "target" = names(targets)[0:length(conserved_targets)],
               "weight" = rep(1, times = length(conserved_targets)))
  }
  
}))

write.table(top500_nohomo_GRN_elegans_conservedempcutoff_plushomo,
            file = "output/GRNs/top500_nohomo_GRN_elegans_conservedempcutoff_plushomo.txt",
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE)

write.table(top500_nohomo_GRN_elegans_conservedempcutoffcontrol_plushomo,
            file = "output/GRNs/top500_nohomo_GRN_elegans_conservedempcutoffcontrol_plushomo.txt",
            sep = "\t",
            row.names = FALSE,
            col.names = TRUE)

# plus homo not helpful. 
  
  
  