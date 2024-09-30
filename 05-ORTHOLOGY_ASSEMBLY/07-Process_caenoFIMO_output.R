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

setwd("~/Cel_GRN_manuscript")

#### DEFINE FUNCTIONS ####

ecdf_fun <- function(x,perc) ecdf(x)(perc)

make.TF.orth.GRN <- function(cutoff,
                             fdrcut,
                             TF_orthology_prob_list = TF_orthology_probs,
                             elegans_targets = elegans_target_table,
                             suffix = NULL,
                             mintargets = 15){
  # cutoff = 1500
  # fdrcut = 0.5
  tempGRN <- do.call(rbind, lapply(names(TF_orthology_probs), function(thisTF){
    # thisTF <- "C04G2.7"
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
  names(TF_orthology_probs)[!names(TF_orthology_probs) %in% unique(tempGRN$source)]

  tempGRN <- tempGRN[!tempGRN$source %in% names(table(tempGRN$source))[table(tempGRN$source) < mintargets], ]
  
  write.table(tempGRN,
              file = paste0("output/GRNs/TF_orthprobs_cut", cutoff, "_fdr", fdrcut, suffix,".txt"),
              sep = "\t",
              row.names = FALSE,
              col.names = TRUE)
  
}

#### LOAD PACKAGES & FUNCTIONS ####

## First specify the packages of interest

packages <- c("stringr",
              "biomaRt",
              "readxl", # for reading xls files (NOTE not xlsx)
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

parasites_BM_list <- readRDS("output/parasites_BM_list.rds")

fullset_TFs_BM <- readRDS("output/fullset_TFs_BM.rds")

CisBP_TFinfo_withmotif <- readRDS("output/CisBP_TFinfo_withmotif.rds")

# load up output from FIMO on command line. 
# This output file has been placed in the output folder prior to running this script

elegans_FIMO_full <- read.table("output/fimo_out/fimo.tsv",
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
blumenthal_s3_sh2exactmatch <- read_excel("input/Supplemental_File3.xls",
                                          sheet = 2,
                                          col_names = FALSE)

blumenthal_s3_sh3partmatch <- read_excel("input/Supplemental_File3.xls",
                                         sheet = 3,
                                         col_names = TRUE)

blumenthal_s3_sh5notmatch <- read_excel("input/Supplemental_File3.xls",
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
briggsae_confirmed_operons <- read.xlsx("input/Jhaveri_et_al_2022_briggsaeoperons_TableS7.xlsx",
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

saveRDS(caeno_calc_score_list, "output/caeno_calc_score_list.rds")
# caeno_calc_score_list <- readRDS("~/Cel_GRN_manuscript/output/caeno_calc_score_list.rds")

#### BUILD ARRAY OF TARGETS ####

elegans_targets <- unique(elegans_FIMO_censored$sequence_name)[!unique(elegans_FIMO_censored$sequence_name) %in% all_downstream_operon_genes]

elegans_target_table <- as.data.frame(lapply(calc_score_by_motif, function(thismotif_hits){

  thismotif_hits <- thismotif_hits[names(thismotif_hits) %in% elegans_targets]
  tempvec <- vector(length = length(elegans_targets))
  names(tempvec) <- elegans_targets
  tempvec[] <- 0
  
  tempvec[match(names(thismotif_hits), names(tempvec))] <- thismotif_hits
  
  tempvec
  
}))

caeno_target_tables <- lapply(names(caeno_calc_score_list), function(thisspecies){

  message(paste0("Doing ", thisspecies))
  
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

# exclude targets that are downstream in elegans or briggsae operons

saveRDS(caeno_target_array,
        "output/caeno_target_array.RDS")

# caeno_target_array <- readRDS("output/caeno_target_array.RDS")

# if I do 10k shuffles it would take about 42h

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
  
  shuffprods <- sapply(1:10000, function(i){
    
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
    
    ecdf_fun(shuffprods[i, ], target_products[i])
    
  })
  
  names(product_vs_shuffle) <- names(target_products)

  adjusted_prod_vs_shuff <- p.adjust(product_vs_shuffle, method = "BH")
  names(adjusted_prod_vs_shuff) <- names(product_vs_shuffle)
  
  return(adjusted_prod_vs_shuff)
  
})

names(TF_orthology_probs) <- dimnames(caeno_target_array)[[2]]
saveRDS(TF_orthology_probs, "output/TF_orthology_probs10k.rds")
# TF_orthology_probs <- readRDS("~/Cel_GRN_manuscript/output/TF_orthology_probs10k.rds")

cutoffs_vec <- c(500, 1000, 1500, 2000, 2500, 3000, 5000, 10000)
fdrvec <- c(0, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)

for(i in 1:length(cutoffs_vec)){
  
  for(j in 1:length(fdrvec)){
    
    message(paste0("Now doing cutoff ", cutoffs_vec[i], " with FDR ", fdrvec[j]))
    make.TF.orth.GRN(cutoff = cutoffs_vec[i],
                     fdrcut = fdrvec[j],
                     TF_orthology_prob_list = TF_orthology_probs,
                     elegans_targets = elegans_target_table)
    
    
  }
  
}

#### look at distributions of overlaps for particular TF across species ####

all_caeno_list <- c(list(elegans = elegans_target_table), caeno_target_tables)
all_caeno_target_array <- array(unlist(all_caeno_list), dim = c(nrow(all_caeno_list[[1]]), ncol(all_caeno_list[[1]]), length(all_caeno_list)))

dimnames(all_caeno_target_array) <- list(elegans_targets,
                                         unique(elegans_FIMO_censored$motif_id),
                                         names(all_caeno_list))

# let's do this test on TFs only present in at least 5 species
distributions <- apply(all_caeno_target_array, 2, function(temp_it){
  # temp_it <- all_caeno_target_array[, 209, ]

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
              table(unlist(tempsamples)),
              alternative = "greater")$p.value

})

distributions_adj <- p.adjust(unlist(distributions), method = "BH")

plot(unlist(distributions_adj))
sum(unlist(distributions_adj) < 0.2, na.rm = TRUE)

# 153/210 TFs present in at least 5 of the species have significant enrichment of overlapping targets.

#### make unfiltered 1500/fdr0.8 net for combination ####

make.TF.orth.GRN(cutoff = 1500,
                 fdrcut = 0.8,
                 TF_orthology_prob_list = TF_orthology_probs,
                 elegans_targets = elegans_target_table,
                 mintargets = 0,
                 suffix = "unfiltered")

