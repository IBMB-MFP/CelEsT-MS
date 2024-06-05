#start with empty workspace

rm(list = ls(all = TRUE))

# turn off scientific notation for plots

options(scipen=10000)

#### set working directory ####

# here create new folder and set working directory within it

dir.create("~/Cel_GRN_manuscript")
setwd("~/Cel_GRN_manuscript")

# create subfolders for input, output and graphics

dir.create("input")

# into input folder, add input files 

dir.create("output")
dir.create("graphics")

dir.create("output/STREME_out")

#### LOAD PACKAGES & FUNCTIONS ####

## First specify the packages of interest

packages <- c("biomaRt",
              "dplyr",
              "stringr",
              "openxlsx",
              "readxl")

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

biocmanager_packages <- c("GenomicRanges", # for intersecting
                          "GenomicFeatures",
                          "TxDb.Celegans.UCSC.ce11.refGene", # for promoter sequences
                          "BSgenome.Celegans.UCSC.ce11",
                          "universalmotif",
                          "Biostrings"
) 

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

github_packages <- c("r-lib/conflicted",
                     "snystrom/memes") 

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

#### DEFINE FUNCTIONS ####

# KUDOS to Josh O'Brien for this function [https://stackoverflow.com/questions/12865218/getting-rid-of-asis-class-attribute]
unAsIs <- function(X){
  if("AsIs" %in% class(X)){
    class(X) <- class(X)[-match("AsIs", class(X))]
  }
  X
}

#### INPUT DATA ####

parasite_mart <- useMart("parasite_mart", dataset = "wbps_gene", host = "https://parasite.wormbase.org", port = 443)

parasites_BM_list <- readRDS("~/Cel_GRN_orthology/output/parasites_BM_list.rds")

STREME_notHOT50_topregions <- readRDS("output/STREME_notHOT50_topregions_ALLMODern.rds")
STREME_notHOT50_topregions_nocutoff <- readRDS("output/STREME_notHOT50nocutoff_topregions_ALLMODern.rds")

allCisBP <- read_meme("output/all_MEME_motifs.txt")
names(allCisBP) <- sapply(allCisBP, function(x){x["name"]})

CisBP_TFinfo_withmotif <- readRDS("output/CisBP_TFinfo_withmotif.rds")

CisBP_direct_TFs <- CisBP_TFinfo_withmotif[CisBP_TFinfo_withmotif$TF_Status == "D", "wormbase_gseq"]

## DEFINE OPERON LIST FOR EXCLUSION 

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

#### RUN TomTom ####

TomTom_notHOT50_pvals <- lapply(names(STREME_notHOT50_topregions[names(STREME_notHOT50_topregions) %in% CisBP_direct_TFs]), function(thisname){

  x <- STREME_notHOT50_topregions[[thisname]]
  
  message(paste0("Doing ", thisname, " which is ", which(names(STREME_notHOT50_topregions)[names(STREME_notHOT50_topregions) %in% CisBP_direct_TFs] == thisname), " of ", length(names(STREME_notHOT50_topregions)[names(STREME_notHOT50_topregions) %in% CisBP_direct_TFs])))
  
  tryCatch({
  if(!any(!is.na(x))){return(NULL)}
  }, error = function(e){
  })
  
  no_of_motifs <- length(unlist(x["motif"]))
  
  sapply(1:no_of_motifs, function(i){

    tempTOMTOM <- runTomTom(
      x[["motif"]][[i]],
      database = allCisBP,
      outdir = "auto",
      thresh = 10,
      min_overlap = 5,
      dist = "ed",
      evalue = TRUE,
      silent = TRUE,
      meme_path = "/opt/local/bin"
    )
    
    tempTOMTOM_df <- data.frame(tempTOMTOM[["tomtom"]])
    
    tempTOMTOM_df[tempTOMTOM_df$match_name == thisname, "match_pval"]
    
  })
  
})

names(TomTom_notHOT50_pvals) <- names(STREME_notHOT50_topregions[names(STREME_notHOT50_topregions) %in% CisBP_direct_TFs])

saveRDS(TomTom_notHOT50_pvals, "output/TomTom_notHOT50_pvals.rds")
# TomTom_notHOT50_pvals <- readRDS("output/TomTom_notHOT50_pvals.rds")

TomTom_notHOT50nocutoff_pvals <- lapply(names(STREME_notHOT50_topregions_nocutoff)[names(STREME_notHOT50_topregions_nocutoff) %in% CisBP_direct_TFs], function(thisname){

  x <- STREME_notHOT50_topregions_nocutoff[[thisname]]
  
  message(paste0("Doing ", thisname, " which is ", which(names(STREME_notHOT50_topregions_nocutoff)[names(STREME_notHOT50_topregions_nocutoff) %in% CisBP_direct_TFs] == thisname), " of ", length(names(STREME_notHOT50_topregions_nocutoff)[names(STREME_notHOT50_topregions_nocutoff) %in% CisBP_direct_TFs])))
  
  tryCatch({
    if(!any(!is.na(x))){return(NULL)}
  }, error = function(e){
  })
  
  no_of_motifs <- length(unlist(x["motif"]))
  
  sapply(1:no_of_motifs, function(i){
    
    tempTOMTOM <- runTomTom(
      x[["motif"]][[i]],
      database = allCisBP,
      outdir = "auto",
      thresh = 10,
      min_overlap = 5,
      dist = "ed",
      evalue = TRUE,
      silent = TRUE,
      meme_path = "/opt/local/bin"
    )
    
    tempTOMTOM_df <- data.frame(tempTOMTOM[["tomtom"]])
    
    tempTOMTOM_df[tempTOMTOM_df$match_name == thisname, "match_pval"]
    
  })
  
})

names(TomTom_notHOT50nocutoff_pvals) <- names(STREME_notHOT50_topregions_nocutoff)[names(STREME_notHOT50_topregions_nocutoff) %in% CisBP_direct_TFs]

saveRDS(TomTom_notHOT50nocutoff_pvals, "output/TomTom_notHOT50nocutoff_pvals.rds")
# TomTom_notHOT50nocutoff_pvals <- readRDS("output/TomTom_notHOT50nocutoff_pvals.rds")

#### ASSESS PERFORMANCE ####

# ok so I am going to have a proper look at how well STREME did in these cases

TomTom_notHOt50_performance <- sapply(TomTom_notHOT50_pvals, function(x){
  
  # x <- TomTom_notHOT_pvals[[10]]
  if(is.null(x)){
    
    output_vec <- c()
    
    output_vec["number_of_motifs"] <- NA
    output_vec["best_position"] <- NA
    output_vec["number_of_hits"] <- NA
    output_vec["best_pval"] <- NA
    
    return(output_vec)
    
  }
  
  number_of_motifs <- length(x)
  position_of_correct <- which(sapply(x, function(y){length(y) == 1}))
  
  best_position_of_correct <- min(position_of_correct) 
  number_of_similar_hits <- length(position_of_correct)
  
  best_pval <- min(unlist(x[position_of_correct]))
  
  output_vec <- c()
  
  output_vec["number_of_motifs"] <- number_of_motifs
  output_vec["best_position"] <- best_position_of_correct
  output_vec["number_of_hits"] <- number_of_similar_hits
  output_vec["best_pval"] <- best_pval
  
  output_vec
  
})

TomTom_notHOt50nocutoff_performance <- sapply(TomTom_notHOT50nocutoff_pvals, function(x){
  
  # x <- TomTom_notHOT_pvals[[10]]
  if(is.null(x)){
    
    output_vec <- c()
    
    output_vec["number_of_motifs"] <- NA
    output_vec["best_position"] <- NA
    output_vec["number_of_hits"] <- NA
    output_vec["best_pval"] <- NA
    
    return(output_vec)
    
  }
  
  number_of_motifs <- length(x)
  position_of_correct <- which(sapply(x, function(y){length(y) == 1}))
  
  best_position_of_correct <- min(position_of_correct) 
  number_of_similar_hits <- length(position_of_correct)
  
  best_pval <- min(unlist(x[position_of_correct]))
  
  output_vec <- c()
  
  output_vec["number_of_motifs"] <- number_of_motifs
  output_vec["best_position"] <- best_position_of_correct
  output_vec["number_of_hits"] <- number_of_similar_hits
  output_vec["best_pval"] <- best_pval
  
  output_vec
  
})

Hot50_performance <- data.frame(t(TomTom_notHOt50_performance))
Hot50nocutoff_performance <- data.frame(t(TomTom_notHOt50nocutoff_performance))

row.names(Hot50_performance)[which(is.na(Hot50_performance$best_position))]
row.names(Hot50nocutoff_performance)[which(is.na(Hot50nocutoff_performance$best_position))]

# missing ones are missing from both sets. So we can exlude the NAs safely.

hot50_vec <- Hot50_performance$best_position
hot50_vec <- hot50_vec[!is.na(hot50_vec)]
hot50_vec[is.infinite(hot50_vec)] <- 0

hot50nocutoff_vec <- Hot50nocutoff_performance$best_position
hot50nocutoff_vec <- hot50nocutoff_vec[!is.na(hot50nocutoff_vec)]
hot50nocutoff_vec[is.infinite(hot50nocutoff_vec)] <- 0


mean(Hot50_performance$number_of_motifs, na.rm = T)
mean(Hot50nocutoff_performance$number_of_motifs, na.rm = T)

sum(!is.infinite(Hot50_performance$best_position) & !is.na(Hot50_performance$best_position))

sum(Hot50_performance$best_position == 1, na.rm = T)
sum(Hot50nocutoff_performance$best_position == 1, na.rm = T)

sum(Hot50nocutoff_performance$number_of_hits != 0, na.rm = T)
sum(Hot50_performance$number_of_hits != 0, na.rm = T)

# others are nearly always better positioned in cutoff set, bar a couple of exceptions.

#### RUN DE NOVO MOTIFS IN OTHER SPECIES ####

cspp_promoter_filelist <- list.files("output/orthologue_promoter_FASTA/")
cspp_promoter_filelist <- cspp_promoter_filelist[str_detect(cspp_promoter_filelist, "fasta")]

cspp_one2one_DNAstringset_list <- lapply(cspp_promoter_filelist, function(x){
  
  readDNAStringSet(paste0("output/orthologue_promoter_FASTA/", x))
  
})

names(cspp_one2one_DNAstringset_list) <- str_extract(cspp_promoter_filelist, "^[a-z]+")

# also must do elegans

elegans_promoters <- promoters(genes(TxDb.Celegans.UCSC.ce11.refGene))

entrez_togseq_BM <- getBM(attributes = c("entrezgene_id", "wormbase_gseq"),
                          filters = "entrezgene_id",
                          values = names(elegans_promoters),
                          mart = parasite_mart)

elegans_promoters_gseq <- elegans_promoters[as.character(entrez_togseq_BM$entrezgene_id)]
mcols(elegans_promoters_gseq)["wormbase_gseq"] <- entrez_togseq_BM[match(names(elegans_promoters_gseq), entrez_togseq_BM$entrezgene_id), "wormbase_gseq"]

# limit to genes with one2one orthologue to at least one gene in other worms
elegans_orthologues <- unique(unlist(lapply(parasites_BM_list, function(x){

  x[x[,6] == "ortholog_one2one", "wormbase_gseq"]
    
})))

elegans_orthologues_promoters <- elegans_promoters_gseq[elegans_promoters_gseq$wormbase_gseq %in% elegans_orthologues]

# exclude downstream operon genes

elegans_orthologues_promoters <- elegans_orthologues_promoters[!elegans_orthologues_promoters$wormbase_gseq %in% all_downstream_operon_genes]

elegans_orthologues_promoters <- trim(elegans_orthologues_promoters)

elegans_orthologues_sequences <- get_sequence(elegans_orthologues_promoters,
                                              BSgenome.Celegans.UCSC.ce11)

names(elegans_orthologues_sequences) <- elegans_orthologues_promoters$wormbase_gseq

cspp_one2one_DNAstringset_list[["elegans"]] <- elegans_orthologues_sequences

FIMO_allspp_modERN_denovomotifs <- lapply(STREME_notHOT50_topregions, function(thisTF){
  
  if(length(thisTF) != 1){
  
  temp_motif_results <- lapply(cspp_one2one_DNAstringset_list, function(thisstringset){

    runFimo(sequences = thisstringset,
            motifs = unAsIs(thisTF$motif),
            parse_genomic_coord = FALSE,
            meme_path = "/opt/local/bin")
    
  })
  
  names(temp_motif_results) <- names(cspp_one2one_DNAstringset_list)
  
  return(temp_motif_results)
  
  } else {
    
    return(NULL)
    
  }
  
})

saveRDS(FIMO_allspp_modERN_denovomotifs, 
        "output/FIMO_allspp_modERN_denovomotifs.rds")

# FIMO_allspp_modERN_denovomotifs <- readRDS("output/FIMO_allspp_modERN_denovomotifs.rds")

#### BUILD TARGET ARRAY ####

FIMO_allspp_modERN_denovomotifs_arrays <- lapply(names(FIMO_allspp_modERN_denovomotifs), function(thisTF){
# thisTF <- names(FIMO_allspp_modERN_denovomotifs)[[1]]
  message(paste0("Now doing ", thisTF, ", which is number ", which(names(FIMO_allspp_modERN_denovomotifs) == thisTF), " of ", length(names(FIMO_allspp_modERN_denovomotifs))))
  
caeno_FIMO_list <- FIMO_allspp_modERN_denovomotifs[[thisTF]]

if(is.null(caeno_FIMO_list)) {
  
  return(NULL)
  
} else {

caeno_sp_by_motif <- lapply(caeno_FIMO_list, function(thissp_FIMO){
# thissp_FIMO <- caeno_FIMO_list[[1]]
  temp_by_motif <- lapply(unique(thissp_FIMO$motif_id), function(x){
    
    as.data.frame(thissp_FIMO[thissp_FIMO$motif_id == x,  ])
    
  })
  
  names(temp_by_motif) <- unique(thissp_FIMO$motif_id)
  
  temp_by_motif
  
})

names(caeno_sp_by_motif) <- names(caeno_FIMO_list)

caeno_calc_score_list <- lapply(names(caeno_sp_by_motif), function(this_sppname){
  # this_sppname <- names(caeno_sp_by_motif)[1]
  message(paste0("Spp: ", this_sppname))

  thisFIMO_list <- caeno_sp_by_motif[[this_sppname]]
  
  temp_calc_score_by_motif <- lapply(thisFIMO_list, function(x){
# this step should last ~ 2min for 5 motifs for C07G2.2
# b4 <- Sys.time()
    # x <- thisFIMO_list[[1]]
    tempout <- sapply(unique(x$seqnames), function(y){

      max(x[x$seqnames == y, "score"])
      
    })
    
    names(tempout) <- unique(x$seqnames)
    tempout
    # print(Sys.time() - b4)
  })
  
  names(temp_calc_score_by_motif) <- names(thisFIMO_list)
  
  temp_calc_score_by_motif_order <- lapply(temp_calc_score_by_motif, function(x){

    x[order(x, decreasing = TRUE)]
    
  })
  
  names(temp_calc_score_by_motif_order) <- names(temp_calc_score_by_motif)
  
  temp_calc_score_by_motif_order
  
})

names(caeno_calc_score_list) <- names(caeno_sp_by_motif)

elegans_targets <- unique(unlist(sapply(caeno_calc_score_list$elegans, names)))
elegans_targets <- elegans_targets[!elegans_targets %in% all_downstream_operon_genes]

elegans_target_table <- as.data.frame(lapply(caeno_calc_score_list$elegans, function(thismotif_hits){
  # thismotif_hits <- caeno_calc_score_list$elegans[[1]]
  tempvec <- vector(length = length(elegans_targets))
  names(tempvec) <- elegans_targets
  tempvec[] <- 0
  
  tempvec[match(names(thismotif_hits), names(tempvec))] <- thismotif_hits
  
  tempvec
  
}))

caeno_target_tables <- lapply(names(caeno_calc_score_list), function(thisspecies){
  # thisspecies <- names(caeno_calc_score_list)[[1]]
  thisspecies_scores <- caeno_calc_score_list[[thisspecies]]
  
  if(thisspecies != "elegans"){
  
  thissp_one2one_object <- readRDS(paste0("output/C", thisspecies, "_one2one_promoterseq.rds"))
  
  thisspecies_scores_ortholog <- as.data.frame(lapply(thisspecies_scores,  function(thisspthisTF){
    # thisspthisTF <- thisspecies_scores[[1]]
    tempvec <- vector(length = length(elegans_targets))
    names(tempvec) <- elegans_targets
    tempvec[] <- 0
    
    thissp_thisTF_gseqnames <- thissp_one2one_object[match(names(thisspthisTF), thissp_one2one_object[, 5]), "wormbase_gseq"]
    
    thissp_thisTF_scores <- thisspthisTF[thissp_thisTF_gseqnames %in% elegans_targets]
    names(thissp_thisTF_scores) <- thissp_thisTF_gseqnames[thissp_thisTF_gseqnames %in% elegans_targets]
    
    tempvec[names(thissp_thisTF_scores)] <- thissp_thisTF_scores
    
    # make NA for ones that are not in the one2one ortholog BM and so were not considered for inclusion in the study
    tempvec[!names(tempvec) %in% thissp_one2one_object$wormbase_gseq] <- NA
    
    # and if this TF is not present in this species, NA it all
    if(!thisTF %in% thissp_one2one_object$wormbase_gseq){
     
      tempvec[] <- NA
       
    }
    
    tempvec
    
  }))
  
  } else {
    
    thisspecies_scores_ortholog <- as.data.frame(lapply(thisspecies_scores,  function(thisspthisTF){
      # thisspthisTF <- thisspecies_scores[[1]]
      tempvec <- vector(length = length(elegans_targets))
      names(tempvec) <- elegans_targets
      tempvec[] <- 0

      tempvec[names(thisspthisTF)] <- thisspthisTF
      
      tempvec
      
    }))
    
  }

  thisspecies_scores_ortholog
  
})

names(caeno_target_tables) <- names(caeno_calc_score_list)

caeno_target_array <- array(unlist(caeno_target_tables), dim = c(nrow(caeno_target_tables[[1]]), ncol(caeno_target_tables[[1]]), 11 ))

dimnames(caeno_target_array) <- list(elegans_targets,
                                     names(caeno_calc_score_list$elegans),
                                     names(caeno_calc_score_list))

return(caeno_target_array)

} # end of else

})

names(FIMO_allspp_modERN_denovomotifs_arrays) <- names(FIMO_allspp_modERN_denovomotifs)

saveRDS(FIMO_allspp_modERN_denovomotifs_arrays,
        "output/FIMO_allspp_modERN_denovomotifs_arrays.rds")
