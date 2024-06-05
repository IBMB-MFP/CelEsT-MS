
library(stringr)
library(openxlsx)

fullset_TFs_BM <- readRDS("output/fullset_TFs_BM.rds")

combo_GRN <- read.table("output/GRNs/allthree_equalweights.txt",
                        header = TRUE)

orthCelEsT <- read.table("output/GRNs/orthCelEsT_equalweights.txt",
                        header = TRUE)

maxCelEsT <- read.table("output/GRNs/orthCelEsTMAXcov_equalweights.txt",
                         header = TRUE)

TF_chip_info <- data.frame(t(readRDS("output/TF_chip_info.rds")))
TF_chip_info <- TF_chip_info[row.names(TF_chip_info) %in% c(fullset_TFs_BM$wormbase_locus, fullset_TFs_BM$wormbase_gseq),]

allTFs_labelused <- readRDS("output/allTFs_labelused.rds")

names(allTFs_labelused) <- allTFs_labelused
names(allTFs_labelused)[1:310] <- fullset_TFs_BM[match(allTFs_labelused[1:310], fullset_TFs_BM$wormbase_locus), "wormbase_gseq"]

CisBP_TFinfo_withmotif <- readRDS("output/CisBP_TFinfo_withmotif.rds")
CisBP_FIMO_censored <- CisBP_TFinfo_withmotif[CisBP_TFinfo_withmotif$wormbase_gseq %in% combo_GRN$source, ]

observations <- read.table("output/benchmark_observations.txt",
                           header = TRUE)

# remove ABI sequenced elt-7

observations <- observations[observations$target_gseq != "elt-7", ]

# add Walhout TFs with >=15 targets
walhoutEV1 <- openxlsx::read.xlsx("~/Downloads/msb167131-sup-0002-datasetev1.xlsx")
# limit to interactions described as high quality
walhoutEV1 <- walhoutEV1[walhoutEV1$`In.high-quality.dataset?` == 'yes', ]

### NEED TF_Chip_info to be better labelled!

# for Supplementary using final network

fullset_TFs_BM[, "in_ChIP?"] <- fullset_TFs_BM$wormbase_gseq %in% allmodERN_TFsONLY_BM$wormbase_gseq
fullset_TFs_BM[, "ChIP_devstages"] <- NA
fullset_TFs_BM[match(names(unlist(TF_chip_info[fullset_TFs_BM[fullset_TFs_BM$`in_ChIP?`, "wormbase_gseq"], "number_of_stages"])), fullset_TFs_BM$wormbase_gseq), "ChIP_devstages"] <- unlist(TF_chip_info[fullset_TFs_BM[fullset_TFs_BM$`in_ChIP?`, "wormbase_gseq"], "number_of_stages"])

fullset_TFs_BM[, "in_Motif?"] <- fullset_TFs_BM$wormbase_gseq %in% unique(CisBP_FIMO_censored$wormbase_gseq)
fullset_TFs_BM[fullset_TFs_BM$`in_Motif?`, "Motif_direct_or_indirect?"] <- CisBP_TFinfo_withmotif[match(fullset_TFs_BM[fullset_TFs_BM$`in_Motif?`, "wormbase_gseq"], CisBP_TFinfo_withmotif$wormbase_gseq), "TF_Status"]
fullset_TFs_BM[fullset_TFs_BM$`in_Motif?`, "MotifID"] <- CisBP_TFinfo_withmotif[match(fullset_TFs_BM[fullset_TFs_BM$`in_Motif?`, "wormbase_gseq"], CisBP_TFinfo_withmotif$wormbase_gseq), "Motif_ID"]

fullset_TFs_BM[, "in_eY1H?"] <- fullset_TFs_BM$wormbase_gseq %in% unique(walhoutEV1$Prey.sequence.Name)
fullset_TFs_BM[, "in_CelEsT?"] <- fullset_TFs_BM$wormbase_gseq %in% unique(combo_GRN$source)

fullset_TFs_BM[, "in_orthCelEsT?"] <- fullset_TFs_BM$wormbase_gseq %in% unique(orthCelEsT$source)
fullset_TFs_BM[, "in_maxCelEsT?"] <- fullset_TFs_BM$wormbase_gseq %in% unique(maxCelEsT$source)

fullset_TFs_BM[, "in_benchmarking_set?"] <- fullset_TFs_BM$wormbase_gseq %in% observations$target_gseq
fullset_TFs_BM[fullset_TFs_BM$`in_benchmarking_set?`, "No_benchmark_experiments"] <- table(observations$target_gseq)[fullset_TFs_BM[fullset_TFs_BM$`in_benchmarking_set?`, "wormbase_gseq"]]

#### NEEDS UNIPROT ANNOTATION 

write.xlsx(fullset_TFs_BM,
           "output/TFs_info_for_supptable.xlsx")
