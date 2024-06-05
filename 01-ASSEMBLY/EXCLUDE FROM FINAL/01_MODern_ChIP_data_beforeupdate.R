#start with empty workspace

rm(list = ls(all = TRUE))

# turn off scientific notation for plots

options(scipen=10000)

#### set working directory ####

# here create new folder and set working directory within it

setwd("~/Cel_GRN_assembly")

dir.create("~/Cel_GRN_manuscript/")
setwd("~/Cel_GRN_manuscript")

# create subfolders for input, output and graphics

dir.create("input")

# into input folder, add input files 

dir.create("output")
dir.create("graphics")

#### DEFINE FUNCTIONS ####

return.manual.targets <- function(tf_name,
                                  peaks_GRanges = chip_peaks_overlap,
                                  cutoff = 500){
  #tf_name <- "ceh-60"
  print(tf_name)
  
  thisTF_peaks_GR <- peaks_GRanges[str_detect(peaks_GRanges$V4, paste0(tf_name, "_"))]
  
  # first the case in which there is only a single replicate
  if(length(unique(thisTF_peaks_GR$V4[str_detect(thisTF_peaks_GR$V4, paste0(tf_name, "_"))])) == 1){
    
    thisTF_peaks_df <- as.data.frame(thisTF_peaks_GR)
    
    thisTF_peaks_collapse <- thisTF_peaks_df %>% group_by(manual_gseq) %>% summarize(sum = sum(V7))
    
    if(nrow(thisTF_peaks_collapse) > cutoff){
      
      top_targets <- base::unname(unlist(thisTF_peaks_collapse[order(thisTF_peaks_collapse$sum, decreasing = TRUE), "manual_gseq"])[1:cutoff])
      
    } else {
      
      top_targets <- base::unname(unlist(thisTF_peaks_collapse[order(thisTF_peaks_collapse$sum, decreasing = TRUE), "manual_gseq"]))
      
    }
    
  }
  
  if(length(unique(thisTF_peaks_GR$V4[str_detect(thisTF_peaks_GR$V4, paste0(tf_name, "_"))])) > 1){
    
    thisTF_stages <- unique(thisTF_peaks_GR$V4[str_detect(thisTF_peaks_GR$V4, paste0(tf_name, "_"))])
    
    thisTF_peaks <- lapply(thisTF_stages, function(thisstage){
      
      thisTF_peaks_GR[str_detect(thisTF_peaks_GR$V4, thisstage), ]
      
    })
    
    allstagetargets <- unique(unlist(lapply(thisTF_peaks, function(zz){zz$manual_gseq})))
    
    number_of_appearances_df <- data.frame("appearances" = sapply(allstagetargets, function(thistarget){
      
      sum(sapply(thisTF_peaks, function(thispeakobject){
        
        thistarget %in% thispeakobject$manual_gseq
        
      }))
      
    }))
    
    # here we look at the targets in each stage, and add up the peaks targeting each target for each stage. 
    # then we return the mean of that sum for each stage overall to give an estimate of signal strength for the target.
    number_of_appearances_df[, "sum_signal"] <- sapply(row.names(number_of_appearances_df), function(thistarget){
      
      thisTF_peaks_subset <- thisTF_peaks[sapply(thisTF_peaks, function(wx){thistarget %in% wx$manual_gseq})]
      
      mean(sapply(thisTF_peaks_subset, function(wy){
        
        wy2 <- as.data.frame(wy[wy$manual_gseq == thistarget, ])
        
        sumtibble <-  wy2 %>% group_by(manual_gseq) %>% summarize(sum = sum(V7))
        
        return(sumtibble$sum)
        
      }))
      
    })
    
    # we order targets primarily by number of stages they are bound at, and then by the average strength of summed peaks across stages
    number_of_appearances_df[] <- number_of_appearances_df[with(number_of_appearances_df, order(appearances, sum_signal, decreasing = TRUE)), ]
    
    if(nrow(number_of_appearances_df) > cutoff){
      
      top_targets <- row.names(number_of_appearances_df)[1:cutoff]
      
    } else {
      
      top_targets <- row.names(number_of_appearances_df)
      
    }
    
  }
  
  return(top_targets)
  
}

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
                          "TxDb.Celegans.UCSC.ce11.refGene" # for promoter sequences
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

#### INPUT DATA ####

# MODern ChIP-seq peaks file downloaded from https://epic.gs.washington.edu/modERN/, tab = Worm / ChIPSeq by Gene / Download Aggregated Peaks
MODern_agg <- read.table("input/annotatedPeak.bed",
                         sep = "\t",
                         header = FALSE,
                         fill = TRUE)

# Extra ENCODE TFs data

ENCODE_filelist <- list.files("~/Cel_GRN_manuscript/input/ENCODE_extraTFs/")

ENCODE_datalist <- lapply(ENCODE_filelist, function(x){
  
  read.table(paste0("~/Cel_GRN_manuscript/input/ENCODE_extraTFs/", x),
             sep = '\t',
             header = FALSE,
             fill = TRUE)
  
})

names(ENCODE_datalist) <- ENCODE_filelist

# process to have gene and stage information in the 4th column as for modERN processed aggregated peaks

ENCODE_metadata <- read.xlsx("input/encode_extra_TFs.xlsx")

ENCODE_datalist <- lapply(ENCODE_filelist, function(thisone){

  thisdata <- ENCODE_datalist[[thisone]]
  thismetadata <- unlist(ENCODE_metadata[match(str_remove(thisone, "\\.bed\\.gz"), ENCODE_metadata$bednarrowpeakfile), ])
  
  thisdata[, 4] <- paste0(thismetadata["Target.gene.symbol"], "_", str_replace_all(unlist(thismetadata["Life.stage"]), pattern = " ", replacement = "_"))
  
  thisdata
  
})

ENCODE_aggregated <- do.call(rbind, ENCODE_datalist)

## Extra TF ChIPs from literature
## Pubmed IDs:
# 31346165
# 36897776
# 30956009
# 34499028
# 28874466
# 25773600

lit_chips_filelist <- list.files("input/additional_literature_chips/")

## REPLACE with supplementary for manuscript
lit_chips_metadata <- read.table("input/additional_literature_chips_filesforpeakcaller.txt",
                                 sep = "\t",
                                 header = TRUE)

lit_chips_datalist <- lapply(lit_chips_filelist, function(thisone){

  thisdata <- read.table(paste0("input/additional_literature_chips/", thisone),
             sep = '\t',
             header = FALSE,
             fill = TRUE)
  
  thisname <- lit_chips_metadata[match(str_remove(thisone, "\\.narrowPeak"), lit_chips_metadata$Accession), "target.TF"]
  
  thisdata[, 4] <- paste0(thisname, "_unknownstage")
  
  thisdata
  
})

lit_chips_aggregated <- do.call(rbind, lit_chips_datalist)

# must rename chromosomes to match the ENCODE files

lit_chips_aggregated[lit_chips_aggregated$V1 != "MtDNA", 1] <- paste0("chr", lit_chips_aggregated[lit_chips_aggregated$V1 != "MtDNA", 1])
lit_chips_aggregated[lit_chips_aggregated$V1 == "MtDNA", 1] <- "chrM"

# combine all

allChIP_agg <- rbind(MODern_agg[, 1:10],
                     ENCODE_aggregated,)

# here is previously processed RNA-seq data from CeNDR, MRN normalised and with transcript counts collapsed to gene level counts
# code for processing can be found in manuscript 'Histone methylation has a direct metabolic role in human cells'
CeNDR_normalised_counts_collapse <- readRDS("output/Cel_CeNDR_normalised_counts_gene_level.rds")

# filter for genes which are detected as expressed in at least 401 samples (of 609)
CeNDR_normalised_counts_collapse_mostlyexpressed <- CeNDR_normalised_counts_collapse[apply(CeNDR_normalised_counts_collapse, 1, function(x){sum(x > 0) > 400}), ] 

saveRDS(CeNDR_normalised_counts_collapse_mostlyexpressed,
        "output/CeNDR_normalised_counts_collapse_mostlyexpressed.rds")

parasite_mart <- useMart("parasite_mart", dataset = "wbps_gene", host = "https://parasite.wormbase.org", port = 443)

wTF3 <- read.xlsx("input/WTF3.xlsx")

# define operon downstream gene for exclusion - use Allen et al (Blumenthal lab) 2011 - 'A global analysis of C. elegans trans splicing' [doi: 10.1101/gr.113811.110]

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

#### NEW ATTEMPT TO DEFINE TFs ####

allTFs_nameused <- unique(str_extract(allChIP_agg$V4, "[a-z]{2,4}-[0-9]{1,3}"))
allTFs_nameuseddot <- unique(str_extract(allChIP_agg$V4, "[a-z]{2,4}-[0-9]{1,3}\\.[0-9]+"))
allTFs_nameuseddot <- allTFs_nameuseddot[!is.na(allTFs_nameuseddot)]

allTFs_nameused_comb <- c(allTFs_nameused, allTFs_nameuseddot)

# add B to lin-15B
allTFs_nameused_comb[allTFs_nameused_comb == "lin-15"] <- "lin-15B"

allTFs_gseqused <- unique(str_extract(allChIP_agg$V4, "[A-Z0-9]+\\.[0-9]{1,3}"))

# remove dot gene names from this vector
allTFs_gseqused <- allTFs_gseqused[sapply(allTFs_gseqused, nchar) != 3]
allTFs_gseqused <- allTFs_gseqused[!is.na(allTFs_gseqused)]

allTFs_nameused_BM <- getBM(attributes = c("wormbase_gseq",
                                           "wormbase_locus",
                                           "wormbase_gene"),
                            filters = "wormbase_locusid",
                            mart = parasite_mart,
                            values = allTFs_nameused_comb)

# missing ones are partial names from dot used in gene name (snpc-1, snpc-3, hmg-1). eliminate. 
allTFs_nameused_comb <- allTFs_nameused_comb[allTFs_nameused_comb %in% allTFs_nameused_BM$wormbase_locus]

allTFs_gseqused_BM <- getBM(attributes = c("wormbase_gseq",
                                           "wormbase_locus",
                                           "wormbase_gene"),
                            filters = "wormbase_gseqname",
                            mart = parasite_mart,
                            values = allTFs_gseqused)

# none missing here.
allTFs_gseqused %in% allTFs_gseqused_BM$wormbase_gseq

allmodERN_TFs_BM <- rbind(allTFs_nameused_BM, 
                          allTFs_gseqused_BM)

# make vector of the label used in all cases
allTFs_labelused <- c(allTFs_nameused_comb, allTFs_gseqused)

# are they all TFs? 

allmodERN_TFs_BM[!allmodERN_TFs_BM$wormbase_gseq %in% wTF3$Sequence.name, ]

allmodERN_TFsONLY_BM <- allmodERN_TFs_BM[allmodERN_TFs_BM$wormbase_gseq %in% wTF3$Sequence.name, ]


allTFs_labelused[!allTFs_labelused %in% newones]
newones[!newones %in% allTFs_labelused]
newones[str_detect(newones, "Y116")]
saveRDS(allmodERN_TFsONLY_BM, 
        "output/allmodERN_TFsONLY_BM.rds")

# 368 TFs in wTF3.0 and ENCODE/literature

#### FIND TARGETS FROM PEAKS FILE ####

# rather than trusting the somewhat mysterious target coluumn from modENCODE, I will define the targets myself such that I know what it represents
# will go for binding in promoter, defined as 1000bp upstream - 200bp downstream

# in order to this, I convert the ChIP peaks into a Genomic Ranges object to overlap with promoters

allChIP_agg_for_GR <- allChIP_agg
colnames(allChIP_agg)[1:3] <-  c("chromosome", "start", "end")

allChIP_agg_GR <- makeGRangesFromDataFrame(allChIP_agg, 
                                          keep.extra.columns = TRUE)

# get promoter sequences 

promoters <- promoters(genes(TxDb.Celegans.UCSC.ce11.refGene), 
                       upstream = 1000, 
                       downstream = 200)

# do some filtering based on the genes we are likely to be able to find. 
# our primary dataset is going to be the CeNDR dataset so we could restrict to genes present in this RNA-seq data
# by restricting to genes with expression in at least 400 samples (see #### INPUT DATA ####), we cut it down to 17000 genes,
# which is already 2.5x less than the total promoters just from the Txdb object
# with overlap between expressed genes and identified promoters, we cut it down to 16032

convert_names_to_entrez <- getBM(mart = parasite_mart,
                                 values = row.names(CeNDR_normalised_counts_collapse_mostlyexpressed),
                                 filters = "wormbase_gseqname",
                                 attributes = c("wormbase_gseq",
                                                "entrezgene_id",
                                                "wormbase_gene",
                                                "wormbase_locus"))

convert_names_to_entrez2 <-  convert_names_to_entrez
convert_names_to_entrez2[convert_names_to_entrez2$wormbase_locus == "", "wormbase_locus"] <- convert_names_to_entrez2[convert_names_to_entrez2$wormbase_locus == "", "wormbase_gseq"]

saveRDS(convert_names_to_entrez2, "output/convert_names_to_entrez2.rds")
# convert_names_to_entrez2 <- readRDS("output/convert_names_to_entrez2.rds")

promoters_censored <- promoters[promoters$gene_id %in% convert_names_to_entrez$entrezgene_id]

convert_operongenes_to_entrez <- getBM(mart = parasite_mart,
                                 values = downstream_operon_genes,
                                 filters = "wormbase_gseqname",
                                 attributes = c("wormbase_gseq",
                                                "entrezgene_id",
                                                "wormbase_gene",
                                                "wormbase_locus"))

promoters_censored_notdowninoperon <- promoters_censored[!promoters_censored$gene_id %in% convert_operongenes_to_entrez$entrezgene_id]

saveRDS(promoters_censored_notdowninoperon, "output/promoters_censored_notdowninoperon.rds")
# promoters_censored_notdowninoperon <- readRDS("output/promoters_censored_notdowninoperon.rds")

promoter_manual_overlaps <- findOverlapPairs(allChIP_agg_GR, promoters_censored_notdowninoperon)

chip_peaks_overlap <- S4Vectors::first(promoter_manual_overlaps)
mcols(chip_peaks_overlap)["manual_target"] <- S4Vectors::second(promoter_manual_overlaps)$gene_id
mcols(chip_peaks_overlap)["manual_genename"] <- convert_names_to_entrez2[match(unlist(mcols(chip_peaks_overlap)["manual_target"]), convert_names_to_entrez2$entrezgene_id), "wormbase_locus" ]
mcols(chip_peaks_overlap)["manual_gseq"] <- convert_names_to_entrez2[match(unlist(mcols(chip_peaks_overlap)["manual_target"]), convert_names_to_entrez2$entrezgene_id), "wormbase_gseq" ]

allTFS_manualtoptargets <- lapply(allTFs_labelused, return.manual.targets, cutoff = 10000)
names(allTFS_manualtoptargets) <- allTFs_labelused

saveRDS(allTFS_manualtoptargets,
        file = "output/MODernENCODE_manualtoptargets_operonexcluded.rds")

# allTFS_manualtoptargets <- readRDS("output/MODernENCODE_manualtoptargets_operonexcluded.rds")

## count stages

TF_chip_info <- sapply(allTFs_labelused, function(tf_name){

  thisTF_peaks_GR <- allChIP_agg_for_GR[str_detect(allChIP_agg_for_GR$V4, paste0(tf_name, "_")), ]
  
  number_of_stages = length(unique(thisTF_peaks_GR$V4[str_detect(thisTF_peaks_GR$V4, paste0(tf_name, "_"))]))
  
  name_of_stages = str_remove(unique(thisTF_peaks_GR$V4[str_detect(thisTF_peaks_GR$V4, paste0(tf_name, "_"))]), paste0(tf_name, "_"))
  
  number_of_targets <- sapply(name_of_stages, function(thisstage){
    
    length(thisTF_peaks_GR$V4[str_detect(thisTF_peaks_GR$V4, paste0(thisstage))])
    
  })
  
  temp_output_list <- list()
  
  temp_output_list[["number_of_stages"]] <- number_of_stages
  
  temp_output_list[["name_of_stages"]] <- name_of_stages
  
  temp_output_list[["number_of_targets"]] <- number_of_targets
  
  return(temp_output_list)
  
})

saveRDS(TF_chip_info, "output/TF_chip_info")

#### HOT target exclusion ####

## Will want to exclude HOT regions. This is Table S10 of Kudron et al.

download.file("https://figshare.com/ndownloader/files/10099716",
              destfile = "input/Kudron_et_al_2018_modERN_TableS10.txt")

# The file has different number of columns because it lists on each line all the experiments which show a peak in the HOT region in question
#  just want the first 4 lines; chrom, start, end, number of peaks, without worrying too much which ones they are
# solution code Kudos to https://stackoverflow.com/questions/68139359/how-to-read-first-four-columns-from-a-file-with-different-number-of-columns-on-e

# read the file as text lines
HOTtxt_lines <- readLines("input/Kudron_et_al_2018_modERN_TableS10.txt")
# split by one or more spaces
HOTtxt <- base::strsplit(HOTtxt_lines, " +")
# keep only the vector elements with more than 0 chars
HOTtxt <- lapply(HOTtxt, function(x) x[sapply(x, nchar) > 0])
# the last line may have a '\n' only, remove it
HOTtxt <- HOTtxt[lengths(HOTtxt) > 0]
# now extract the first 4 elements of each vector
HOTtxt <- lapply(HOTtxt, '[', 1:4)
# and rbind to data.frame
HOT_df <- do.call(rbind.data.frame, HOTtxt)
names(HOT_df) <- c("chrom", "start", "end", "peakno")

HOT_df$peakno <- as.numeric(HOT_df$peakno)

## Kudron define HOT regions by number of experiments but some of them are separate experiments with the same TF
# and I would rather the different TFs as replicate / diff stage for same TF binding in same place should not be any kind of worry

HOT_exp <- lapply(HOTtxt_lines, function(x){
  
  x_split <- base::strsplit(x, " +")
  
  if(unlist(x_split)[4] == 1){
    
    return(NULL)
    
  }
  
  x_temp <- unlist(x_split)[5:length(unlist(x_split))]
  
  ## Here we WRAP in str_remove for two freaky exceptions that otherwise dont make much sense; remove the fem--2 (presumably refers to him-8 background) and remove xtl1186_
  # x_TFs <- str_remove(str_remove(str_remove(x_temp, "ce_"), "[A-Z]{2,3}[0-9]+_"), "_[A-Za-z0-9]+_[A-Z]{2}:.*$")
  x_TFs <- str_remove(str_remove(str_remove(str_remove(str_remove(x_temp, "ce_"), "[A-Z]{2,3}[0-9]+_"), "fem-2_"), "xtl1186_"), "_[A-Za-z0-9]+_[A-Z]{2}:.*$")
  
  unique(x_TFs)
  
})

HOT_df[, "distinctpeak_no"] <- sapply(HOT_exp, length)

# needs mitochondrial DNA to be same as in the other one
HOT_df[HOT_df$chrom == "chrMtDNA", "chrom"] <- "chrM"

# for my purposes want to define HOT peaks as.... 
# with more than 50 separate TFs of the 217 in Kudron et al.
sum(HOT_df$distinctpeak_no > 50)

HOT_cluster_GR <- makeGRangesFromDataFrame(HOT_df[HOT_df$distinctpeak_no > 50, 1:3])

allChIP_agg_peaksummits_HOToverlaps <- findOverlaps(allChIP_agg_GR, HOT_cluster_GR)

# Remove peaks from agg file that are in HOT regions

unique(queryHits(allChIP_agg_peaksummits_HOToverlaps))

# 598325 peaks overlap HOT regions; around 51.5%

allChIP_agg_HOTexcl_GR <- allChIP_agg_GR[-queryHits(allChIP_agg_peaksummits_HOToverlaps)]

HOTexcl_promoter_manual_overlaps <- findOverlapPairs(allChIP_agg_HOTexcl_GR, promoters_censored_notdowninoperon)

HOTexcl_chip_peaks_overlap <- S4Vectors::first(HOTexcl_promoter_manual_overlaps)
mcols(HOTexcl_chip_peaks_overlap)["manual_target"] <- S4Vectors::second(HOTexcl_promoter_manual_overlaps)$gene_id
mcols(HOTexcl_chip_peaks_overlap)["manual_genename"] <- convert_names_to_entrez2[match(unlist(mcols(HOTexcl_chip_peaks_overlap)["manual_target"]), convert_names_to_entrez2$entrezgene_id), "wormbase_locus" ]
mcols(HOTexcl_chip_peaks_overlap)["manual_gseq"] <- convert_names_to_entrez2[match(unlist(mcols(HOTexcl_chip_peaks_overlap)["manual_target"]), convert_names_to_entrez2$entrezgene_id), "wormbase_gseq" ]

allTFS_manualtoptargets_HOTexcl <- lapply(allTFs_labelused, return.manual.targets, cutoff = 10000, peaks_GRanges = HOTexcl_chip_peaks_overlap)
names(allTFS_manualtoptargets_HOTexcl) <- allTFs_labelused

saveRDS(allTFS_manualtoptargets_HOTexcl,
        file = "output/MODernENCODE_manualtoptargets_operon&HOTexcluded.rds")

# HOT exclusion leads to 50 TFs with fewer than 20 targets
sum(sapply(allTFS_manualtoptargets_HOTexcl, length) < 20)
