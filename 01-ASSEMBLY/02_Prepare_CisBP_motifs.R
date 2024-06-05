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
dir.create("graphics")

# make a subfolder in output for depositing the motifs converted to uniprobe format for use with FIMO
dir.create("output/uniprobe_conversions")

#### DEFINE FUNCTIONS ####

# this is a function from the now-defunct RADami package
write.DNAStringSet <- function(x,
                               format= c('phylip', 'fasta'),
                               padding = 30,
                               filename = 'DNAStringSetOut.phy',
                               fastaPrefix = '>') {
  
    # writes a sequence matrix to phylip or fasta format
  
    x.width <- width(x)[1]
    
    x <- as.character(x)
    
    if(format[1] == 'phylip') {
      
      for(i in 1:length(x)) x[i] <- paste(names(x)[i], paste(rep(" ", (padding - nchar(names(x)[i]))), collapse = ''), x[i], sep = '')
      
      writeLines(c(paste(length(x), x.width), x), filename)
      
    }
    
    if(format[1] == 'fasta') {
      
      out <- as.character(matrix(c(paste(fastaPrefix, names(x), sep = ''), x), nrow = 2, byrow = T))
      
      writeLines(out, filename)
      
    }
    
    return(0)
    
  }

#### LOAD PACKAGES & FUNCTIONS ####

## First specify the packages of interest

packages <- c("biomaRt",
              "stringr")

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
                          "TxDb.Celegans.UCSC.ce11.refGene", # for promoter ranges
                          "BSgenome.Celegans.UCSC.ce11") # for promoter sequences

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

if(!require(memes)){
  
  remotes::install_github("snystrom/memes")
  library(memes)
  
}

#### INPUT DATA ####

# set biomaRt site to parasite_mart for worms
parasite_mart <- useMart("parasite_mart", dataset = "wbps_gene", host = "https://parasite.wormbase.org", port = 443)

# C. elegans species archive downloaded from CisBP TFBS database [http://cisbp.ccbr.utoronto.ca/bulk.php]
# incudes these information files

CisBP_TFinfo <- read.table("input/Caenorhabditis_elegans_2023_02_08_5_54_am/TF_Information.txt",
                           sep = "\t",
                           header = TRUE)

# select only ones with good evidence.
# don't want motifs inferred from ChIP-seq as I want this to be a line of evidence fully independent from the 
# MODern ChIP-seq approach
# so select motifs based on protein binding microarrays or SELEX
CisBP_TFinfo_withmotif <- CisBP_TFinfo[CisBP_TFinfo$Motif_Type %in% c("PBM", "SELEX") , ]

# one WormBase Gene ID is not listed correctly, for HLH27
# identifier can be copied from anther entry for the same gene

CisBP_TFinfo_withmotif[CisBP_TFinfo_withmotif$DBID == "HLH27", "DBID"] <- CisBP_TFinfo_withmotif[CisBP_TFinfo_withmotif$TF_Name == "hlh-27", "DBID"]
CisBP_TFinfo_withmotif[CisBP_TFinfo_withmotif$TF_Name == "HLH27", "TF_Name"] <- "hlh-27"

# we have 4 genes with duplications. 
# HLH-1 has mutant motifs; don't want that. We keep the WT motif only
# efl-2 has Weiracuh 2014 and Narasimhan 2015; let's take the one from the later publication (same lab)
# nhr-182 has two from the same publication! then I suppose it doesnt matter much
# hlh-27 has a direct and an inferred motif. Let's keep the direct one

CisBP_TFinfo_withmotif[duplicated(CisBP_TFinfo_withmotif$TF_Name)|duplicated(CisBP_TFinfo_withmotif$TF_Name, fromLast = TRUE), ]
CisBP_TFinfo_withmotif[duplicated(CisBP_TFinfo_withmotif$DBID)|duplicated(CisBP_TFinfo_withmotif$DBID, fromLast = TRUE), ]

CisBP_TFinfo_withmotif <- CisBP_TFinfo_withmotif[!str_detect(CisBP_TFinfo_withmotif$DBID.1, "HLH-1_"), ]
CisBP_TFinfo_withmotif <- CisBP_TFinfo_withmotif[!str_detect(CisBP_TFinfo_withmotif$DBID.1, "HLH1"), ]

CisBP_TFinfo_withmotif <- CisBP_TFinfo_withmotif[!(CisBP_TFinfo_withmotif$TF_Name == "efl-2" & CisBP_TFinfo_withmotif$MSource_Year == "2014"), ]
CisBP_TFinfo_withmotif <- CisBP_TFinfo_withmotif[!duplicated(CisBP_TFinfo_withmotif$TF_ID), ]

CisBP_TFinfo_withmotif <- CisBP_TFinfo_withmotif[!(CisBP_TFinfo_withmotif$TF_Name == "hlh-27" & CisBP_TFinfo_withmotif$TF_Status == "I"), ]

# want to give them gseq identifiers for easier integration across datasets
CisBP_motifs_BM <- getBM(attributes = c("wormbase_gseq",
                                                "wormbase_locus",
                                                "wormbase_gene"),
                                 filters = "wbps_gene_id",
                                 mart = parasite_mart,
                                 values = CisBP_TFinfo_withmotif$DBID)

CisBP_TFinfo_withmotif[, "wormbase_gseq"] <- CisBP_motifs_BM[match(CisBP_TFinfo_withmotif$DBID, CisBP_motifs_BM$wormbase_gene), "wormbase_gseq"]

# save the updated table
saveRDS(CisBP_TFinfo_withmotif,
        "output/CisBP_TFinfo_withmotif.rds")

CisBP_motifs <- lapply(CisBP_TFinfo_withmotif$Motif_ID, function(x){
  
  read.table(paste0("input/Caenorhabditis_elegans_2023_02_08_5_54_am/pwms_all_motifs/", x, ".txt"))
  
})

names(CisBP_motifs) <- CisBP_motifs_BM[match(CisBP_TFinfo_withmotif$DBID, CisBP_motifs_BM$wormbase_gene), "wormbase_gseq"]

#### CONVERT CisBP MOTIFS FOR USE WITH FIMO ####

# convert all of these motifs to uniprobe for use with FIMO on the command line

for(i in 1:length(CisBP_motifs)){

  thisonetranspose <- t(CisBP_motifs[[i]][, 2:ncol(CisBP_motifs[[i]])])
  thisonetranspose[, 1] <- paste0(thisonetranspose[, 1], ":")
  output_mat <- rbind(c(names(CisBP_motifs)[i], rep("", times = (ncol(thisonetranspose)-1))), thisonetranspose)
  
  write.table(output_mat,
              paste0("output/uniprobe_conversions/", names(CisBP_motifs)[i], ".txt"),
              sep = "\t",
              col.names = FALSE,
              row.names = FALSE,
              quote = FALSE
  )
  
}

#### PREPARE PROMOTER SEQUENCE FASTA FOR USE WITH FIMO ####

promoters_censored <- readRDS("output/promoters_censored_notdowninoperon.rds")

convert_names_to_entrez2 <- readRDS("output/convert_names_to_entrez2.rds")

Cel_promoters_censored_seq <- get_sequence(promoters_censored,
                                           BSgenome.Celegans.UCSC.ce11)

# rename it according to entregene_ids
names(Cel_promoters_censored_seq) <- convert_names_to_entrez2[match(promoters_censored$gene_id, convert_names_to_entrez2$entrezgene_id), "wormbase_gseq"]

write.DNAStringSet(Cel_promoters_censored_seq,
                   format = "fasta",
                   filename = "output/Cel_promoters_censored.fasta")

