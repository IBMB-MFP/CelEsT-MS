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

dir.create("~/Cel_GRN_manuscript")
setwd("~/Cel_GRN_manuscript")

# create subfolders for input, output and graphics

dir.create("input")
dir.create("input/sequences")

# into input folder, add input files 

dir.create("output")
dir.create("output/BSGenome_pkgs")

dir.create("graphics")

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

# this is my adaptation of readToList from chromseq, which seems to not work as it should

readToList.MFP <- function(id, text, con){
  
  pos <- match(id, text)
  
  tex = list()
  
  for (i in 1:(length(pos) - 1)) {
    
    message(paste0("Reading contig ", i, " of ", length(id)))
    
    tex[[i]] = readLines(con, n = pos[i + 1])[pos[i]:(pos[i + 1] - 1)]
    
    if (i == (length(pos) - 1)) {
      
      tempread <- readLines(con, n = -1)
      
      tex[[i + 1]] <- readLines(con, n = -1)[pos[i + 1]:length(tempread)]
      
    }
    
  }
  
  return(tex)
  
}

collapse_number_ranges <- function(numbervec){
  
  numbervec <- numbervec[order(numbervec)]
  
  diffs <- c(1, diff(numbervec))
  
  start_indexes <- c(1, which(diffs > 1))
  end_indexes <- c(start_indexes - 1, length(numbervec))
  
  coloned <- paste(numbervec[start_indexes], numbervec[end_indexes], sep=":")
  
  for(i in 1:length(coloned)) {
    
    coloned_split <- unlist(str_split(coloned[i], ":"))
    
    if(coloned_split[1] == coloned_split[2]){
      
      coloned[i] <- coloned_split[1]
      
    }
    
  }
  
  paste0(coloned, collapse=",")
  
}

#### LOAD PACKAGES & FUNCTIONS ####

## First specify the packages of interest

packages <- c("biomaRt",
              "stringr",
              "chromseq",
              "conflicted")

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
                          "BSgenome.Celegans.UCSC.ce11",
                          "BSgenome") # for promoter sequences

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

caeno_genomes <- read.table("input/caenorhabditis_genomes.txt",
                            sep = "\t",
                            header = TRUE)

subFas_list <- list()

for(i in 1:nrow(caeno_genomes)){

  tempFASTA_split <- unlist(str_split(caeno_genomes[i, "Seq_link"], "/"))
  tempFASTA_name <- tempFASTA_split[length(tempFASTA_split)]
  
  if(!file.exists(paste0("input/", tempFASTA_name))){
    
    download.file(url = caeno_genomes[i, "Seq_link"],
                  destfile = paste0("input/", tempFASTA_name))
    
  }
  
  # preliminary read in
  temp_sp_readL <- readLines(paste0("input/", tempFASTA_name))
  
  # identify subentries in FASTA file
  temp_sp_subfas <- subFasID(temp_sp_readL)
  
  subFas_list[[i]] <- temp_sp_subfas
  
}

names(subFas_list) <- caeno_genomes$specific_name

saveRDS(subFas_list, "output/subFas_list.rds")
# subFas_list <- readRDS("output/subFas_list.rds")

#### C brenneri ####

brenner_contig_numbers <- as.numeric(str_remove(subFas_list$brenneri, '^>Cbre_Contig'))[order(as.numeric(str_remove(subFas_list$brenneri, '^>Cbre_Contig')))]

# >Cbre_contig is unique prefix
subFas_list$brenneri[!str_detect(subFas_list$brenneri, ">Cbre_Contig")]

brenner_contigs <-  collapse_number_ranges(numbervec = brenner_contig_numbers[!is.na(brenner_contig_numbers)])

brenner_contigs_collapse <- paste0("c(", paste(brenner_contigs, collapse = ","), ")")

caeno_genomes[caeno_genomes$specific_name == "brenneri", "seqnames_exp"] <- paste0("c(paste0('Cbre_Contig', ", brenner_contigs_collapse, "))")

#### C briggsae ####

# briggsae then has chromosomes + contigs

brigg_contig_numbers <- as.numeric(str_remove(subFas_list$briggsae, '^>cb25\\.NA_'))[order(as.numeric(str_remove(subFas_list$briggsae, '^>cb25\\.NA_')))]
brigg_contigs <- collapse_number_ranges(numbervec = brigg_contig_numbers[!is.na(brigg_contig_numbers)])
brigg_contigs_collapse <- paste0("sprintf('%03d', c(", paste(brigg_contigs, collapse = ","), "))")

brigg_fpc_numbers <- as.numeric(str_remove(subFas_list$briggsae, '^>cb25\\.fpc'))[order(as.numeric(str_remove(subFas_list$briggsae, '^>cb25\\.fpc')))]
brigg_fpcs <- collapse_number_ranges(numbervec = brigg_fpc_numbers[!is.na(brigg_fpc_numbers)])
brigg_fpcs_collapse <- paste0("sprintf('%04d', c(", paste(brigg_fpcs, collapse = ","), "))")

caeno_genomes[caeno_genomes$specific_name == "briggsae", "seqnames_exp"] <- paste0("c(paste0('cb25.NA_', ", brigg_contigs_collapse, "),'I', 'II', 'III', 'IV', 'V', 'X', 'cb25.fpc2310b', paste0('cb25.fpc', ", brigg_fpcs_collapse, "))"
                                                                                   )
# briggsae has 1, II, III, IV, V, X? yes

#### C inopinata ####

caeno_genomes[caeno_genomes$specific_name == "inopinata", "seqnames_exp"] <- paste0("c(paste0('Sp34_Chr', 1:5), 'Sp34_ChrX')")

#### C latens ####

latens_contig_numbers <- as.numeric(str_remove(str_remove(subFas_list$latens, '^>scaffold_'), " .*"))[order(as.numeric(str_remove(str_remove(subFas_list$latens, '^>scaffold_'), " .*")))]

latens_contigs <- collapse_number_ranges(numbervec = latens_contig_numbers[!is.na(latens_contig_numbers)])

latens_contigs_collapse <- paste0("c(", paste(latens_contigs, collapse = ","), ")")

caeno_genomes[caeno_genomes$specific_name == "latens", "seqnames_exp"] <- paste0("c(paste0('scaffold_', ", latens_contigs_collapse, "))")
                        
#### C nigoni ####

# nigoni has 6 chromosomes plus a number of unplaced contigs. must do both
# also they all have .1 suffix which must be removed to make ranges but later replaced in expressions
nigoni_chrom_numbers <- as.numeric(str_remove(str_remove(subFas_list$nigoni, '^>CM00'), " .*"))[order(as.numeric(str_remove(str_remove(subFas_list$nigoni, '^CM00'), " .*")))]
nigoni_chrom_numbers <- nigoni_chrom_numbers[!is.na(nigoni_chrom_numbers)]
nigoni_chrom_numbers <- str_remove(nigoni_chrom_numbers, "\\.1")
nigoni_chroms <-collapse_number_ranges(as.numeric(nigoni_chrom_numbers))

nigoni_chroms_collapse <- paste0("c(", paste(nigoni_chroms, collapse = ","), ")")

nigoni_contig_numbers <- as.numeric(str_remove(str_remove(subFas_list$nigoni, '^>PDUG01000'), " .*"))[order(as.numeric(str_remove(str_remove(subFas_list$nigoni, '^PDUG01000'), " .*")))]
nigoni_contig_numbers <- nigoni_contig_numbers[!is.na(nigoni_contig_numbers)]
nigoni_contig_numbers <- str_remove(nigoni_contig_numbers, "\\.1")
nigoni_contigs <-collapse_number_ranges(as.numeric(nigoni_contig_numbers))

nigoni_contigs_collapse <- paste0("sprintf('%03d', c(", paste(nigoni_contigs, collapse = ","), "))")

caeno_genomes[caeno_genomes$specific_name == "nigoni", "seqnames_exp"] <- paste0("c(paste0('CM00', ", nigoni_chroms_collapse, ", '.1'), paste0('PDUG01000', ", nigoni_contigs_collapse, ",'.1'))")

#### C remanei ####

# all have same prefix
any(!str_detect(subFas_list$remanei, ">Crem_Contig"))

remanei_contig_numbers <- as.numeric(str_remove(str_remove(subFas_list$remanei, '^>Crem_Contig'), " .*"))[order(as.numeric(str_remove(str_remove(subFas_list$remanei, '^>Crem_Contig'), " .*")))]

remanei_contigs <-  collapse_number_ranges(numbervec = remanei_contig_numbers)

remanei_contigs_collapse <- paste0("c(", paste(remanei_contigs, collapse = ","), ")")

caeno_genomes[caeno_genomes$specific_name == "remanei", "seqnames_exp"] <- paste0("c(paste0('Crem_Contig', ", remanei_contigs_collapse, "))")

#### C sinica ####

# loadsa contigs here in this species....
# all 5 digits, needs sprintf

# all have the same prefix, which is Csp5_scaffold
any(!str_detect(subFas_list$sinica, ">Csp5_scaffold_"))

sinica_contig_numbers <- as.numeric(str_remove(str_remove(subFas_list$sinica, '^>Csp5_scaffold_'), " .*"))[order(as.numeric(str_remove(str_remove(subFas_list$sinica, '^Csp5_scaffold_'), " .*")))]

sinica_contigs <-  collapse_number_ranges(numbervec = sinica_contig_numbers)

sinica_contigs_collapse <- paste0("sprintf('%05d', c(", paste(sinica_contigs, collapse = ","), "))")

caeno_genomes[caeno_genomes$specific_name == "sinica", "seqnames_exp"] <- paste0("c(paste0('Csp5_scaffold_', ", sinica_contigs_collapse, "))")

#### C tribulationis ####

# loadsa contigs here in this species too....
# all 5 digits, needs sprintf

# all have the same prefix, which is CSP40.scaffold
any(!str_detect(subFas_list$tribulationis, ">CSP40.scaffold"))

tribulationis_contig_numbers <- as.numeric(str_remove(str_remove(subFas_list$tribulationis, '^>CSP40.scaffold'), " .*"))[order(as.numeric(str_remove(str_remove(subFas_list$tribulationis, '^CSP40.scaffold'), " .*")))]

tribulationis_contigs <-  collapse_number_ranges(numbervec = tribulationis_contig_numbers)

tribulationis_contigs_collapse <- paste0("sprintf('%05d', c(", paste(tribulationis_contigs, collapse = ","), "))")

caeno_genomes[caeno_genomes$specific_name == "tribulationis", "seqnames_exp"] <- paste0("c(paste0('CSP40.scaffold', ", tribulationis_contigs_collapse, "))")

#### C tropicalis ####

# all have the same prefix, which is >Scaffold
any(!str_detect(subFas_list$tropicalis, ">Scaffold"))

tropicalis_contig_numbers <- as.numeric(str_remove(str_remove(subFas_list$tropicalis, '^>Scaffold'), " .*"))[order(as.numeric(str_remove(str_remove(subFas_list$tropicalis, '^>Scaffold'), " .*")))]

tropicalis_contigs <-  collapse_number_ranges(numbervec = tropicalis_contig_numbers)

tropicalis_contigs_collapse <- paste0("c(", paste(tropicalis_contigs, collapse = ","), ")")

caeno_genomes[caeno_genomes$specific_name == "tropicalis", "seqnames_exp"] <- paste0("c(paste0('Scaffold', ", tropicalis_contigs_collapse, "))")

#### C zanzibari ####

# all have the same prefix, which is >CSP26.scaffold
# all 5 digits, needs sprintf
any(!str_detect(subFas_list$zanzibari, ">CSP26.scaffold"))

zanzibari_contig_numbers <- as.numeric(str_remove(str_remove(subFas_list$zanzibari, '^>CSP26.scaffold'), " .*"))[order(as.numeric(str_remove(str_remove(subFas_list$zanzibari, '^CSP26.scaffold'), " .*")))]

zanzibari_contigs <-  collapse_number_ranges(numbervec = zanzibari_contig_numbers)

zanzibari_contigs_collapse <- paste0("sprintf('%05d', c(", paste(zanzibari_contigs, collapse = ","), "))")

caeno_genomes[caeno_genomes$specific_name == "zanzibari", "seqnames_exp"] <- paste0("c(paste0('CSP26.scaffold', ", zanzibari_contigs_collapse, "))")

write.table(caeno_genomes,
            "input/caenorhabditis_genomes_seqnameexp.txt",
            sep = "\t",
            col.names = TRUE,
            row.names = FALSE )

