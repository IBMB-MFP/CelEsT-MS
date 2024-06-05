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

dir.create("~/Cel_GRN_manuscript/")
setwd("~/Cel_GRN_manuscript/")

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

caeno_genomes <- read.table("input/caenorhabditis_genomes_seqnameexp.txt",
                            sep = "\t",
                            header = TRUE)

subFas_list <- readRDS("output/subFas_list.rds")

#### BUILD SEED FILES ####

for(i in 1:nrow(caeno_genomes)){

  # according to Herve Pages, if we have many sequence contigs(e.g. in Brenneri) they need to go in mseqnames instead of seqnames. 
  # Thus I will do this for species with more than 50 seq names (which in fact is all but inopinata)
  if(length(subFas_list[[caeno_genomes[i, "specific_name"]]]) < 50){
  
  tempseed_vec <- c(str_replace_all(str_replace_all(paste0("Package: BSGenome.C", caeno_genomes[i, "specific_name"], ".Wormbase.", caeno_genomes[i, "Assembly"]), pattern = "-", replacement = ""), pattern = "_", replacement = ""),
                    paste0("Title: Full genomic sequences for Caenorhabditis ", caeno_genomes[i, "specific_name"], " (WormBase version ", caeno_genomes[i, "Assembly"], ")"),
                    paste0("Description: Full genomic sequences for Caenorhabditis ", caeno_genomes[i, "specific_name"], " (WormBase version ", caeno_genomes[i, "Assembly"], ")"),
                    "Version: 1.0",
                    paste0("organism: ", caeno_genomes[i, "Species"]),
                    paste0("common_name: ", caeno_genomes[i, "Species"]),
                    paste0("genome: ", caeno_genomes[i, "genome"]),
                    "provider: WormBase ParaSite",
                    "provider_version: WBPS18",
                    paste0("release_date:", caeno_genomes[i, "last_updated"]),
                    paste0("source_url: ", caeno_genomes[i, "source_url"]),
                    paste0("organism_biocview: ", str_replace(caeno_genomes[i, "Species"], pattern = " ", replacement = "_")),
                    paste0("BSgenomeObjname: C", caeno_genomes[i, "specific_name"]),
                    paste0("seqnames: ", caeno_genomes[i, "seqnames_exp"]),
                    "circ_seqs: character(0)",
                    paste0("seqs_srcdir: ~/Cel_GRN_orthology/input/sequences/", caeno_genomes[i, "specific_name"]))
  
  }
  
  if(length(subFas_list[[caeno_genomes[i, "specific_name"]]]) >= 50){

    tempseed_vec <- c(str_replace_all(str_replace_all(paste0("Package: BSGenome.C", caeno_genomes[i, "specific_name"], ".Wormbase.", caeno_genomes[i, "Assembly"]), pattern = "-", replacement = ""), pattern = "_", replacement = ""),
                      paste0("Title: Full genomic sequences for Caenorhabditis ", caeno_genomes[i, "specific_name"], " (WormBase version ", caeno_genomes[i, "Assembly"], ")"),
                      paste0("Description: Full genomic sequences for Caenorhabditis ", caeno_genomes[i, "specific_name"], " (WormBase version ", caeno_genomes[i, "Assembly"], ")"),
                      "Version: 1.0",
                      paste0("organism: ", caeno_genomes[i, "Species"]),
                      paste0("common_name: ", caeno_genomes[i, "Species"]),
                      paste0("genome: ", caeno_genomes[i, "genome"]),
                      "provider: WormBase ParaSite",
                      "provider_version: WBPS18",
                      paste0("release_date:", caeno_genomes[i, "last_updated"]),
                      paste0("source_url: ", caeno_genomes[i, "source_url"]),
                      paste0("organism_biocview: ", str_replace(caeno_genomes[i, "Species"], pattern = " ", replacement = "_")),
                      paste0("BSgenomeObjname: C", caeno_genomes[i, "specific_name"]),
                      "seqnames: character(0)",
                      paste0("mseqnames: ", caeno_genomes[i, "seqnames_exp"]),
                      "circ_seqs: character(0)",
                      paste0("seqs_srcdir: ~/Cel_GRN_orthology/input/sequences/", caeno_genomes[i, "specific_name"]))
    
  }

  write.table(tempseed_vec,
              file = paste0("input/sequences/", caeno_genomes[i, "specific_name"], "/C", caeno_genomes[i, "specific_name"], "_seed.dcf"),
              sep = "\n",
              row.names = FALSE,
              col.names = FALSE,
              quote = FALSE)

  message(paste0("Now forging data package for ", caeno_genomes[i, "specific_name"]))
  
  forgeBSgenomeDataPkg(paste0("input/sequences/", caeno_genomes[i, "specific_name"], "/C", caeno_genomes[i, "specific_name"], "_seed.dcf"),
                       destdir = paste0("output/BSGenome_pkgs/", caeno_genomes[i, "specific_name"]))
  
}