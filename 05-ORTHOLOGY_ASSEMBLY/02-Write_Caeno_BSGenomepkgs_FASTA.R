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
  
  # read in from file with subentries to compile list
  temp_sp_readL_list <- readToList.MFP(id = temp_sp_subfas,
                                text = temp_sp_readL,
                                con = paste0("input/", tempFASTA_name))
  
  # sort alphabetically
  temp_sp_readL_list2 <- temp_sp_readL_list[match(sort(sapply(temp_sp_readL_list, function(x){x[1]})), sapply(temp_sp_readL_list, function(x){x[1]}))]
  
  # remove chr / contig sgtart end from FASTA subheader
  temp_sp_readL_list3 <- lapply(temp_sp_readL_list2, function(x){
    
    x[1] <- str_remove(x[1], " .*")
    return(x)
    
    })
  
  # name list by subheader names
  names(temp_sp_readL_list3) <- sapply(temp_sp_readL_list3, function(x){x[1]})
  
  # make folder in input for sequences and for the data package
  dir.create(paste0("input/sequences/", caeno_genomes[i, "specific_name"]))
  dir.create(paste0("output/BSGenome_pkgs/", caeno_genomes[i, "specific_name"]))
  
  # write sequences to separate FASTAs
  for(j in 1:length(temp_sp_readL_list3)){

    writeLines(text = temp_sp_readL_list3[[j]],
               con = paste0("input/sequences/", caeno_genomes[i, "specific_name"], "/", str_remove(names(temp_sp_readL_list3)[j], "^>"), ".fa"),
               sep = "\n")
    
  }
 
}
