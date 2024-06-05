#start with empty workspace

rm(list = ls(all = TRUE))

# # clear all loaded packages
# invisible(lapply(paste0("package:", names(sessionInfo()$otherPkgs)),
#                  detach,
#                  character.only = TRUE, unload = TRUE))

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
              "conflicted",
              "BSGenome.Cbrenneri.Wormbase.Cbrenneri6.0.1b",
              "BSGenome.Cbriggsae.Wormbase.CB4",
              "BSGenome.Cinopinata.Wormbase.Sp34v7",
              "BSGenome.Clatens.Wormbase.CaeLat1.0",
              "BSGenome.Cnigoni.Wormbase.nigoni.pc2016.07.14",
              "BSGenome.Cremanei.Wormbase.Cremanei15.0.1",
              "BSGenome.Csinica.Wormbase.1.0",
              "BSGenome.Ctribulationis.Wormbase.CTRIBv1",
              "BSGenome.Ctropicalis.Wormbase.Caenorhabditissp11JU13733.0.1",
              "BSGenome.Czanzibari.Wormbase.CZANZv1"
              )

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

caeno_genomes[order(caeno_genomes$specific_name), "BSGenome_packagename"] <- c("BSGenome.Cbrenneri.Wormbase.Cbrenneri6.0.1b",
                                                             "BSGenome.Cbriggsae.Wormbase.CB4",
                                                             "BSGenome.Cinopinata.Wormbase.Sp34v7",
                                                             "BSGenome.Clatens.Wormbase.CaeLat1.0",
                                                             "BSGenome.Cnigoni.Wormbase.nigoni.pc2016.07.14",
                                                             "BSGenome.Cremanei.Wormbase.Cremanei15.0.1",
                                                             "BSGenome.Csinica.Wormbase.1.0",
                                                             "BSGenome.Ctribulationis.Wormbase.CTRIBv1",
                                                             "BSGenome.Ctropicalis.Wormbase.Caenorhabditissp11JU13733.0.1",
                                                             "BSGenome.Czanzibari.Wormbase.CZANZv1")

# define orthology relationships

parasite_mart <- useMart("parasite_mart", dataset = "wbps_gene", host = "https://parasite.wormbase.org", port = 443)

Cel_genes <- genes(TxDb.Celegans.UCSC.ce11.refGene)$gene_id

Cel_genes_BM <- getBM(mart = parasite_mart,
                      values = Cel_genes,
                      filters = "entrezgene_id",
                      attributes = c("wormbase_gseq",
                                     "entrezgene_id",
                                     "wormbase_gene",
                                     "wormbase_locus"
                      ))

Cel_genes_21urcensored <- Cel_genes_BM[!str_detect(Cel_genes_BM$wormbase_locus, "21ur") & !str_detect(Cel_genes_BM$wormbase_locus, "mir-"), ]

parasite_attributes_list <- list(brenneri = c("cabrenprjna20035_gene",
                  "cabrenprjna20035_orthology_type",
                  "cabrenprjna20035_homolog_perc_id"),
     inopinata = c("cainopprjdb5687_gene",
                   "cainopprjdb5687_orthology_type",
                   "cainopprjdb5687_homolog_perc_id"),
     tropicalis = c("catropprjna53597_gene",
                    "catropprjna53597_orthology_type",
                    "catropprjna53597_homolog_perc_id"),
     latens = c("calateprjna248912_gene",
                "calateprjna248912_orthology_type",
                "calateprjna248912_homolog_perc_id"),
     remanei = c("caremaprjna53967_gene",
                 "caremaprjna53967_orthology_type",
                 "caremaprjna53967_homolog_perc_id"),
     nigoni = c("canigoprjna384657_gene",
                "canigoprjna384657_orthology_type",
                "canigoprjna384657_homolog_perc_id"),
     briggsae = c("cabrigprjna10731_gene",
                  "cabrigprjna10731_orthology_type",
                  "cabrigprjna10731_homolog_perc_id"),
     sinica = c("casiniprjna194557_gene",
                "casiniprjna194557_orthology_type",
                "casiniprjna194557_homolog_perc_id"),
     zanzibari = c("cazanzprjeb12596_gene",
                   "cazanzprjeb12596_orthology_type",
                   "cazanzprjeb12596_homolog_perc_id"),
     tribulationis = c("catribprjeb12608_gene",
                       "catribprjeb12608_orthology_type",
                       "catribprjeb12608_homolog_perc_id"))

parasites_BM_list <- lapply(parasite_attributes_list, function(x){

  temp_orthology_BM <- getBM(mart = parasite_mart,
                                  values = Cel_genes_21urcensored$entrezgene_id,
                                  filters = "entrezgene_id",
                                  attributes = c("wormbase_gseq",
                                                 "entrezgene_id",
                                                 "wormbase_gene",
                                                 "wormbase_locus", x))
  
})

names(parasites_BM_list) <- names(parasite_attributes_list)

saveRDS(parasites_BM_list,
        "output/parasites_BM_list.rds")

# parasites_BM_list<- readRDS("output/parasites_BM_list.rds")

#### DOWNLOAD ANNOTATIONS ####

# Annotation links are found in the 'caeno_genomes' file

# create directory to deposit FASTA sequences

dir.create("output/orthologue_promoter_FASTA")

for(i in 1:nrow(caeno_genomes)){
    
temp_sp_specificname <-  caeno_genomes[i, "specific_name"] 
  
tempGFF_split <- unlist(str_split(caeno_genomes[i, "Annot_link"], "/"))
tempGFF_name <- tempGFF_split[length(tempGFF_split)]

if(!file.exists(paste0("input/", tempGFF_name))){
  
  download.file(url = caeno_genomes[i, "Annot_link"],
                destfile = paste0("input/", tempGFF_name))
  
}


# now can make the txdb

if(!file.exists(paste0("output/Txdb.C", temp_sp_specificname))){
  
  message(paste0("Creating ", temp_sp_specificname, " Txdb object now"))
  
  # first I need to read it in to get the sequence lengths. later I will make another txdb with the sequence information
  temp_sp_txdb <- makeTxDbFromGFF(file = paste0("input/", tempGFF_name),
                                  format = "gff3")
  
  temp_sp_seqlengths_vec <- sapply(temp_sp_txdb$user_seqlevels, function(x){

    templength <- length(get(caeno_genomes[i, "BSGenome_packagename"])[[paste0(x)]])
    
    if(templength > 1){
      return(templength)
    } else {
      return(width(get(caeno_genomes[i, "BSGenome_packagename"])[[paste0(x)]]))
    }
    
  })
  
  temp_sp_seqinfo <- Seqinfo(seqnames = temp_sp_txdb$user_seqlevels,
                             seqlengths = temp_sp_seqlengths_vec,
                             isCircular = rep(FALSE, times = length(temp_sp_seqlengths_vec)))
  
temp_sp_txdb <- makeTxDbFromGFF(file = paste0("input/", tempGFF_name),
                               format = "gff3",
                               chrominfo = temp_sp_seqinfo)

saveDb(temp_sp_txdb,
       file = paste0("output/Txdb.C", temp_sp_specificname))

} else { 
  
message(paste0("Loading ", temp_sp_specificname, " Txdb object now"))
  
temp_sp_txdb <- loadDb(paste0("output/Txdb.C", temp_sp_specificname))
  
}

# assign(paste0("Txdb.C", temp_sp_specificname, ".WBPS"), temp_sp_txdb)

message(paste0("Defining ", temp_sp_specificname, " promoters now"))

temp_sp_promoters <- trim(promoters(genes(temp_sp_txdb),
                                   upstream = 1000,
                                   downstream = 200))

# filter out promoters without full range

temp_sp_promoters <- temp_sp_promoters[width(temp_sp_promoters) == 1200]

# filter out ones that are not one to one orthologues of worm genes

temp_sp_onetoone_genes <- parasites_BM_list[[temp_sp_specificname]][parasites_BM_list[[temp_sp_specificname]][parasite_attributes_list[[temp_sp_specificname]][2]] == "ortholog_one2one", ]

temp_sp_onetoone_40pcid <- temp_sp_onetoone_genes[temp_sp_onetoone_genes[parasite_attributes_list[[temp_sp_specificname]][3]] > 40, ]

temp_sp_promoters_filt <- temp_sp_promoters[temp_sp_promoters$gene_id %in% unlist(temp_sp_onetoone_40pcid[parasite_attributes_list[[temp_sp_specificname]][1]]), ]

message(paste0("Getting ", temp_sp_specificname, " sequences now"))

temp_sp_promoters_filt_seq <- get_sequence(temp_sp_promoters_filt,
                                     get(caeno_genomes[i, "BSGenome_packagename"]))

names(temp_sp_promoters_filt_seq) <- temp_sp_promoters_filt$gene_id

message(paste0("Writing ", temp_sp_specificname, " sequences now"))

write.DNAStringSet(temp_sp_promoters_filt_seq,
                   format = "fasta",
                   filename = paste0("output/orthologue_promoter_FASTA/", temp_sp_specificname,"_one2onepromoters_seq.fasta"))

saveRDS(temp_sp_onetoone_40pcid[match(names(temp_sp_promoters_filt_seq), temp_sp_onetoone_40pcid[, parasite_attributes_list[[temp_sp_specificname]][1]]), ],
        paste0("output/C", temp_sp_specificname, "_one2one_promoterseq.rds"))

}

          