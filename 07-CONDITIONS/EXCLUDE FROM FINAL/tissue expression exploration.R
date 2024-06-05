library("openxlsx")
library("DESeq2")
library(biomaRt)
library(decoupleR)

parasite_mart <- useMart("parasite_mart", dataset = "wbps_gene", host = "https://parasite.wormbase.org", port = 443)


cao <- read.xlsx("~/Downloads/aam8940_cao_sm_tables_s1_to_s14.xlsx",
                 sheet = 3,
                 startRow = 2)



row.names(cao) <- cao[, 1]
cao <- cao[, 3:ncol(cao)]

cao[] <- lapply(cao, round)

col_tissue <- as.matrix(colnames(cao))

rownames(col_tissue) <- colnames(cao)
colnames(col_tissue) <- c("Tissue")

col_tissue[, "Tissue"] <- factor(as.character(col_tissue[, "Tissue"]))

# need to create a DESeq2 object. Design set to ~1 allows for use of estimateSizeFactors
tempdds <- DESeqDataSetFromMatrix(countData = cao, colData = col_tissue, design = ~ 1)
tempdds <- estimateSizeFactors(tempdds)

testdds <- DESeqDataSetFromMatrix(countData = cao, colData = col_tissue, design = ~Tissue)
testdds <- estimateSizeFactors(testdds)
colData(testdds)
testdds <- estimateDispersions(tempdds)

design(testdds) <- ~0+Tissue

testdds_test <- nbinomWaldTest(testdds)

contnumber = -1/6
contrastvec <- c(1, rep(contnumber, times = 6))
res <- results(testdds, contrast = contrastvec)


?DESeq

View(as.data.frame(res))

packageVersion("DESeq2")

genesBM <- getBM(attributes = c("wormbase_locus", "wormbase_gene", "wormbase_gseq"),
                 filter = "wbps_gene_id",
                 mart = parasite_mart,
                 values = row.names(cao))

cao_names <- genesBM[match(row.names(cao), genesBM$wormbase_gene), "wormbase_gseq"]
cao_names[duplicated(cao_names)|duplicated(cao_names, fromLast = TRUE)]


cao <- cao[-which(duplicated(cao_names)|duplicated(cao_names, fromLast = TRUE)), ]
cao_names <- cao_names[-which(duplicated(cao_names)|duplicated(cao_names, fromLast = TRUE))]



row.names(cao) <- cao_names

CelEsT <- read.table("~/Cel_GRN_manuscript/output/GRNs/allthree_equalweights.txt",
                     sep = "\t",
                     header = TRUE)
?decouple
cao_decouple <- decoupleR::decouple(mat = cao, 
                    network = CelEsT,
                    .source = "source",
                    .target = "target",
                    statistics = "mlm",
                    args = list(mlm = list(.mor = "weight")))


cao_decouple_mlm <- cao_decouple[cao_decouple$statistic == "mlm", ]

zscores <- sapply(unique(cao_decouple_mlm$source), function(x){

  thisone_vec <- unlist(cao_decouple_mlm[cao_decouple_mlm$source == x, "score"])
  names(thisone_vec) <- unlist(cao_decouple_mlm[cao_decouple_mlm$source == x, "condition"])
  
  (thisone_vec - mean(thisone_vec)) / sd(thisone_vec)
  
})

zscores <- data.frame(t(zscores))

rawscores <- sapply(unique(cao_decouple_mlm$source), function(x){
  
  thisone_vec <- unlist(cao_decouple_mlm[cao_decouple_mlm$source == x, "score"])
  names(thisone_vec) <- unlist(cao_decouple_mlm[cao_decouple_mlm$source == x, "condition"])

  thisone_vec
  
})

rawscores <- data.frame(t(rawscores))

#### eliminate lowly expressed TFs

maxTPMs <- apply(cao, 1, max)

expressed <- names(maxTPMs[maxTPMs > 10])

nonexpressedTFs <- unique(cao_decouple$source)[!unique(cao_decouple$source) %in% expressed]
zscores_cens <- zscores[!row.names(zscores) %in% nonexpressedTFs, ]

# convertnames <- genesBM[match(row.names(zscores_cens), genesBM$wormbase_gseq), "wormbase_locus"]
# convertnames[convertnames == ""] <- row.names(zscores_cens)[convertnames == ""]
# 
# row.names(zscores_cens) <- convertnames


# how often do z scores agree with expression specificity



colSums(do.call(rbind, lapply(colnames(zscores_cens), function(x){

  thisonevec <- zscores_cens[, x]
  names(thisonevec) <- row.names(zscores_cens)
  
  max_exp <- apply(cao[names(thisonevec[thisonevec > 1.75]), ], 1,  function(y){which(y == max(y))})
  
  c("match" = sum(colnames(zscores_cens)[max_exp] == x), "nomatch" = sum(colnames(zscores_cens)[max_exp] != x))
  
})))

colSums(do.call(rbind, lapply(colnames(zscores_cens), function(x){
  
  thisonevec <- zscores_cens[, x]
  names(thisonevec) <- row.names(zscores_cens)
  
  median_exp <- apply(cao[names(thisonevec[thisonevec > 1.75]), ], 1, median)
  
  c("match" = sum(cao[names(median_exp), x] > median_exp), "nomatch" = sum(cao[names(median_exp), x] <= median_exp))
  
})
))

colSums(do.call(rbind, lapply(colnames(zscores_cens), function(x){

  thisonevec <- zscores_cens[, x]
  names(thisonevec) <- row.names(zscores_cens)

  max_exp <- apply(cao[names(thisonevec[thisonevec < -1.75]), ], 1,  function(y){which(y == max(y))})

  c("match" = sum(colnames(zscores_cens)[unlist(max_exp)] == x), "nomatch" = sum(colnames(zscores_cens)[unlist(max_exp)] != x))
  
})))



library(Hmisc)

rcorr(as.matrix(cao))


