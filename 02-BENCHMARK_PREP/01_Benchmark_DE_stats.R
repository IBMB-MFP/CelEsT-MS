#start with empty workspace

rm(list = ls(all = TRUE))                                               

# # clear all loaded packages
# invisible(lapply(paste0("package:", names(sessionInfo()$otherPkgs)),
#                  detach,
#                  character.only = TRUE, unload = TRUE))

# turn off scientific notation for plots

options(scipen = 10000)

#### set working directory ####

# here create new folder and set working directory within it

setwd("~/Cel_GRN_manuscript/")

dir.create("input")
dir.create("output")

#### DEFINE FUNCTIONS ####

# Kudos to https://stackoverflow.com/questions/18813526/check-whether-all-elements-of-a-list-are-in-equal-in-r
all.identical.list <- function(l) identical(unname(l[-length(l)]), unname(l[-1]))

chip_or_fimo_subset <- function(cutoff = 100,
                                targetlist,
                                TFlist = NULL,
                                add_sign = FALSE){
  
  if(!is.null(TFlist)){
  targetlist_benchmarkTFs <- targetlist[names(targetlist) %in% TFlist]
  } else {
    targetlist_benchmarkTFs <- targetlist
  }
  
  targetlist_benchmarkcutoff <- reshape2::melt(lapply(targetlist_benchmarkTFs, function(x){x[1:min(cutoff, length(x))]}))
  names(targetlist_benchmarkcutoff) <- c("target", "source")
  
  # remove autoregulation
  targetlist_benchmarkcutoff <- targetlist_benchmarkcutoff[!(targetlist_benchmarkcutoff$target == targetlist_benchmarkcutoff$source), ]
  
  # prepend X to targets with gseq ID beginning with a number 
  # (because R will have prepended it automatically in expression data where gseq IDs are column names)
  targetlist_benchmarkcutoff[str_detect(targetlist_benchmarkcutoff$target, "^[0-9]"), "target"] <- paste0("X", targetlist_benchmarkcutoff[str_detect(targetlist_benchmarkcutoff$target, "^[0-9]"), "target"])
  
  targetlist_benchmarkcutoff[, "weight"] <- 1
  
  if(isTRUE(add_sign)){
    
    for(i in 1:length(negative_regulators)){
      
      targetlist_benchmarkcutoff[targetlist_benchmarkcutoff$TF == negative_regulators[i], "weight"] <- -1
      
    }
    
  }
  
  targetlist_benchmarkcutoff <- targetlist_benchmarkcutoff[, c("source", "target", "weight")]
  
  return(targetlist_benchmarkcutoff)
  
}

randomise_GRN <- function(GRN){
  
  all_unique_targets <- unique(GRN$target)
  
  temp_TF_randomlist <- lapply(unique(GRN$source), function(thisTF){
    
    thisTF_no_of_targets <- nrow(GRN[GRN$source == thisTF, ])
    
    random_sample <- sample(all_unique_targets, 
                            size = thisTF_no_of_targets,
                            replace = FALSE)
    
    random_weights <- sample(GRN$weight, 
                             size = thisTF_no_of_targets,
                             replace = FALSE)
    
    data.frame("source" = thisTF,
               "target" = random_sample,
               "weight" = random_weights)
    
  })
  
  do.call(rbind, temp_TF_randomlist)
  
}

calc.score.by.motif <- function(motif_list = separate_by_motif,
                                homotypic_binding = TRUE,
                                returns_parameter = 3) {
  
  
  outlist <- lapply(motif_list, function(x){
    
    sapply(unique(x$sequence_name), function(y){
      
      thisonedata <- x[x$sequence_name == y, ]
      
      n <- nrow(thisonedata)
      
      if(n == 1){
        
        return(thisonedata$score)
        
      } else {
        
        if(homotypic_binding == TRUE){
          
          thisonedata <- thisonedata[order(thisonedata$score, decreasing = TRUE), ]
          
          sumscore <- 0
          
          for(i in 1:n){
            
            thisscore <- thisonedata[i, 'score']
            sumscore <- sumscore + thisscore * (1 / returns_parameter^(i-1))
            
          }
          
          return(sumscore)
          
        } else {
          
          return(thisonedata$score[1])
          
        }
        
      }
      
    })
    
  })
  
  outlist_order <- lapply(outlist, function(x){
    
    x[order(x, decreasing = TRUE)]
    
  })
  
  return(outlist_order)
  
}

# KUDOS to Kamil from biostars [https://www.biostars.org/p/171766/]
counts_to_tpm <- function(counts, featureLength, meanFragmentLength) {

  # Ensure valid arguments.
  stopifnot(length(featureLength) == nrow(counts))
  stopifnot(length(meanFragmentLength) == ncol(counts))
  
  # Compute effective lengths of features in each library.
  effLen <- do.call(cbind, lapply(1:ncol(counts), function(i) {
    featureLength - meanFragmentLength[i] + 1
  }))

  # Exclude genes with length less than the mean fragment length.
  idx <- apply(effLen, 1, function(x) min(x) > 1)
  temp_counts <- counts[idx,]
  temp_effLen <- effLen[idx,]
  temp_featureLength <- featureLength[idx]

  # Process one column at a time.
  tpm <- do.call(cbind, lapply(1:ncol(temp_counts), function(i) {
    rate = log(temp_counts[,i]) - log(temp_effLen[,i])
    denom = log(sum(exp(rate)))
    exp(rate - denom + log(1e6))
  }))
  
  # Copy the row and column names from the original matrix.
  colnames(tpm) <- colnames(temp_counts)
  rownames(tpm) <- rownames(temp_counts)
  
  return(tpm)
  
}

#### LOAD PACKAGES & FUNCTIONS ####

## First specify the packages of interest

packages <- c("stringr",
              "biomaRt",
              "XML",
              "openxlsx",
              "edgeR",
              "splines")

## Now load or install&load all
package.check <- lapply(packages, function(y){
  
  if (!require(y, character.only = TRUE)) {
    
    install.packages(y, dependencies = TRUE)
    
    library(y, character.only = TRUE)
    
  }
  
})

# Here define packages which need to be loaded through biocmanager

biocmanager_packages <- c("DESeq2",
                          "decoupleR") 

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
                     "LBMC/RAPToR",
                     "LBMC/wormRef")

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

parasite_mart <- useMart("parasite_mart", dataset = "wbps_gene", host = "https://parasite.wormbase.org", port = 443)

benchmark_SRA <- read.table("input/benchmark_SRA_all.txt",
                            header = TRUE,
                            sep = "\t")

benchmark_SRA_BM <- getBM(mart = parasite_mart,
                          values = benchmark_SRA$target_name,
                          filters = "wormbase_locusid",
                          attributes = c("wbps_gene_id",
                                         "wormbase_gseq",
                                         "wormbase_locus"))

benchmark_SRA$target_gseq <- benchmark_SRA_BM[match(benchmark_SRA$target_name, benchmark_SRA_BM$wormbase_locus), "wormbase_gseq"]

benchmark_SRA[, "ExperimentNo"] <- paste0("Experiment", sprintf("%02d", 1:nrow(benchmark_SRA)))

#### PREPARE COUNTS FOR BENCHMARK SET (SRA EXPERIMENTS) ####

FCfilelist <- list.files("input/benchmark_counts/")

counts_list <- lapply(FCfilelist, function(x){
  
  read.table(paste0("input/benchmark_counts/", x),
             sep = "\t",
             skip = 1,
             header = TRUE)
  
})

counts_mat <- sapply(counts_list, function(x){
  
  x[, 7]
  
})

colnames(counts_mat) <- str_remove(FCfilelist, "_FCcounts.txt")
row.names(counts_mat) <- counts_list[[1]]$Geneid

#### RAW DE STATISTICS ####

allstudiesDE <- apply(benchmark_SRA, 1, function(thisstudy){
  
  ## For DEBUGGING
  # thisstudy <- unlist(benchmark_SRA[benchmark_SRA$Accession == "GSE173580", ][1, ])
  
  tryCatch({
    print(thisstudy["Accession"])
    
    control_samples <- unlist(str_split(thisstudy["control_samples"], ", "))
    treatment_samples <- unlist(str_split(thisstudy["treatment_samples"], ", "))
    
    thisdata <- counts_mat[, c(control_samples, treatment_samples)]
    
    tempcoldata <- data.frame(group = c(rep("control", times = length(control_samples)), rep("treatment", times = length(treatment_samples))))
    row.names(tempcoldata) <- c(control_samples, treatment_samples)
    
    tempdds <-  DESeqDataSetFromMatrix(countData = thisdata, colData = tempcoldata, design = ~ group)
    
    tempdds <- DESeq(tempdds)
    
    temp_res <- results(tempdds)
    
    return(as.data.frame(temp_res))
    
  }, error = function(e){
    
    message(paste0(thisstudy["Accession"], " has produced an error"))
    
    return(NULL)
    
  })
  
})

t_values <- lapply(allstudiesDE, function(x){x[, "stat"]})

names(t_values) <- benchmark_SRA$ExperimentNo

t_values <- t_values[!sapply(t_values, is.null)]

t_values <- as.data.frame(t_values)

t_values_t <- data.frame(t(t_values))
colnames(t_values_t) <- row.names(allstudiesDE[[1]])

# convert to gseq
genenameBM <- getBM(mart = parasite_mart,
                    values = colnames(t_values_t),
                    filters = "wbps_gene_id",
                    attributes = c("wbps_gene_id",
                                   "wormbase_gseq"))

colnames(t_values_t) <- genenameBM[match(colnames(t_values_t), genenameBM$wbps_gene_id), "wormbase_gseq"]
t_values_t <- t_values_t[, !is.na(colnames(t_values_t))]

observations <- benchmark_SRA
row.names(observations) <- observations$ExperimentNo

OEexperiments <- observations[observations$mode == "overexpression", "ExperimentNo"]

t_values_OE <- t_values_t
t_values_OE[OEexperiments, ] <- -(t_values_OE[OEexperiments, ])

# t_values_OE[is.na(t_values_OE)] <- 0
t_values_OE[is.na(t_values_OE)] <- NaN


write.table(t_values_OE,
            "output/benchmark_DEstats_RAW.txt",
            sep = "\t",
            row.names = TRUE,
            col.names = TRUE)

write.table(observations,
            "output/benchmark_observations.txt",
            sep = "\t",
            row.names = TRUE,
            col.names = TRUE)

observations_for_TS1 <- observations
observations_for_TS1 <- observations_for_TS1[!str_detect(observations_for_TS1$Note, "abi"), ]

colnames(observations_for_TS1) <- str_replace_all(colnames(observations_for_TS1), "target", "TF")

row.names(observations_for_TS1) <- NULL

write.xlsx(observations_for_TS1,
           "output/observations_for_TS1.xlsx",
           )

# how many samples is it?

length(unlist(str_split(str_remove_all(observations_for_TS1$control_samples, " "), ","))) +
length(unlist(str_split(str_remove_all(observations_for_TS1$treatment_samples, " "), ",")))

# 301 controls
# 296 treatment
# 597 total





#### AGE CORRECTION WITH RAPToR ####

# prepare refs

ref_vec <- c("Cel_embryo",
             "Cel_larval",
             "Cel_larv_YA",
             "Cel_YA_1",
             "Cel_YA_2",
             "Cel_aging_1")

refdata_list <- lapply(ref_vec, function(x){
  
  prepare_refdata(x, 'wormRef', 5000)
  
})

names(refdata_list) <- ref_vec

allstudiesDE_RAPToR <- apply(benchmark_SRA, 1, function(thisstudy){

  ## FOR DEBUGGING
  # thisstudy <- unlist(benchmark_SRA[90, ])
  
  tryCatch({
    
    print(thisstudy["Accession"])
    
    control_samples <- unlist(str_split(thisstudy["control_samples"], ", "))
    treatment_samples <- unlist(str_split(thisstudy["treatment_samples"], ", "))
    
    thisdata <- counts_mat[, c(control_samples, treatment_samples)]

    tempcoldata <- data.frame(group = c(rep("control", times = length(control_samples)), rep("treatment", times = length(treatment_samples))))
    row.names(tempcoldata) <- c(control_samples, treatment_samples)

    # First I need to convert the counts to TPM 
    
    my_temp_lengths <- sapply(row.names(thisdata), function(x){

      mean(wormRef::Cel_genes[wormRef::Cel_genes$wb_id == x, "transcript_length"])
      
    })

    my_temp_lengths <- my_temp_lengths[!is.na(my_temp_lengths)]
    
    gene_there <- base::intersect(names(my_temp_lengths), row.names(thisdata))
    
    thisdata_TPM <- counts_to_tpm(counts = thisdata[gene_there, ],
                                  featureLength = my_temp_lengths[gene_there],
                                  meanFragmentLength = as.numeric(rep(thisstudy["AvgSpotLen"], times = ncol(thisdata))))

    # remove genes with zero expression across all samples
    thisdata_TPM_expressed <- thisdata_TPM[apply(thisdata_TPM, 1, function(x){any(x != 0)}), ]
    
    # transform for later use with RAPToR as per RAPToR vignette
    thisdata_TPM_expressed_norm <- limma::normalizeBetweenArrays(thisdata_TPM_expressed, method = "quantile")
    thisdata_TPM_expressed_norm_log <- log1p(thisdata_TPM_expressed_norm) # log1p(x) = log(x + 1)
    
    dev_stage <- thisstudy["dev_stage"]
    
    if(is.na(dev_stage)){
      
      my_ref <- refdata_list[["Cel_larv_YA"]]
      
    } else {
      
      if(dev_stage == "embryo"){
        
        my_ref <- refdata_list[["Cel_embryo"]]
        
      }
      
      if(dev_stage == "L1"){
        
        my_ref <- refdata_list[["Cel_larval"]]
        
      }
      
      if(dev_stage %in% c(paste0("L", 2:4), "mixed_stage", "unclear", "Larv", "L4/YA")){
        
        my_ref <- refdata_list[["Cel_larv_YA"]]
        
      }
      
      if(dev_stage %in% c("YA", paste0("Ad_D", 0:2))) {
        
        my_ref <- refdata_list[["Cel_YA_2"]]
        
      }
      
      if(dev_stage %in% c(paste0("Ad_D", 3:15))) {
        
        my_ref <- refdata_list[["Cel_aging_1"]]
        
      }
      
      # as we don't have an appropriate reference to stage dauer samples, we will exempt them from this treatment
      if(dev_stage == "dauer"){
        
        tempdds <-  DESeqDataSetFromMatrix(countData = thisdata, colData = tempcoldata, design = ~ group)
        
        tempdds <- DESeq(tempdds)
        
        temp_res <- results(tempdds)
        
        return(as.data.frame(temp_res))
        
      }
      
    } # end of else
    
    temp_sample_ae <- RAPToR::ae(samp = thisdata_TPM_expressed_norm_log,                         # input gene expression matrix
                                 refdata = my_ref)
    
    if(dev_stage == "embryo"){
      interpolation_offset <- 60
    } else {
      interpolation_offset <- 1
    }
    
    interpolation_range <- range(temp_sample_ae$age.estimates[, 1]) + c(-interpolation_offset, interpolation_offset)
    
    interpolation_index <- get_refTP(my_ref, ae_values = interpolation_range)
    interpolation_index <- interpolation_index[1]:interpolation_index[2]
    
    interpolated_ref_expression <- get_refTP(my_ref,
                                             ae_values = interpolation_range,
                                             return.idx = FALSE)
    
    interp_time <- my_ref$time[interpolation_index]
    interp_tpm <- my_ref$interpGE[, interpolation_index]
    
    wormRef_gene_lengths <- sapply(row.names(interp_tpm), function(x){
      
      mean(wormRef::Cel_genes[wormRef::Cel_genes$wb_id == x, "transcript_length"])
      
    })
    
    libsize <- 25e6

    interp_count <- t((t (exp(interp_tpm) - 1)/ 1e6) * (libsize / median(wormRef_gene_lengths))) * wormRef_gene_lengths
    interp_count[interp_count < 0] <- 0
    interp_count <- round(interp_count)

    # now I have the age estimates. I need to do GLM with edgeR with glmFit. Include batch (sample vs reference data), variable of interest (mutation, group reference with control) and developmental time modelled with splines ( )
    # may need to determine optimal number of spline df by fitting models and using elbow plot (with automatic knee finder)
    
    # inital step; fit with fixed df for splines.
    # then figure out how to try multiple models and choose optimal.
    
    # ok now will try to implement this having understood it better

    df_SSQ <- sapply(1:8, function(df_param){
      
      sum(residuals(lm(t(interp_tpm) ~ splines::ns(interp_time, df = df_param))) ^2)
      
    })

    # here find the point at which the curve plateaus appropriately
    # threshold is arbitrary but 0.1 seems too high.
    threshold <- 0.01
    diff_1 <- diff(df_SSQ)
    
    temp_plateau <- which.max((diff_1 / diff_1[1]) < threshold) 
    
    ## Having found the plateau, proceed to fit model with appropriate spline df 
    
    thisdata <- thisdata[apply(thisdata, 1, max) > 5, ]
    
    combine_count <- do.call(cbind, format_to_ref(thisdata, interp_count)[1:2])
    
    combine_coldata <- data.frame(time = c(temp_sample_ae$age.estimates[, 1], interp_time),
                                  strain = c(tempcoldata[row.names(temp_sample_ae$age.estimates), 1], rep("control", ncol(interp_count))),
                                  batch = rep(c("sample", "reference"), c(ncol(thisdata), ncol(interp_count))))
    
    combine_coldata$strain <- factor(combine_coldata$strain, levels = c("control", "treatment"))
    combine_coldata$batch <- factor(combine_coldata$batch, levels = c("sample", "reference"))
    
    # here this was missing from the vignette but I think its necessary)
    row.names(combine_coldata) <- c(row.names(temp_sample_ae$age.estimates), colnames(interp_count))
    
    # estimate dispersions on sample data alone
    
    dd0 <- DESeqDataSetFromMatrix(countData = combine_count[, 1:length(c(control_samples, treatment_samples))],
                                  colData = combine_coldata[1:length(c(control_samples, treatment_samples)), ],
                                  design = ~ strain)

    dd0 <- estimateSizeFactors(dd0)
    dd0 <- estimateDispersions(dd0, fitType = "local") 
    d0 <- dispersions(dd0) # store dispersions
    d0[is.na(d0)] <- 0 # remove NAs
    
    formula <- paste0("~ splines::ns(time, df = ", temp_plateau, ") + batch + strain")
    
    # I have the problem that if I use a straight numerical number 3 for the spline df, it works, but if I feed it an object that equals 3 (ie temp_plateau), it doesnt work.
    # ok this appears to be solved by defining formula and using as.formula
    dd1 <-  DESeqDataSetFromMatrix(
      countData = combine_count,
      colData = combine_coldata,
      design = as.formula(formula) # use df found above
    )
    
    dd1 <- estimateSizeFactors(dd1)
    # inject dispersions from sample-only model
    dispersions(dd1) <- d0
    
    dd1 <- nbinomWaldTest(dd1)

    return(results(dd1, contrast = c("strain", "treatment", "control")))

  }, error = function(e){
    
    message(paste0(thisstudy["Accession"], " has produced an error"))
    
    return(NULL)
    
  })
  
})

names(allstudiesDE_RAPToR) <- benchmark_SRA$ExperimentNo

allstudiesDE_RAPToR <- allstudiesDE_RAPToR[!sapply(allstudiesDE_RAPToR, is.null)]

allstudiesDE_stats <- lapply(allstudiesDE_RAPToR, function(x){

  tempout <- as.data.frame(x)[, "stat"]
  names(tempout) <- row.names(as.data.frame(x))
  
  tempout
  
})

allgenesthereDE_new <- unique(unlist(lapply(allstudiesDE_stats, names)))

allstudiesDE_RAPToR_df <- data.frame(lapply(allstudiesDE_stats, function(x){
  
  tempvec <- x[match(allgenesthereDE_new, names(x))]
  names(tempvec) <- allgenesthereDE_new
  
  tempvec
  
}))

allstudiesDE_RAPToR_t <- data.frame(t(allstudiesDE_RAPToR_df))

# convert to gseq
genenameBM <- getBM(mart = parasite_mart,
                    values = colnames(allstudiesDE_RAPToR_t),
                    filters = "wbps_gene_id",
                    attributes = c("wbps_gene_id",
                                   "wormbase_gseq"))

colnames(allstudiesDE_RAPToR_t) <- genenameBM[match(colnames(allstudiesDE_RAPToR_t), genenameBM$wbps_gene_id), "wormbase_gseq"]
allstudiesDE_RAPToR_t <- allstudiesDE_RAPToR_t[, !is.na(colnames(allstudiesDE_RAPToR_t))]

row.names(benchmark_SRA) <- benchmark_SRA$ExperimentNo
observations_RAPToR <- benchmark_SRA[, c(ncol(benchmark_SRA), 1:(ncol(benchmark_SRA)-1))]

OEexperiments_RAPToR <- observations_RAPToR[observations_RAPToR$mode == "overexpression", "ExperimentNo"]

allstudiesDE_RAPToR_OE <- allstudiesDE_RAPToR_t
allstudiesDE_RAPToR_OE[OEexperiments_RAPToR, ] <- -(allstudiesDE_RAPToR_OE[OEexperiments_RAPToR, ])

# allstudiesDE_RAPToR_OE[is.na(allstudiesDE_RAPToR_OE)] <- NaN
allstudiesDE_RAPToR_OE[is.na(allstudiesDE_RAPToR_OE)] <- 0

write.table(allstudiesDE_RAPToR_OE,
            "output/benchmark_DEstats_RAPToR.txt",
            sep = "\t",
            row.names = TRUE,
            col.names = TRUE)



