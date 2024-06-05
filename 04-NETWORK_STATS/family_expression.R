
library(stringr)
library(DESeq2)

MODern_allTFs_data <- read.xlsx("input/Kudron2024_biorXiv_Supp1.xlsx",
                                sheet = 6)

# bcl-11 is missing correct sequence name

MODern_allTFs_data[MODern_allTFs_data$Sequence_name == "bcl-11", "Sequence_name"] <- "F13H6.1"

wTF3 <- read.xlsx("input/WTF3.xlsx")

wTF3[, "simplified_fam"] <- str_remove(str_remove(str_remove(str_remove(wTF3$DBD, " - [0-9]{1,2} finger.*$"), " - [0-9]{1,2} domain.*$"), " x[0-9]{1,2}$"), " $")

allTFs <- unique(c(MODern_allTFs_data$Sequence_name, wTF3$Sequence.name))

fullset_TFs_BM <- readRDS("output/fullset_TFs_BM.rds")

fam_table <- table(wTF3$simplified_fam)
fam_table[order(fam_table, decreasing = TRUE)]

combo_GRN <- read.table("output/GRNs/allthree_equalweights.txt",
                        header = TRUE)

fullset_TFs_BM [match(unique(combo_GRN$source)[!unique(combo_GRN$source) %in% MODern_allTFs_data$Sequence_name], fullset_TFs_BM $wormbase_gseq), ]

MODern_allTFs_data[MODern_allTFs_data$Sequence_name %in% missingTFs, ]



nonCelEsT_TFs <- unique(c(MODern_allTFs_data[!MODern_allTFs_data$Sequence_name %in% CelEsT_TFs, "Sequence_name"],
                          wTF3$Sequence.name[!wTF3$Sequence.name %in% CelEsT_TFs]))

CelEsT_df <- data.frame("wormbase_gseq" = CelEsT_TFs,
                        "family" = wTF3[match(CelEsT_TFs, wTF3$Sequence.name), "simplified_fam"])

nonCelEsT_df <- data.frame("wormbase_gseq" = nonCelEsT_TFs)
                           
nonCelEsT_df[!is.na(nonCelEsT_df$wormbase_gseq), "family"] <- wTF3[match(nonCelEsT_df$wormbase_gseq[!is.na(nonCelEsT_df$wormbase_gseq)], wTF3$Sequence.name), "simplified_fam"]



table(CelEsT_df$family)[order(table(CelEsT_df$family), decreasing = TRUE)]
table(nonCelEsT_df$family)[order(table(nonCelEsT_df$family), decreasing = TRUE)]



##### EXPRESSION #####


library("stringr")

all_series_filelist <- list.files("~/Cel_GRN_exploration/input/dev_timecourse_counts/")

# change colnames to description
all_sra_dev <- read.table("~/Cel_GRN_exploration/all_dev_timecourse_SRA.txt",
                          sep = "\t",
                          header = TRUE)

meeuse_2020_sra_dev <- all_sra_dev[all_sra_dev$Batch %in% c("meeuse_late", "meeuse_early"), ]

meeuse_filelist <- all_series_filelist[str_remove(all_series_filelist, "_FCcounts.txt") %in% meeuse_2020_sra_dev$Run]

meeuse_counts_list <- lapply(meeuse_filelist, function(x){
  
  read.table(paste0("~/Cel_GRN_exploration/input/dev_timecourse_counts/", x),
             sep = "\t",
             skip = 1,
             header = TRUE)
  
})

meeuse_counts_mat <- sapply(meeuse_counts_list, function(x){
  
  x[, 7]
  
})

colnames(meeuse_counts_mat) <- str_remove(meeuse_filelist, "_FCcounts.txt")
row.names(meeuse_counts_mat) <- meeuse_counts_list[[1]]$Geneid

meeuse_2020_sra_times <- str_remove(meeuse_2020_sra_dev$Developmental_stage, "wild-type \\(N2\\) at [0-9]{2}C\\: ")
meeuse_2020_sra_times <- str_remove(meeuse_2020_sra_times, " hours$")
meeuse_2020_sra_times <- str_remove(meeuse_2020_sra_times, " hrs$")
meeuse_2020_sra_times <- str_remove(meeuse_2020_sra_times, " h$")
meeuse_2020_sra_times <- str_remove(meeuse_2020_sra_times, "h$")

names(meeuse_2020_sra_times) <- meeuse_2020_sra_dev$Run

colnames(meeuse_counts_mat) <- meeuse_2020_sra_times[colnames(meeuse_counts_mat)]

colnames(meeuse_counts_mat)[duplicated(colnames(meeuse_counts_mat))] <- paste0(colnames(meeuse_counts_mat)[duplicated(colnames(meeuse_counts_mat))], "b")

# DESeq2 wants a colData object. Not actually used for the normalisation. Here we can use the sample IDs with tissue.
dev_col_data <- as.matrix(colnames(meeuse_counts_mat))
colnames(dev_col_data) <- "time"
row.names(dev_col_data) <- colnames(meeuse_counts_mat)

meeuse_counts_mat_dds <- DESeqDataSetFromMatrix(countData = meeuse_counts_mat, colData = dev_col_data, design = ~ 1)

# this function estimates the scaling factors from the samples from the median of ratios wrt to the geometric mean for each gene across samples
meeuse_counts_mat_dds <- estimateSizeFactors(meeuse_counts_mat_dds)

# put the counts normalised by the scaling factors in a new object
meeuse_normalised_counts <- counts(meeuse_counts_mat_dds, normalized = TRUE)

dev_BM <- getBM(attributes = c("wormbase_gene",
                               "wormbase_gseq",
                               "wormbase_locus"),
                filters = "wbps_gene_id",
                values = row.names(meeuse_normalised_counts),
                mart = parasite_mart)

templook <- dev_BM[match(row.names(meeuse_normalised_counts), dev_BM$wormbase_gene), "wormbase_gseq"]
templook[(duplicated(templook)|duplicated(templook, fromLast = TRUE))]

meeuse_normalised_counts <- data.frame(meeuse_normalised_counts)
meeuse_normalised_counts[, "gseq"] <- dev_BM[match(row.names(meeuse_normalised_counts), dev_BM$wormbase_gene), "wormbase_gseq"]

meeuse_normalised_counts <- meeuse_normalised_counts[!is.na(meeuse_normalised_counts$gseq), ]

# duplicated genes. stops setting gseq as row.names
# weirdly it does seem to be two separate genes with different locus ids and WBgene IDs on Wormbasd but with the same gseq ID. 
# anyway, get rid of both

meeuse_normalised_counts <- meeuse_normalised_counts[!(duplicated(meeuse_normalised_counts$gseq)|duplicated(meeuse_normalised_counts$gseq, fromLast = TRUE)), ]

row.names(meeuse_normalised_counts) <- meeuse_normalised_counts$gseq

# get rid of gseq column now
meeuse_normalised_counts <- meeuse_normalised_counts[, -ncol(meeuse_normalised_counts)]

CelEsT_TFs_normcounts <- meeuse_normalised_counts[CelEsT_TFs, ]

nonCelEsT_TFs_normcounts <- meeuse_normalised_counts[nonCelEsT_TFs, ]

CelEsT_TFs_max <- apply(CelEsT_TFs_normcounts, 1, max)
nonCelEsT_TFs_max <- apply(nonCelEsT_TFs_normcounts, 1, max)

max_plot_df <- rbind(data.frame("max" = log1p(CelEsT_TFs_max),
                               "group" = "YES"),
                     data.frame("max" = log1p(nonCelEsT_TFs_max),
                               "group" = "NO"))

max_plot_df$group <- factor(max_plot_df$group, levels = c("YES", "NO"))

shapiro.test(log10(CelEsT_TFs_max + 1))
shapiro.test(log10(nonCelEsT_TFs_max + 1))

t.test(unlist(log10(CelEsT_TFs_max + 1)), unlist(log10(nonCelEsT_TFs_max + 1)))$p.value

ggplot(max_plot_df, aes(x = group, y = max, fill = group)) + 
  geom_violin() + 
  scale_fill_manual(values = c("purple", "yellow")) +
  geom_jitter(width = 0.1,
              alpha = 0.2,
              colour = "grey") + 
  geom_boxplot(width = 0.07, 
               outlier.shape = NA,
               alpha = 0) +
  theme_classic() + 
  theme(legend.position = "none",
        axis.text = element_text(colour = "black")) + 
  xlab("TF in CelEsT?") + 
  ylab("Max."~log[10]~"expression + 1 (pseudocounts)")

boxplot(log10(CelEsT_TFs_max), log10(nonCelEsT_TFs_max))


CelEsT_TFs_mean <- apply(CelEsT_TFs_normcounts, 1, mean)
nonCelEsT_TFs_mean <- apply(nonCelEsT_TFs_normcounts, 1, mean)

boxplot(unlist(log10(CelEsT_TFs_mean + 1)), log10(nonCelEsT_TFs_mean + 1))
t.test(unlist(log10(CelEsT_TFs_mean + 1)), unlist(log10(nonCelEsT_TFs_mean + 1)))


CelEsT_TFs_mean[order(CelEsT_TFs_mean)]

mean_plot_df <- rbind(data.frame("mean" = log1p(CelEsT_TFs_mean),
                                "group" = "YES"),
                     data.frame("mean" = log1p(nonCelEsT_TFs_mean),
                                "group" = "NO"))

mean_plot_df$group <- factor(mean_plot_df$group, levels = c("YES", "NO"))

ggplot(mean_plot_df, aes(x = group, y = mean, fill = group)) + 
  geom_violin() + 
  scale_fill_manual(values = c("purple", "yellow")) +
  geom_jitter(width = 0.1,
              alpha = 0.2,
              colour = "grey") + 
  geom_boxplot(width = 0.07, 
               outlier.shape = NA,
               alpha = 0) +
  theme_classic() + 
  theme(legend.position = "none",
        axis.text = element_text(colour = "black")) + 
  xlab("TF in CelEsT?") + 
  ylab("Mean."~log[10]~"expression + 1 (pseudocounts)")
