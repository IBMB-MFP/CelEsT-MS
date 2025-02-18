#start with empty workspace

rm(list = ls(all = TRUE))

# turn off scientific notation for plots

options(scipen=10000)

#### set working directory ####

# here create new folder and set working directory within it

dir.create("~/Cel_GRN_revisions/")
setwd("~/Cel_GRN_revisions/")

# create subfolders for input, output and graphics

dir.create("input")

# into input folder, add input files 

dir.create("output")

dir.create("output/GRNs")

dir.create("graphics")

#### DEFINE FUNCTIONS ####

#### LOAD PACKAGES & FUNCTIONS ####

## First specify the packages of interest

packages <- c("biomaRt",
              "openxlsx",
              "FSA",
              "ggplot2",
              "DescTools",
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

biocmanager_packages <- c("TxDb.Celegans.UCSC.ce11.refGene",
                          "DESeq2",
                          "GenomicRanges") 

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

if(!file.exists("input/aam8940_cao_sm_tables_s1_to_s14.xlsx")){
  
download.file(url = "https://www.science.org/doi/suppl/10.1126/science.aam8940/suppl_file/aam8940_cao_sm_tables_s1_to_s14.xlsx",
              destfile = "input/aam8940_cao_sm_tables_s1_to_s14.xlsx")
  
}

Cao_celltypes <- read.xlsx("input/aam8940_cao_sm_tables_s1_to_s14.xlsx",
                           sheet = 4
)

colnames(Cao_celltypes) <- Cao_celltypes[1, ]
Cao_celltypes <- Cao_celltypes[2:nrow(Cao_celltypes), ]

row.names(Cao_celltypes) <- Cao_celltypes[, 1]
Cao_celltypes <- Cao_celltypes[, 3:ncol(Cao_celltypes)]

Cao_celltypes[] <- lapply(Cao_celltypes, as.numeric)

Cao_tissuetypes <- read.xlsx("input/aam8940_cao_sm_tables_s1_to_s14.xlsx",
                             sheet = 3
)

colnames(Cao_tissuetypes) <- Cao_tissuetypes[1, ]
Cao_tissuetypes <- Cao_tissuetypes[2:nrow(Cao_tissuetypes), ]

row.names(Cao_tissuetypes) <- Cao_tissuetypes[, 1]
Cao_tissuetypes <- Cao_tissuetypes[, 3:ncol(Cao_tissuetypes)]

Cao_tissuetypes[] <- lapply(Cao_tissuetypes, as.numeric)

Cao_tissuespecific <- read.xlsx("~/Cel_GRN_revisions/input/aam8940_cao_sm_tables_s1_to_s14.xlsx",
                                sheet = 7
)

colnames(Cao_tissuespecific) <- Cao_tissuespecific[1, ]
Cao_tissuespecific <- Cao_tissuespecific[2:nrow(Cao_tissuespecific), ]

CelEsT <- read.table("output/CelEsT_annotated.txt",
                     header = TRUE)

# annotate CelEsT
motif1000 <- read.table("~/Cel_GRN_manuscript/output/GRNs/FIMO_nohomo_1000.txt",
                        header = TRUE)

ChIP_nocut <- read.table("~/Cel_GRN_manuscript/output/GRNs/allChIP_10000_HOTexcl.txt",
                         header = TRUE)

eY1H <- read.table("~/Cel_GRN_manuscript/output/GRNs/walhoutGRN_unfiltered.txt",
                   header = TRUE)

CelEsT[, "in_motif"] <- paste0(CelEsT$source, "_",CelEsT$target) %in% paste0(motif1000$source, "_", motif1000$target)
CelEsT[!CelEsT$source %in% unique(motif1000$source), "in_motif"] <- NA

CelEsT[, "in_ChIP"] <- paste0(CelEsT$source, "_",CelEsT$target) %in% paste0(ChIP_nocut$source, "_", ChIP_nocut$target)
CelEsT[!CelEsT$source %in% unique(ChIP_nocut$source), "in_ChIP"] <- NA  

CelEsT[, "in_eY1H"] <- paste0(CelEsT$source, "_",CelEsT$target) %in% paste0(eY1H$source, "_", eY1H$target)
CelEsT[!CelEsT$source %in% unique(eY1H$source), "in_eY1H"] <- NA  

CelEsT_targets <- unique(CelEsT$target)

# note Captcha filter might prevent automated download.

if(!file.exists("input/pnas.1608162113.sd02.xlsx")){ 
  
  download.file(url = "http://www.pnas.org/lookup/suppl/doi:10.1073/pnas.1608162113/-/DCSupplemental/pnas.1608162113.sd02.xlsx",
                destfile = "input/pnas.1608162113.sd02.xlsx")
  
}

Ahringer_domains_L3_active <- read.xlsx("input/pnas.1608162113.sd02.xlsx",
                                        sheet = 4)

Ahringer_domains_L3_active <- rbind(colnames(Ahringer_domains_L3_active),
                                    Ahringer_domains_L3_active)
colnames(Ahringer_domains_L3_active) <- c("chr", "start", "end")

Ahringer_domains_L3_border <- read.xlsx("input/pnas.1608162113.sd02.xlsx",
                                        sheet = 5)

Ahringer_domains_L3_border <- rbind(colnames(Ahringer_domains_L3_border),
                                    Ahringer_domains_L3_border)
colnames(Ahringer_domains_L3_border) <- c("chr", "start", "end")

Ahringer_domains_L3_regulated <- read.xlsx("input/pnas.1608162113.sd02.xlsx",
                                           sheet = 6)

Ahringer_domains_L3_regulated <- rbind(colnames(Ahringer_domains_L3_regulated),
                                       Ahringer_domains_L3_regulated)
colnames(Ahringer_domains_L3_regulated) <- c("chr", "start", "end")

# set biomaRt site
parasite_mart <- useMart("parasite_mart", dataset = "wbps_gene", host = "https://parasite.wormbase.org", port = 443)

#### GINI COEFFICIENTS ####

allGenes_celltype_Gini <- apply(Cao_celltypes[, ], 1, DescTools::Gini)

allGenes_tissuetype_Gini <- apply(Cao_tissuetypes[, ], 1, DescTools::Gini)

motif_targets <- CelEsT[CelEsT$in_motif, "target"]
motif_targets <- motif_targets[!is.na(motif_targets)]

ChIP_targets <- CelEsT[CelEsT$in_ChIP, "target"]
ChIP_targets <- ChIP_targets[!is.na(ChIP_targets)]

eY1H_targets <- CelEsT[CelEsT$in_eY1H, "target"]
eY1H_targets <- eY1H_targets[!is.na(eY1H_targets)]

gini_df <- data.frame(gini = c(allGenes_tissuetype_Gini[wormRef::Cel_genes[match(ChIP_targets, Cel_genes$sequence_name), "wb_id"]],
                               allGenes_tissuetype_Gini[wormRef::Cel_genes[match(motif_targets, Cel_genes$sequence_name), "wb_id"]],
                               allGenes_tissuetype_Gini[wormRef::Cel_genes[match(eY1H_targets, Cel_genes$sequence_name), "wb_id"]],
                               allGenes_tissuetype_Gini[wormRef::Cel_genes[match(CelEsT_targets, Cel_genes$sequence_name), "wb_id"]])
                      ,
                      group = c(rep("ChIP", times = length(ChIP_targets)),
                                rep("motif", times = length(motif_targets)),
                                rep("eY1H", times = length(eY1H_targets)),
                                rep("unique", times = length(CelEsT_targets))))

gini_df$group <- factor(gini_df$group, levels = c("unique", "ChIP", "motif"))

saveRDS(gini_df, 
        "plotdata/ChIP_motif_gini_df.txt")

pdf("graphics/Tissue_gini_violin.pdf",
    height = 2.5,
    width = 2.5)

ggplot(gini_df[gini_df$group %in% c("ChIP", "motif", "unique"), ], aes(x = group, y = gini, fill = group)) + 
  geom_violin() +
  theme_classic() + 
  scale_fill_manual(values = c("grey", "orange", "magenta")) +
  ylab("Cross-tissue Gini coefficient") + 
  xlab("Data type") +
  theme(axis.text = element_text(colour = "black"),
        legend.position = "none")

dev.off()

dunnresults <- dunnTest(gini ~ group,
                        data = gini_df[gini_df$group %in% c("ChIP", "motif", "all"), ]
)

dunnresults$res

wilcox.test(gini_df[gini_df$group %in% c("ChIP"), "gini"], gini_df[gini_df$group %in% c("motif"), "gini"])

#### Ahringer chromatin domains ####

domains_vec <- c(rep("active", times = nrow(Ahringer_domains_L3_active)),
                 rep("border", times = nrow(Ahringer_domains_L3_border)),
                 rep("regulated", times = nrow(Ahringer_domains_L3_regulated)))

coords <- rbind(Ahringer_domains_L3_active,
                Ahringer_domains_L3_border[, 1:3],
                Ahringer_domains_L3_regulated)

domain_df <- cbind(coords, domains_vec)

domain_df[, 1] <- paste0("chr", domain_df[, 1])

domain_GR <- GenomicRanges::makeGRangesFromDataFrame(domain_df,
                                                     keep.extra.columns = TRUE)

# add 1 for 1-based coordinates (Ahringer supplementary in 0-based / BED format)
start(domain_GR) <- start(domain_GR) + 1

domain_GR <- sortSeqlevels(domain_GR)
domain_GR <- sort(domain_GR)

genes <- genes(TxDb.Celegans.UCSC.ce11.refGene)

gene_domain_overlap <- GenomicRanges::findOverlaps(genes, domain_GR, ignore.strand = TRUE)
domains_hit <- domain_GR$domains_vec[subjectHits(gene_domain_overlap)]
genes_hit <- genes$gene_id[queryHits(gene_domain_overlap)]

gene_domains_df <- data.frame(gene = genes_hit,
                              domain = domains_hit)

genes_convert <- getBM(mart = parasite_mart,
                       filters = "entrezgene_id",
                       value = gene_domains_df$gene,
                       attributes = c("wormbase_gene", "wormbase_gseq", "description", "wormbase_locus", "entrezgene_id"))

gene_domains_df[, c("wormbase_gseq", "wormbase_locus")] <- genes_convert[match(gene_domains_df$gene, genes_convert$entrezgene_id), c("wormbase_gseq", "wormbase_locus")]

domainplot_df <- data.frame(domain = c(gene_domains_df[match(ChIP_targets, gene_domains_df$wormbase_gseq), "domain"],
                                       gene_domains_df[match(motif_targets, gene_domains_df$wormbase_gseq), "domain"],
                                       gene_domains_df[match(eY1H_targets, gene_domains_df$wormbase_gseq), "domain"],
                                       gene_domains_df[match(CelEsT_targets, gene_domains_df$wormbase_gseq), "domain"]),
                            group = c(rep("ChIP", times = length(ChIP_targets)),
                                      rep("motif", times = length(motif_targets)),
                                      rep("eY1H", times = length(eY1H_targets)),
                                      rep("all", times = length(CelEsT_targets))),
                            genename = c(ChIP_targets, motif_targets, eY1H_targets, CelEsT_targets))

domain_freq_df <- sapply(unique(domainplot_df$group), function(x){table(domainplot_df[domainplot_df$group == x, "domain"])})

domain_frac_df <- t(t(domain_freq_df) / colSums(domain_freq_df))

# limit to ChIP and motif (and active and regulated) for plot
domain_frac_plot_df <- domain_frac_df[c("active", "regulated"), c("ChIP", "motif", "all")]

domain_frac_plot_long <- reshape2::melt(domain_frac_plot_df)

domain_frac_plot_long$Var2 <- factor(domain_frac_plot_long$Var2, levels = c("all", "ChIP", "motif")) 

saveRDS(domain_frac_plot_long, 
        "plotdata/ChIP_motif_chomdomain_fractions.txt")

pdf("graphics/Chromatin_domains_nolegend.pdf",
    height = 2.5,
    width = 2.5)

ggplot(data = domain_frac_plot_long, aes(x = Var1, y = value, fill = Var2)) + 
  geom_bar(position="dodge", stat="identity", colour = "black") + 
  theme_classic() + 
  theme(axis.text = element_text(colour = "black"),
        legend.position = "none"
  ) + 
  ylab("Fraction of target genes") +
  xlab("Chromatin domain type") + 
  coord_cartesian(ylim = c(0, 1)) + 
  scale_fill_manual(values = c("grey", "orange", "magenta"))

dev.off()

domain_chisq_test_list <- lapply(1:3, function(x){chisq.test(domain_freq_df[, c(x, 4)])})
names(domain_chisq_test_list) <- colnames(domain_freq_df)[1:3]

#### tissue-specific expression ####

Cao_tissuespecific$qval <- as.numeric(Cao_tissuespecific$qval)

# as in Meeuse 2020 will use ratio >5 and p < 0.05 to define tissue specific genes

Cao_tissuespecific <- Cao_tissuespecific[Cao_tissuespecific$ratio > 5 & Cao_tissuespecific$qval < 0.05, ]

Celltypes <- names(table(Cao_tissuespecific$max.cell.type))

## pool by larger tissue type

Celltypes_tissuecounterpart <- c("Neurons",
                                 "Body_wall_muscle",
                                 rep("Neurons", 3),
                                 "none",
                                 "Gonad",
                                 "none",
                                 rep("Neurons",2),
                                 "Gonad",
                                 "none",
                                 "Intestine",
                                 "Hypodermis",
                                 rep("Neurons", 2),
                                 rep("Pharynx",3),
                                 "Neurons",
                                 "none",
                                 "Hypodermis",
                                 "none",
                                 "Glia",
                                 "Gonad",
                                 "Neurons")

names(Celltypes_tissuecounterpart) <- Celltypes

# number of tissue-specific genes
tissue_ns <- sapply(unique(Celltypes_tissuecounterpart), function(x){length(Cel_genes[match(Cao_tissuespecific[Cao_tissuespecific$max.cell.type %in% names(Celltypes_tissuecounterpart[Celltypes_tissuecounterpart == x]), "gene_id"], Cel_genes$wb_id), "sequence_name"])})

tissue_ns <- tissue_ns[names(tissue_ns) != "none"]

ChIP_fisher <- lapply(unique(Celltypes_tissuecounterpart), function(x){

  temp_ct <- rbind(
    c(sum(ChIP_targets %in% Cel_genes[match(Cao_tissuespecific[Cao_tissuespecific$max.cell.type %in% names(Celltypes_tissuecounterpart[Celltypes_tissuecounterpart == x]), "gene_id"], Cel_genes$wb_id), "sequence_name" ]), 
      sum(!ChIP_targets %in% Cel_genes[match(Cao_tissuespecific[Cao_tissuespecific$max.cell.type %in% names(Celltypes_tissuecounterpart[Celltypes_tissuecounterpart == x]), "gene_id"], Cel_genes$wb_id), "sequence_name" ]))
    ,
    c(sum(unique(ChIP_targets) %in% Cel_genes[match(Cao_tissuespecific[Cao_tissuespecific$max.cell.type %in% names(Celltypes_tissuecounterpart[Celltypes_tissuecounterpart == x]), "gene_id"], Cel_genes$wb_id), "sequence_name" ]),
      sum(!unique(ChIP_targets) %in% Cel_genes[match(Cao_tissuespecific[Cao_tissuespecific$max.cell.type %in% names(Celltypes_tissuecounterpart[Celltypes_tissuecounterpart == x]), "gene_id"], Cel_genes$wb_id), "sequence_name" ]))
  )
  
  fisher.test(temp_ct)
  
})

names(ChIP_fisher) <- unique(Celltypes_tissuecounterpart)

motif_fisher <- lapply(unique(Celltypes_tissuecounterpart), function(x){
  
  temp_ct <- rbind(
    c(sum(motif_targets %in% Cel_genes[match(Cao_tissuespecific[Cao_tissuespecific$max.cell.type %in% names(Celltypes_tissuecounterpart[Celltypes_tissuecounterpart == x]), "gene_id"], Cel_genes$wb_id), "sequence_name" ]), 
      sum(!motif_targets %in% Cel_genes[match(Cao_tissuespecific[Cao_tissuespecific$max.cell.type %in% names(Celltypes_tissuecounterpart[Celltypes_tissuecounterpart == x]), "gene_id"], Cel_genes$wb_id), "sequence_name" ]))
    ,
    c(sum(unique(motif_targets) %in% Cel_genes[match(Cao_tissuespecific[Cao_tissuespecific$max.cell.type %in% names(Celltypes_tissuecounterpart[Celltypes_tissuecounterpart == x]), "gene_id"], Cel_genes$wb_id), "sequence_name" ]),
      sum(!unique(motif_targets) %in% Cel_genes[match(Cao_tissuespecific[Cao_tissuespecific$max.cell.type %in% names(Celltypes_tissuecounterpart[Celltypes_tissuecounterpart == x]), "gene_id"], Cel_genes$wb_id), "sequence_name" ]))
  )
  
  fisher.test(temp_ct)
  
})

names(motif_fisher) <- unique(Celltypes_tissuecounterpart)

eY1H_fisher <- lapply(unique(Celltypes_tissuecounterpart), function(x){
  
  temp_ct <- rbind(
    c(sum(eY1H_targets %in% Cel_genes[match(Cao_tissuespecific[Cao_tissuespecific$max.cell.type %in% names(Celltypes_tissuecounterpart[Celltypes_tissuecounterpart == x]), "gene_id"], Cel_genes$wb_id), "sequence_name" ]), 
      sum(!eY1H_targets %in% Cel_genes[match(Cao_tissuespecific[Cao_tissuespecific$max.cell.type %in% names(Celltypes_tissuecounterpart[Celltypes_tissuecounterpart == x]), "gene_id"], Cel_genes$wb_id), "sequence_name" ]))
    ,
    c(sum(unique(eY1H_targets) %in% Cel_genes[match(Cao_tissuespecific[Cao_tissuespecific$max.cell.type %in% names(Celltypes_tissuecounterpart[Celltypes_tissuecounterpart == x]), "gene_id"], Cel_genes$wb_id), "sequence_name" ]),
      sum(!unique(eY1H_targets) %in% Cel_genes[match(Cao_tissuespecific[Cao_tissuespecific$max.cell.type %in% names(Celltypes_tissuecounterpart[Celltypes_tissuecounterpart == x]), "gene_id"], Cel_genes$wb_id), "sequence_name" ]))
  )
  
  fisher.test(temp_ct)
  
})

names(eY1H_fisher) <- unique(Celltypes_tissuecounterpart)

odds_ratio_table <- data.frame(ChiP = sapply(ChIP_fisher, function(x){x$estimate}),
                               motif = sapply(motif_fisher, function(x){x$estimate}))

odds_ratio_table <- odds_ratio_table[!str_detect(row.names(odds_ratio_table), "none"), ]

OR_melt <- reshape2::melt(odds_ratio_table)

pval_table <- data.frame(ChiP = sapply(ChIP_fisher, function(x){x$p.value}),
                         motif = sapply(motif_fisher, function(x){x$p.value}))

pval_table <- pval_table[!str_detect(row.names(pval_table), "none"), ]

confint_min_table <- data.frame(ChiP = sapply(ChIP_fisher, function(x){x$conf.int[[1]]}),
                                motif = sapply(motif_fisher, function(x){x$conf.int[[1]]}))

confint_min_table <- confint_min_table[!str_detect(row.names(confint_min_table), "none"), ]
min_melt <- reshape2::melt(confint_min_table)

confint_max_table <- data.frame(ChiP = sapply(ChIP_fisher, function(x){x$conf.int[[2]]}),
                                motif = sapply(motif_fisher, function(x){x$conf.int[[2]]}))

confint_max_table <- confint_max_table[!str_detect(row.names(confint_max_table), "none"), ]
max_melt <- reshape2::melt(confint_max_table)

pval_table <- pval_table[!str_detect(row.names(pval_table), "none"), ]

pval_melt <- reshape2::melt(pval_table)

OR_melt[, "pval"] <- pval_melt[, 2]
OR_melt[, "adjpval"] <- p.adjust(OR_melt$pval, method = "BH")

OR_melt[, "minuslog10pval"] <- -log10(OR_melt$adjpval)

OR_melt[, "tissue"] <- rep(row.names(pval_table), 2)
OR_melt[, "significant"] <- OR_melt$adjpval < 0.1

OR_melt[, "pointcol"] <- "grey"
OR_melt[OR_melt$significant, "pointcol"] <- "red"

OR_melt[, "upper"] <- max_melt$value
OR_melt[, "lower"] <- min_melt$value

OR_melt[, "alpha"] <- 0.5
OR_melt[OR_melt$significant, "alpha"] <- 1

OR_melt$tissue <- paste0(OR_melt$tissue, " (", tissue_ns, ")")

OR_melt$tissue <- str_replace_all(OR_melt$tissue, "Body_wall_muscle", "BW muscle")

saveRDS(OR_melt, 
        "plotdata/ChIP_motif_tissuespecific_genes.txt")

# for a plot with different coloured facet labels.. its complicated. 
# to save time, just make one orange and one magenta and combine as appropriate in illustrator

pdf("graphics/ChIP_motif_tissuespeific_ORplot_orangelabel.pdf",
    width = 6,
    height = 3.5)

ggplot(OR_melt, aes(x = tissue, y = value)) + 
  geom_point(aes(size = minuslog10pval), col = OR_melt$pointcol) + 
  ylab("Odds ratio") + 
  xlab("Tissue") + 
  theme_classic() +
  geom_errorbar(aes(ymin = lower, ymax = upper)) +
  coord_cartesian(ylim = c(0, 2)) + 
  geom_hline(yintercept = 1) + 
  facet_wrap(~ variable) + 
  theme(axis.text.x = element_text(angle = 45,
                                   vjust = 1,
                                   hjust = 1,
                                   colour = "black"),
        axis.text.y = element_text(colour = "black"),
        strip.text = element_text(size = 12, color = "black", face = "bold"),
        strip.background = element_rect(fill = c("orange"))
        # , legend.position = "none"
  )

dev.off()

pdf("graphics/ChIP_motif_tissuespeific_ORplot_legend.pdf",
    width = 6,
    height = 3.5)

ggplot(OR_melt, aes(x = tissue, y = value)) + 
  geom_point(aes(size = minuslog10pval/4), col = OR_melt$pointcol) + 
  ylab("Odds ratio") + 
  xlab("Tissue") + 
  theme_classic() +
  geom_errorbar(aes(ymin = lower, ymax = upper)) +
  coord_cartesian(ylim = c(0, 2)) + 
  geom_hline(yintercept = 1) + 
  facet_wrap(~ variable) + 
  theme(axis.text.x = element_text(angle = 45,
                                   vjust = 1,
                                   hjust = 1,
                                   colour = "black"),
        strip.text = element_text(size = 12, color = "white", face = "bold"),
        strip.background = element_rect(fill = "black")) + 
  scale_size_continuous(range  = c(0.1, 10), 
                        limits = c(0, 7), 
                        breaks = c(0, 0.75, 1.5, 2.25))

dev.off()

# so ChIP favours genes in active chromatin regions, expressed either ubiquitously or in tissues with a lot of chromatin.


