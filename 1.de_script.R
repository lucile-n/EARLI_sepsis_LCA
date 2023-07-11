######################################################
# 1.de_script.R
# created on Oct 22 2020
# lucile.neyton@ucsf.edu
######################################################

rm(list = ls())
setwd("/Users/lucileneyton/Box Sync/EARLI_VALID/de_analysis_sepsis/")

# load libraries
library(ggplot2)
library(edgeR)
library(DESeq2)
library(ggfortify)
library(biomaRt)
library(ggrepel)
library(plyr)
library(GGally)
library(corrplot)

# set ggplot2 theme
theme_update(text = element_text(family = "Helvetica", size=8),
             axis.text.x = element_text(family = "Helvetica", size=8),
             axis.text.y = element_text(family = "Helvetica", size=8),
             axis.title.x = element_text(family = "Helvetica", size=8),
             axis.title.y = element_text(family = "Helvetica", size=8),
             legend.text = element_text(family = "Helvetica", size=8),
             legend.title = element_text(family = "Helvetica", size=8),
             strip.text.x = element_text(family = "Helvetica", size=8),
             strip.text.y = element_text(family = "Helvetica", size=8),
             plot.title = element_text(family = "Helvetica", size=8))
fig_width <- 6
fig_height <- 6

# set paths
data_path <- "/Users/lucileneyton/Box Sync/EARLI_VALID/data/"
results_path <-
  "/Users/lucileneyton/Box Sync/EARLI_VALID/de_analysis_sepsis/results/"

#########################
# DATA LOADING
#########################
# list data files
# downloaded on Oct 5 2020
cnt_data_path <- paste(data_path, "raw/earli_counts_kallisto.csv", sep = "")
# LCA labels file
meta_data_path <- paste(data_path, "processed/covs_df_sepsis.csv", sep = "")

# read data in
cnt_data <- read.csv(cnt_data_path, row.names = 1)
meta_data <- read.csv(meta_data_path, row.names = 1)

#########################
# DATA PREPROCESSING
#########################
# make sure our variable of interest is a factor
meta_data$LCA_label <- as.factor(meta_data$LCA_label)

# drop rows of meta data pertaining to samples without sequencing data
rownames(meta_data) <- paste("EARLI", meta_data$Barcode, sep="_")
meta_data <- meta_data[rownames(meta_data) %in% colnames(cnt_data), ]

# filter count data to exclude samples without a label
cnt_data <- cnt_data[, rownames(meta_data)]

# 5391 is male sex at birth
meta_data[meta_data$X==5391, "Gender1"] <- "Male"

# save the output for future analyses
write.csv(meta_data, paste(data_path, "processed/covs_df_sepsis_cleaned.csv",
                           sep = ""))
write.csv(cnt_data, paste(data_path, "processed/cnt_data.csv",
                           sep = ""))

# standard scaling on Age
meta_data$agephi <- scale(meta_data$agephi)

#########################
# DE ANALYSIS
#########################
# cpm for cibersort
cpm_data <- cpm(cnt_data)

# save VST data
write.csv(t(cpm_data),
          paste(data_path, "processed/cpm_data.csv", sep = ""))

# build DESeq2 object
dds <- DESeqDataSetFromMatrix(countData = cnt_data, colData = meta_data,
                              design = ~ LCA_label + agephi + Gender1)

# choose the reference level for the factor of interest
dds$LCA_label <- relevel(dds$LCA_label, ref = "Hypo")

# run DESeq
dds <- DESeq(dds)

# transform data for corr plots
vsd <- vst(dds, blind = TRUE)

# extract results
res <- results(dds, contrast = c("LCA_label", "Hyper", "Hypo"))

# add gene symbols
ensembl <- useEnsembl(
  biomart = "ensembl", dataset = "hsapiens_gene_ensembl",
  version = 99
) # version 99

ensembl_res <- getBM(
  values = rownames(res),
  filters = "ensembl_gene_id",
  attributes = c("ensembl_gene_id", "hgnc_symbol"),
  mart = ensembl
)

# drop duplicates -> we only want one match per ENSG
ensembl_res <- ensembl_res[!duplicated(ensembl_res$ensembl_gene_id), ]

# save mapping
write.csv(ensembl_res, paste0(data_path, "processed/genes_mapping2.csv"), row.names = F)

# add gene symbols and descriptions to DESeq's results
res$hgnc_symbol <- ensembl_res$hgnc_symbol[match(rownames(res),
                                                 ensembl_res$ensembl_gene_id)]

# sort the genes from lowest to highest given adjusted p-values
res <- res[order(res$padj, decreasing = F), ]

# replace NA values with 1s
res$padj[is.na(res$padj)] <- 1
sig_results <- data.frame(res[res$padj < 0.05, ])

# save the output as a CSV file
write.csv(sig_results, paste(results_path, "DGEA_results_sepsis.csv", sep = ""))
write.csv(data.frame(res), paste(results_path, "DGEA_results_all_sepsis.csv", sep = ""))

# generate a volcano plot
# only display top 25 gene symbols
res_df <- data.frame(res)
res_df$sig <- res_df$padj < 0.05

gene_symbols_italic <- res_df[1:25, "hgnc_symbol"]
gene_symbols_italic <- sapply(gene_symbols_italic, function(x) paste0("italic('", paste0(x, "')")))

pdf(paste(results_path, "volcano_plot.pdf", sep = ""), width = 0.5*fig_width, height = 0.5*fig_height)
p <- ggplot(res_df, aes(log2FoldChange, -log10(pvalue))) +
  geom_point(aes(col = sig)) +
  scale_color_manual(values = c("black", "red")) +
  theme_bw() +
  theme(legend.position = "none", text = element_text(family = "Helvetica", size=8),
        axis.text.x = element_text(family = "Helvetica", size=8),
        axis.text.y = element_text(family = "Helvetica", size=8),
        axis.title.x = element_text(family = "Helvetica", size=8),
        axis.title.y = element_text(family = "Helvetica", size=8),
        legend.text = element_text(family = "Helvetica", size=8),
        legend.title = element_text(family = "Helvetica", size=8),
        strip.text.x = element_text(family = "Helvetica", size=8),
        strip.text.y = element_text(family = "Helvetica", size=8),
        plot.title = element_text(family = "Helvetica", size=8)) +
  geom_text_repel(data = res_df[1:25, ], size=8*0.35,
                  aes(label = gene_symbols_italic),
                  max.overlaps=20, parse=T) +
  xlab("Log2 Fold Change") +
  ylab("-Log10 P-value")
print(p)
dev.off()

#########################
# DE ANALYSIS - WITH ANC CORRECTION
#########################
# for ANC/WBC values measured on samples > 12 hours from PAXgene tube
# exclude samples with NA values
meta_data_filt <- meta_data[!is.na(meta_data$frac_anc), ]
cnt_data_filt <- cnt_data[, rownames(meta_data_filt)]

# build DESeq2 object
dds_anc <- DESeqDataSetFromMatrix(countData = cnt_data_filt, colData = meta_data_filt,
                              design = ~ LCA_label + agephi + Gender1 + frac_anc)

# choose the reference level for the factor of interest
dds_anc$LCA_label <- relevel(dds_anc$LCA_label, ref = "Hypo")

# run DESeq
dds_anc <- DESeq(dds_anc)

# extract results
res_anc <- results(dds_anc, contrast = c("LCA_label", "Hyper", "Hypo"))

# add gene symbols and descritptions to DESeq's results
res_anc$hgnc_symbol <- ensembl_res$hgnc_symbol[match(rownames(res_anc),
                                                     ensembl_res$ensembl_gene_id)]

# replace NA values with 1s
res_anc$padj[is.na(res_anc$padj)] <- 1

# sort the genes from lowest to highest given adjusted p-values
res_anc <- res_anc[order(res_anc$padj, decreasing = F), ]

sig_results_anc <- data.frame(res_anc[res_anc$padj < 0.05, ])

# compare with results not including ANC values as covariate
# 0.92 correlation
cor.test(res$stat, res_anc[rownames(res), "stat"], method = "spearman")

#########################
# BIOMARKERS VS. GENE EXPRESSSION
#########################
# plasma biomarkers of interest
covs_for_comp <- c("bilirubin_fin", "PAI1", "IL6", "sTNFr1", "RAGE", 
                   "IL8", "ANG2", "Proc_C", "ICAM", "naminsaps", "hco3minsaps", 
                   "creatininemaxsaps", "albuminminsaps", "glucoseminsaps", 
                   "wbcmaxsaps", "hctminsaps", "plateletsminsaps")
meta_data_for_comp <- meta_data[, covs_for_comp]

# rename variables for vis
from_ <- covs_for_comp
to_ <- c("Bilirubin", "PAI-1", "Interleukin-6", "sTNFr1", "RAGE", 
         "Interleukin-8", "ANG-2", "Protein-C", "ICAM", "Sodium", "Bicarbonate",
         "Creatinine", "Albumin", "Glucose", "White cell count", "Hematocrit",
         "Platelets")

colnames(meta_data_for_comp) <- mapvalues(colnames(meta_data_for_comp), 
                                          from_, to_)

covs_for_comp <- colnames(meta_data_for_comp)

# add LCA
meta_data_for_comp$LCA_label <- meta_data$LCA_label

vsd_assay <- assay(vsd)
rownames(vsd_assay)[rownames(vsd_assay)=="ENSG00000242574"] <- "HLA-DMB"
rownames(vsd_assay)[rownames(vsd_assay)=="ENSG00000169896"] <- "ITGAM"
rownames(vsd_assay)[rownames(vsd_assay)=="ENSG00000100985"] <- "MMP9"
rownames(vsd_assay)[rownames(vsd_assay)=="ENSG00000136634"] <- "IL10"

# test
meta_data_for_comp[, c("Albumin", "Interleukin-6", "Interleukin-8")] <- log10(meta_data_for_comp[, c("Albumin", "Interleukin-6", "Interleukin-8")] +1)

# glm to get coefficients
for (biom_ in colnames(meta_data_for_comp[, c("Albumin", "Interleukin-6", "Interleukin-8")])){
  for (gene_ in colnames(t(vsd_assay[c("HLA-DMB", "ITGAM", "MMP9", "IL10") , ]))){
    glm_data <- cbind(meta_data_for_comp[, c(biom_, "LCA_label")], t(vsd_assay[c("HLA-DMB", "ITGAM", "MMP9", "IL10") , ])[, gene_])
    colnames(glm_data)[3] <- gene_
    
    print(biom_)
    print(gene_)
  
    # for all
    print(summary(glm(as.formula(paste0(paste0(paste0("`", gene_), "`~`"), paste0(biom_,"`"))), data = glm_data))$coefficients[2, c(1, 4)])
    
    # hypo
    print(summary(glm(as.formula(paste0(paste0(paste0("`", gene_), "`~`"), paste0(biom_,"`"))), data = glm_data[glm_data$LCA_label=="Hypo", ]))$coefficients[2, c(1, 4)])
    # hyper
    print(summary(glm(as.formula(paste0(paste0(paste0("`", gene_), "`~`"), paste0(biom_,"`"))), data = glm_data[glm_data$LCA_label=="Hyper",]))$coefficients[2, c(1, 4)])
    
    # interaction
    print(summary(glm(as.formula(paste0(paste0(paste0("`", gene_), "`~LCA_label*`"), paste0(biom_,"`"))), data = glm_data))$coefficients[4, c(1, 4)])
    
  }
}


ggpairs_plot <- ggpairs(cbind(meta_data_for_comp[, c("Albumin", "Interleukin-6", "Interleukin-8")], t(vsd_assay[c("HLA-DMB", "ITGAM", "MMP9", "IL10") , ])), 
        aes(color=meta_data_for_comp$LCA_label),
        lower = list(continuous = "blank"), 
        upper = list(continuous = wrap("smooth", size=1)), 
        diag = list(continous = "blank"),
        axisLabels = "show") + 
  scale_color_manual(breaks=c("Hyper", "Hypo"), values = c("#EE7D31","#4472C4")) +
  theme_bw()

ggpairs_plot$plots <- ggpairs_plot$plots[1:21]

pdf(paste0(results_path, "ggpairs_plot.pdf"), width = 7, height = 7)
print(ggpairs_plot +
        theme(text = element_text(family = "Helvetica", size=8),
              axis.text.x = element_text(family = "Helvetica", size=8),
              axis.text.y = element_text(family = "Helvetica", size=8),
              axis.title.x = element_text(family = "Helvetica", size=8),
              axis.title.y = element_text(family = "Helvetica", size=8),
              legend.text = element_text(family = "Helvetica", size=8),
              legend.title = element_text(family = "Helvetica", size=8),
              strip.text.x = element_text(family = "Helvetica", size=8, face = "italic"),
              strip.text.y = element_text(family = "Helvetica", size=8),
              plot.title = element_text(family = "Helvetica", size=8)))
dev.off()

#########################
# GSEA - IPA
#########################
# read-in IPA data
load(paste0(results_path, "earli-hypervshypo-ipa.Rda"))
hypervshypo_cp <- earliphenotypes.ipa$`Canonical Pathways`

# filter out NAs for Z-scores
hypervshypo_cp <- hypervshypo_cp[!is.na(hypervshypo_cp$zScore), ]

# as per IPA's guidelines, 2 as activation score + 0.001 as a p-value threshold
hypervshypo_cp_sig <- hypervshypo_cp[abs(hypervshypo_cp$zScore)>2, ]

write.csv(hypervshypo_cp_sig, paste(results_path, "hypervshypo_cp_sig.csv", sep = ""))

# hyper vs hypo figure
top_n <- 10

cp_sel <- hypervshypo_cp_sig[order(abs(hypervshypo_cp_sig[, "zScore"]), decreasing = T), "Ingenuity Canonical Pathways"][1:top_n]

hypervshypo_cp_sig_full <- hypervshypo_cp_sig
hypervshypo_cp_sig <- hypervshypo_cp_sig[hypervshypo_cp_sig$`Ingenuity Canonical Pathways` %in% cp_sel, ]

hypervshypo_cp_sig <- hypervshypo_cp_sig[order(hypervshypo_cp_sig[, "zScore"], decreasing = T), c("Ingenuity Canonical Pathways", "geneNames", "zScore")]
hypervshypo_cp_sig$zscore_abs <- abs(hypervshypo_cp_sig$zScore)
hypervshypo_cp_sig$zscore_sign[hypervshypo_cp_sig$zScore > 0] <- "+"
hypervshypo_cp_sig$zscore_sign[hypervshypo_cp_sig$zScore < 0] <- "-"

hypervshypo_cp_sig$zscore_sign <- as.factor(hypervshypo_cp_sig$zscore_sign)
hypervshypo_cp_sig$zscore_sign <- relevel(hypervshypo_cp_sig$zscore_sign, "+")

hypervshypo_cp_sig$cp <- hypervshypo_cp_sig$`Ingenuity Canonical Pathways`
hypervshypo_cp_sig$cp <- factor(hypervshypo_cp_sig$cp, levels=hypervshypo_cp_sig$cp)

# dot plot
pdf(paste(results_path, "ipa_cp_plot.pdf", sep = ""), width = fig_width*0.7, height = fig_height*0.5)
print(ggplot(hypervshypo_cp_sig, aes(x=zScore, y=cp)) +
        geom_point(aes(color=zscore_sign), size=3, show.legend=F) +
        scale_color_manual(breaks=c("+", "-"), values = c("#EE7D31","#4472C4")) +
        theme_bw() +
        theme(text = element_text(family = "Helvetica", size=8),
              axis.text.x = element_text(family = "Helvetica", size=8),
              axis.text.y = element_text(family = "Helvetica", size=8),
              axis.title.x = element_text(family = "Helvetica", size=8),
              axis.title.y = element_text(family = "Helvetica", size=8),
              legend.text = element_text(family = "Helvetica", size=8),
              legend.title = element_text(family = "Helvetica", size=8),
              strip.text.x = element_text(family = "Helvetica", size=8),
              strip.text.y = element_text(family = "Helvetica", size=8),
              plot.title = element_text(family = "Helvetica", size=8)) +
        xlab("Z-score") +
        ylab("Ingenuity Canonical Pathways") +
        scale_y_discrete(limits = rev(levels(hypervshypo_cp_sig$cp))))
dev.off()

hypervshypo_cp_sig_full <- hypervshypo_cp_sig_full[order(hypervshypo_cp_sig_full[, "zScore"], decreasing = T), c("Ingenuity Canonical Pathways", "geneNames", "zScore")]
hypervshypo_cp_sig_full$zscore_abs <- abs(hypervshypo_cp_sig_full$zScore)
hypervshypo_cp_sig_full$zscore_sign[hypervshypo_cp_sig_full$zScore > 0] <- "+"
hypervshypo_cp_sig_full$zscore_sign[hypervshypo_cp_sig_full$zScore < 0] <- "-"

hypervshypo_cp_sig_full$zscore_sign <- as.factor(hypervshypo_cp_sig_full$zscore_sign)
hypervshypo_cp_sig_full$zscore_sign <- relevel(hypervshypo_cp_sig_full$zscore_sign, "+")

hypervshypo_cp_sig_full$cp <- hypervshypo_cp_sig_full$`Ingenuity Canonical Pathways`
hypervshypo_cp_sig_full$cp <- factor(hypervshypo_cp_sig_full$cp, levels=hypervshypo_cp_sig_full$cp)

# dot plot
pdf(paste(results_path, "ipa_cp_full_plot.pdf", sep = ""), width = fig_width*0.8, height = fig_height*0.7)
print(ggplot(hypervshypo_cp_sig_full, aes(x=zScore, y=cp)) +
        geom_point(aes(color=zscore_sign), size=3, show.legend=F) +
        scale_color_manual(breaks=c("+", "-"), values = c("#EE7D31","#4472C4")) +
        theme_bw() +
        theme(text = element_text(family = "Helvetica", size=8),
              axis.text.x = element_text(family = "Helvetica", size=8),
              axis.text.y = element_text(family = "Helvetica", size=8),
              axis.title.x = element_text(family = "Helvetica", size=8),
              axis.title.y = element_text(family = "Helvetica", size=8),
              legend.text = element_text(family = "Helvetica", size=8),
              legend.title = element_text(family = "Helvetica", size=8),
              strip.text.x = element_text(family = "Helvetica", size=8),
              strip.text.y = element_text(family = "Helvetica", size=8),
              plot.title = element_text(family = "Helvetica", size=8)) +
        xlab("Z-score") +
        ylab("Ingenuity Canonical Pathways") +
        scale_y_discrete(limits = rev(levels(hypervshypo_cp_sig_full$cp))))
dev.off()

#########################
# GSEA - REACTOME
#########################
library(fgsea)

# load gene set DB
m_dbs <- list(reactome=msigdbr(species="Homo sapiens", "C2", "CP:REACTOME"))

m_db <- m_dbs[["reactome"]]    

m_db_list <- lapply(unique(m_db$gs_name), function(x) unname(unlist(m_db[m_db$gs_name == x, "ensembl_gene"])))
names(m_db_list) <- unique(m_db$gs_name)

# sort genes given log fc values
top_table_sorted <- res_df[order(res_df[, "log2FoldChange"], decreasing = T), ]
gene_ranks <- top_table_sorted[, "log2FoldChange"]
names(gene_ranks) <- rownames(top_table_sorted)
gene_ranks <- gene_ranks[!is.na(gene_ranks)]

# GSEA
set.seed(123)
gsea_res <- fgsea(m_db_list, stats=gene_ranks, 
                  maxSize=500, eps=0)

# filter
gsea_res$padj[is.na(gsea_res$padj)] <- 1

gsea_res_df <- gsea_res[(gsea_res$padj<0.05), ]

gsea_res_df$pathway <- as.factor(gsea_res_df$pathway)

write.csv(gsea_res_df[, c("pathway", "padj", "NES")], paste0(results_path, "gsea_res_df_reactome.csv"))

# dot plot
gsea_res_df_pos <- gsea_res_df[gsea_res_df$NES > 0, ]
gsea_res_df_neg <- gsea_res_df[gsea_res_df$NES < 0, ]

gsea_res_df_pos <- gsea_res_df_pos[order(gsea_res_df_pos$padj, decreasing = F), ]
gsea_res_df_neg <- gsea_res_df_neg[order(gsea_res_df_neg$padj, decreasing = F), ]

gsea_res_df_pos_top10 <- as.character(gsea_res_df_pos[1:(min(nrow(gsea_res_df_pos), 5)), "pathway"]$pathway)
gsea_res_df_neg_top10 <- as.character(gsea_res_df_neg[1:(min(nrow(gsea_res_df_neg), 5)), "pathway"]$pathway)

gsea_res_df <- gsea_res_df[order(gsea_res_df$NES, decreasing=T), ]
gsea_res_df_filt <- gsea_res_df[gsea_res_df$pathway %in% c(gsea_res_df_pos_top10, gsea_res_df_neg_top10), ]

gsea_res_df_filt$sign[gsea_res_df_filt$NES > 0] <- "+"
gsea_res_df_filt$sign[gsea_res_df_filt$NES < 0] <- "-"

gsea_res_df_filt$sign <- factor(gsea_res_df_filt$sign, levels=c("+", "-"))

gsea_res_df_filt$pathway <- str_replace(gsea_res_df_filt$pathway, paste0(toupper("reactome"), "_"), "")
gsea_res_df_filt$pathway <- tools::toTitleCase(tolower(str_replace_all(gsea_res_df_filt$pathway, "_", " ")))

gsea_res_df_filt$pathway <- factor(gsea_res_df_filt$pathway, levels=gsea_res_df_filt$pathway)

# dot plot
pdf(paste(results_path, "reactome_cp_plot.pdf", sep = ""), width = fig_width*0.6, height = fig_height*0.5)
print(ggplot(gsea_res_df_filt, aes(x=NES, y=pathway)) +
        geom_point(aes(color=sign), size=3, show.legend=F) +
        scale_color_manual(breaks=c("+", "-"), values = c("#EE7D31","#4472C4")) +
        theme_bw() +
        theme(text = element_text(family = "Helvetica", size=8),
              axis.text.x = element_text(family = "Helvetica", size=8),
              axis.text.y = element_text(family = "Helvetica", size=8),
              axis.title.x = element_text(family = "Helvetica", size=8),
              axis.title.y = element_text(family = "Helvetica", size=8),
              legend.text = element_text(family = "Helvetica", size=8),
              legend.title = element_text(family = "Helvetica", size=8),
              strip.text.x = element_text(family = "Helvetica", size=8),
              strip.text.y = element_text(family = "Helvetica", size=8),
              plot.title = element_text(family = "Helvetica", size=8)) +
        xlab("NES") +
        ylab("Reactome Pathways") +
        scale_y_discrete(limits = rev(levels(gsea_res_df_filt$pathway))))
dev.off()

#########################
# BIOMARKERS - IPA
#########################
# colors
pal_ <- colorRampPalette(rev(c("#EE7D31",
                               "#FFFFFF",
                               "#4472C4")))(200)

# for genes without symbol, use ENSG id
ensembl_res$hgnc_symbol[ensembl_res$hgnc_symbol==""] <- ensembl_res$ensembl_gene_id[ensembl_res$hgnc_symbol==""]

# for pathways of interest
poi_names <- hypervshypo_cp_sig_full$`Ingenuity Canonical Pathways`

clin_biom <- c("Bilirubin", "Sodium", "Bicarbonate", "Creatinine", "Albumin", "Glucose", "White cell count", "Hematocrit", "Platelets" )
research_biom <- c("PAI-1", "Interleukin-6", "sTNFr1", "RAGE", "Interleukin-8", "ANG-2", "Protein-C", "ICAM")

for (poi_name in poi_names){
  poi_genes <- strsplit(hypervshypo_cp_sig_full[hypervshypo_cp_sig_full$`Ingenuity Canonical Pathways`==poi_name, "geneNames"], ",")[[1]]
  poi_genes <- unique(poi_genes)
  poi_genes <- poi_genes[poi_genes!=""]
  
  poi_vsd <- vsd_assay[rownames(vsd_assay) %in% ensembl_res$ensembl_gene_id[match(poi_genes, ensembl_res$hgnc_symbol)], ]
  rownames(poi_vsd) <- ensembl_res$hgnc_symbol[match(rownames(poi_vsd),
                                                     ensembl_res$ensembl_gene_id)]
  
  if (dim(poi_vsd)[1] > 25){
    sig_results_tmp <- sig_results[sig_results$hgnc_symbol %in% rownames(poi_vsd), ]
    sig_genes_tmp <- sig_results_tmp[order(sig_results_tmp$padj), "hgnc_symbol"][1:25]
    
    poi_vsd <- poi_vsd[sig_genes_tmp, ]
  }
  
  if (dim(poi_vsd)[1] > 1){
    # correlation matrices
    poi_cov <- data.frame(row.names = rownames(poi_vsd))
    for (cov_ in colnames(meta_data_for_comp)){
      if (cov_!="LCA_label"){
        cov_data <- meta_data_for_comp[, cov_, drop=F]
        poi_cov[, cov_] <- apply(poi_vsd, 1, function(x) cor(x, cov_data, method = "spearman", use = "complete.obs"))
      }
      
    }

    # pathway genes vs biomarkers
    poi_name_output <- gsub(" ", "_", poi_name)
    
    pdf(paste(results_path, paste(poi_name_output, "_ipa_corrplot.pdf", sep=""), sep = ""), 
        width = 0.8*fig_width, height = 0.8*fig_height, fonts = "Helvetica", pointsize = 8)   
    par(mfrow=c(1, 2))
    print(corrplot(as.matrix(poi_cov[, clin_biom]), col=pal_, mar=c(1.2,1,0.9,1) + 0.1))
    print(corrplot(as.matrix(poi_cov[, research_biom]), col=pal_, mar = c(0.3,1,1.1,1) + 0.9))
    dev.off()
  }
}






