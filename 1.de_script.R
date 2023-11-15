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
library(WGCNA)

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
  ylab("-Log10 p-value")
print(p)
dev.off()

#########################
# COMPARISON WITH BOS PAPER
#########################
# list genes from paper
top_53_genes <- c(
  "MMP8", "OLFM4", "RETN", "GPR84", "CEACAM1", "LCN2", "ZDHHC19",
  "ANKRD22", "CD177", "HP", "TCN1", "PFKFB2", "CD24", "LTF", "DHRS9",
  "CHCHD7", "EXOSC4", "SMPDL3A", "PLAC8", "PDSS1", "TFRC", "MS4A4A",
  "BPI", "CLEC5A", "ASPRV1", "PILRA", "TNFAIP2", "CRISPLD2", "IFIT2",
  "ATP2B1", "MAK", "FGL2", "LOC728392", "REM2", "CAMK1D", "VCAN",
  "BTNL8", "ST6GALNAC2", "CECR1", "CPVL", "RNASE6", "CFD", "RPS6KA5",
  "TREM1", "DPEP2", "KRT23", "HCAR2", "SULF2", "HAL", "TGFBI",
  "ADGRE3", "RBP7", "MME"
)

# keep genes in common
top_53_genes_tmp <- intersect(top_53_genes, sig_results$hgnc_symbol)
res_53_genes <- sig_results[sig_results$hgnc_symbol %in% top_53_genes_tmp, ]
rownames(res_53_genes) <- res_53_genes$hgnc_symbol

# order DE analysis results given results from paper
res_53_genes <- res_53_genes[top_53_genes_tmp, ]

# add sign
res_53_genes$sign <- res_53_genes$log2FoldChange < 0

# to keep order for the plot
res_53_genes$hgnc_symbol <- factor(res_53_genes$hgnc_symbol,
                                   levels = res_53_genes$hgnc_symbol)

# plot log fold change values and sign
pdf(paste(results_path, "bos_comp.pdf", sep = ""), width = fig_width*0.7, height = fig_height)
print(ggplot(res_53_genes, aes(x = hgnc_symbol,
                               y = log2FoldChange, fill = sign)) +
        geom_bar(stat = "identity") +
        coord_flip() +
        scale_x_discrete(limits = rev(levels(res_53_genes$hgnc_symbol))) +
        theme_bw()+
        theme(legend.position = "none") + ylab("Log2 Fold Change") + 
        xlab("HGNC Symbol"))
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

# with APACHE and SOFA scores
meta_data_tmp <- meta_data

meta_data_tmp$sofa[is.na(meta_data_tmp$sofa)] <- median(meta_data_tmp$sofa, na.rm = T)
meta_data_tmp$apacheiii[is.na(meta_data_tmp$apacheiii)] <- median(meta_data_tmp$apacheiii, na.rm = T)
meta_data_tmp$wbcmaxsaps[is.na(meta_data_tmp$wbcmaxsaps)] <- median(meta_data_tmp$wbcmaxsaps, na.rm = T)

meta_data_tmp$sofa_scaled <- scale(meta_data_tmp$sofa)
meta_data_tmp$apacheiii_scaled <- scale(meta_data_tmp$apacheiii)
meta_data_tmp$wbcmaxsaps_scaled <- scale(meta_data_tmp$wbcmaxsaps)

# build DESeq2 object
dds_sofa <- DESeqDataSetFromMatrix(countData = cnt_data, colData = meta_data_tmp,
                                  design = ~ LCA_label + agephi + Gender1 + sofa_scaled)

# choose the reference level for the factor of interest
dds_sofa$LCA_label <- relevel(dds_sofa$LCA_label, ref = "Hypo")

# run DESeq
dds_sofa <- DESeq(dds_sofa)

# extract results
res_sofa <- results(dds_sofa, contrast = c("LCA_label", "Hyper", "Hypo"))

# replace NA values with 1s
res_sofa$padj[is.na(res_sofa$padj)] <- 1

# sort the genes from lowest to highest given adjusted p-values
res_sofa <- res_sofa[order(res_sofa$padj, decreasing = F), ]

# compare with results not including sofa values as covariate
# 0.92 correlation
cor.test(res$stat, res_sofa[rownames(res), "stat"], method = "spearman")

# build DESeq2 object
dds_apacheiii <- DESeqDataSetFromMatrix(countData = cnt_data, colData = meta_data_tmp,
                                   design = ~ LCA_label + agephi + Gender1 + apacheiii_scaled)

# choose the reference level for the factor of interest
dds_apacheiii$LCA_label <- relevel(dds_apacheiii$LCA_label, ref = "Hypo")

# run DESeq
dds_apacheiii <- DESeq(dds_apacheiii)

# extract results
res_apacheiii <- results(dds_apacheiii, contrast = c("LCA_label", "Hyper", "Hypo"))

# replace NA values with 1s
res_apacheiii$padj[is.na(res_apacheiii$padj)] <- 1

# sort the genes from lowest to highest given adjusted p-values
res_apacheiii <- res_apacheiii[order(res_apacheiii$padj, decreasing = F), ]

# compare with results not including apache iii values as covariate
# 0.92 correlation
cor.test(res$stat, res_apacheiii[rownames(res), "stat"], method = "spearman")

# build DESeq2 object
dds_wbc <- DESeqDataSetFromMatrix(countData = cnt_data, colData = meta_data_tmp,
                                        design = ~ LCA_label + agephi + Gender1 + wbcmaxsaps_scaled)

# choose the reference level for the factor of interest
dds_wbc$LCA_label <- relevel(dds_wbc$LCA_label, ref = "Hypo")

# run DESeq
dds_wbc <- DESeq(dds_wbc)

# extract results
res_wbc <- results(dds_wbc, contrast = c("LCA_label", "Hyper", "Hypo"))

# replace NA values with 1s
res_wbc$padj[is.na(res_wbc$padj)] <- 1

# sort the genes from lowest to highest given adjusted p-values
res_wbc <- res_wbc[order(res_wbc$padj, decreasing = F), ]

# compare with results not including apache iii values as covariate
# 0.92 correlation
cor.test(res$stat, res_wbc[rownames(res), "stat"], method = "spearman")

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
pdf(paste(results_path, "ipa_cp_full_plot.pdf", sep = ""), width = fig_width*0.8, height = fig_height)
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
# WGCNA
#########################
wgcna_data <- t(assay(vsd))
wgcna_meta_data <- meta_data

# covs to display
covs_filt <- c("Gender1", "Race1", "agephi", "diabetes", "PrimaryDiag1", "BASEABG_FiO2_1","BASEABG_PAO2_1","BASEABG_PaCO2_1","PFRatio_Brus","SFRatio_Brus","maxtemp","lowestsbp_today",
               "hrmax","mechvent","minutevent","peep","plateaupressure","meanairwaypressure","sofa", "positiveculture", "ali_risk",
               "sepsiscat","BERLIN_ARDS",
               "VFD","onpressorssaps","rrmaxsaps","urineoutputsaps", "naminsaps", "hco3minsaps","creatininemaxsaps",                
               "albuminminsaps", "glucoseminsaps", "wbcmaxsaps","hctminsaps","plateletsminsaps", "sapsii", "apacheiii", "bmi", "bilirubin_fin", "PAI1", "IL6","sTNFr1","RAGE",                                 
               "IL8","ANG2","Proc_C","ICAM","in_hospital_death",
               "LCA_label", "frac_anc",
               "Caucasian_Ethnicity1","PrimaryDiagCat1","bacteremia",
               "viruspos")

cat_vars <- c("Gender1", "Race1", "PrimaryDiag1", "positiveculture", "ali_risk", "BERLIN_ARDS", "onpressorssaps", "mechvent", "diabetes",
              "sepsiscat", "LCA_label", "Caucasian_Ethnicity1","PrimaryDiagCat1", "in_hospital_death", "bacteremia", "viruspos")

covs_filt <- covs_filt[!(covs_filt %in% cat_vars)]

#wgcna_meta_data[, cat_vars] <- apply(wgcna_meta_data[, cat_vars], 2, function(x) as.numeric(as.factor(x)))

wgcna_meta_data <- wgcna_meta_data[, covs_filt]

powers_ = c(seq(4,10,by=1), seq(12,20, by=2))

# network topology analysis function for each set
power_table <- pickSoftThreshold(wgcna_data, powerVector=powers_, verbose = 2)

plot_cols <- c(2,5,6,7)
col_names <- c("Scale Free Topology Model Fit", "Mean connectivity", "Median connectivity",
             "Max connectivity")
colors_ <- c("black", "red")

# choose power
sizeGrWindow(9, 5)
par(mfrow = c(1,2))
cex1 = 0.9

# diag plots
plot(power_table$fitIndices[,1], -sign(power_table$fitIndices[,3])*power_table$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"))
text(power_table$fitIndices[,1], -sign(power_table$fitIndices[,3])*power_table$fitIndices[,2],
     labels=powers_,cex=cex1,col="red")
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")

plot(power_table$fitIndices[,1], power_table$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(power_table$fitIndices[,1], power_table$fitIndices[,5], labels=powers_, cex=cex1,col="red")

# choose 12 from plots
net_ <- blockwiseModules(wgcna_data, power = 12, TOMType = "unsigned", minModuleSize = 30,
                         reassignThreshold = 0, mergeCutHeight = 0.25,
                         numericLabels = TRUE, pamRespectsDendro = FALSE,
                         saveTOMs = TRUE,
                         saveTOMFileBase = "femaleMouseTOM",
                         verbose = 3, maxBlockSize=20000)

# define ns
n_genes <- ncol(wgcna_data)
n_samples <-  nrow(wgcna_data)

mod_colors <- labels2colors(net_$colors)

mod_eig <- moduleEigengenes(wgcna_data, mod_colors)$eigengenes
mods_ <- orderMEs(mod_eig)
mod_trait_cor <- cor(mods_, wgcna_meta_data, use = "p")
mod_trait_pval <- corPvalueStudent(mod_trait_cor, nSamples)

# mod covs relationships
sizeGrWindow(10,6)

text_mat <- paste(signif(mod_trait_cor, 2), "\n(",
                   signif(mod_trait_pval, 1), ")", sep = "")

dim(text_mat) <- dim(mod_trait_cor)
par(mar = c(6, 8.5, 3, 3))

# plot
pdf(paste(results_path, "wgcna_mods_covs.pdf", sep = ""), width = 12, height = 10)
print(labeledHeatmap(Matrix = mod_trait_cor,
               xLabels = names(wgcna_meta_data),
               yLabels = names(mods_),
               ySymbols = names(mods_),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = text_mat,
               setStdMargins = T,
               cex.text = 0.5,
               zlim = c(-1,1), 
               main = paste("Module-trait relationships")))
dev.off()

# with LCA
meta_data$LCA_label <- factor(meta_data$LCA_label, levels=c("Hypo", "Hyper"))

lca_pvals <- apply(mods_, 2, function(x) summary(glm(meta_data$LCA_label~x, family="binomial"))$coefficients["x", "Pr(>|z|)"])
lca_est <- apply(mods_, 2, function(x) summary(glm(meta_data$LCA_label~x, family="binomial"))$coefficients["x", "Estimate"])

lca_adjp <- p.adjust(lca_pvals, "BH")

# in_hospital_death bacteremia viruspos
meta_data$in_hospital_death <- factor(meta_data$in_hospital_death, levels=c("0", "1"))
meta_data$bacteremia <- factor(meta_data$bacteremia, levels=c("0", "1"))

death_pvals <- apply(mods_, 2, function(x) summary(glm(meta_data$in_hospital_death~x, family="binomial"))$coefficients["x", "Pr(>|z|)"])
death_est <- apply(mods_, 2, function(x) summary(glm(meta_data$in_hospital_death~x, family="binomial"))$coefficients["x", "Estimate"])
death_adjp <- p.adjust(death_pvals, "BH")

bact_pvals <- apply(mods_, 2, function(x) summary(glm(meta_data$bacteremia~x, family="binomial"))$coefficients["x", "Pr(>|z|)"])
bact_est <- apply(mods_, 2, function(x) summary(glm(meta_data$bacteremia~x, family="binomial"))$coefficients["x", "Estimate"])
bact_adjp <- p.adjust(bact_pvals, "BH")

est_df <- data.frame(lca_est, death_est, bact_est)

lca_adjp_fmt <- lca_adjp
lca_adjp_fmt[lca_adjp<0.001] <- "P<.001"
lca_adjp_fmt[lca_adjp<0.01 & lca_adjp>=0.001] <- "P<.01"
lca_adjp_fmt[lca_adjp>=0.01] <- format(round(as.numeric(lca_adjp[lca_adjp>=0.01]), 2), nsmall = 2)
lca_adjp_fmt <- paste0(paste(format(round(as.numeric(lca_est), 2), nsmall = 2), lca_adjp_fmt, sep="\n("), ")")

death_adjp_fmt <- death_adjp
death_adjp_fmt[death_adjp<0.001] <- "P<.001"
death_adjp_fmt[death_adjp<0.01 & death_adjp>=0.001] <- "P<.01"
death_adjp_fmt[death_adjp>=0.01] <- format(round(as.numeric(death_adjp[death_adjp>=0.01]), 2), nsmall = 2)
death_adjp_fmt <- paste0(paste(format(round(as.numeric(death_est), 2), nsmall = 2), death_adjp_fmt, sep="\n("), ")")

bact_adjp_fmt <- bact_adjp
bact_adjp_fmt[bact_adjp<0.001] <- "P<.001"
bact_adjp_fmt[bact_adjp<0.01 & bact_adjp>=0.001] <- "P<.01"
bact_adjp_fmt[bact_adjp>=0.01] <- format(round(as.numeric(bact_adjp[bact_adjp>=0.01]), 2), nsmall = 2)
bact_adjp_fmt <- paste0(paste(format(round(as.numeric(bact_est), 2), nsmall = 2), bact_adjp_fmt, sep="\n("), ")")

text_df <- data.frame(lca_adjp_fmt, death_adjp_fmt, bact_adjp_fmt)

pdf(paste(results_path, "wgcna_mods_LCA_label.pdf", sep = ""), width = 7, height = 10)
print(labeledHeatmap(Matrix = est_df[, c(1:3)],
                     xLabels = c("LCA_label", "in_hospital_death", "bacteremia"),
                     yLabels = names(mods_),
                     ySymbols = names(mods_),
                     colorLabels = FALSE,
                     colors = greenWhiteRed(50),
                     textMatrix = text_df[, c(1:3)],
                     setStdMargins = T,
                     cex.text = 0.5,
                     zlim = c(-20,20),
                     main = paste("Module-trait relationships")))
dev.off()

# higher in Hyper
green_genes <- ensembl_res$hgnc_symbol[match(colnames(wgcna_data)[mod_colors=="green"],
                              ensembl_res$ensembl_gene_id)]

green_df <- data.frame(colnames(wgcna_data)[mod_colors=="green"], green_genes, as.data.frame(cor(wgcna_data, mods_, use = "p"))[colnames(wgcna_data)[mod_colors=="green"], c("MEgreen")])

genes_oi <- c("KY", "LTF", "TCN1", "GPI", "NKG7", "ENPP4", "CLINT1", "CD63", "PSMA6", "MTHFD2", "BPI", "CD24", "RETN", "MMP8")

green_df[green_df$green_genes %in% genes_oi, ]

#########################
# COMPARISON WITH STAR
######################### 
cnt_data_STAR <- read.csv("/Users/lucileneyton/Box Sync/EARLI_plasma/data/raw/EARLI_star_pc_and_lincRNA_genecounts.qc.tsv", row.names = 1, sep = "\t")

# filter and format
cnt_data_STAR <- cnt_data_STAR[ , meta_data$HOST_PAXgene_filename]
colnames(cnt_data_STAR) <- rownames(meta_data)

cnt_data_STAR <- cnt_data_STAR[!duplicated(sapply(rownames(cnt_data_STAR), function(x) strsplit(x, "\\.")[[1]][1])), ]
rownames(cnt_data_STAR) <- sapply(rownames(cnt_data_STAR), function(x) strsplit(x, "\\.")[[1]][1])

cnt_data_STAR <- cnt_data_STAR[rownames(cnt_data_STAR) %in% rownames(cnt_data), ]

# build DESeq2 object
dds_STAR <- DESeqDataSetFromMatrix(countData = cnt_data_STAR, colData = meta_data,
                              design = ~ LCA_label + agephi + Gender1)

# choose the reference level for the factor of interest
dds_STAR$LCA_label <- relevel(dds_STAR$LCA_label, ref = "Hypo")

# run DESeq
dds_STAR <- DESeq(dds_STAR)

# extract results
res_STAR <- results(dds_STAR, contrast = c("LCA_label", "Hyper", "Hypo"))

# sort the genes from lowest to highest given adjusted p-values
res_STAR <- res_STAR[order(res_STAR$padj, decreasing = F), ]

# replace NA values with 1s
res_STAR$padj[is.na(res_STAR$padj)] <- 1

kallisto_results <- data.frame(res)

dim(res_STAR[res_STAR$padj<0.05 , ])
dim(kallisto_results[kallisto_results$padj<0.05, ])
length(intersect(rownames(res_STAR[res_STAR$padj<0.05 , ]), rownames(kallisto_results[kallisto_results$padj<0.05, ])))

cor.test(kallisto_results[intersect(rownames(kallisto_results), rownames(res_STAR)), ]$log2FoldChange,
         res_STAR[intersect(rownames(kallisto_results), rownames(res_STAR)), ]$log2FoldChange, method="spearman")

plot(kallisto_results[intersect(rownames(kallisto_results), rownames(res_STAR)), ]$log2FoldChange,
     res_STAR[intersect(rownames(kallisto_results), rownames(res_STAR)), ]$log2FoldChange,
     xlab="kallisto", ylab="STAR", main="Log Fold-changes")




