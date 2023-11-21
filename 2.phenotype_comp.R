######################################################
# 2.phenotype_comp.R
# created on Nov 19 2020
# lucile.neyton@ucsf.edu
######################################################

rm(list = ls())
setwd("/Users/lucileneyton/Box Sync/EARLI_VALID/de_analysis_sepsis/")

# load libraries
library(ggplot2) 
library(DESeq2)
library(binom)
library(fgsea)
library(ggpubr)

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

file_prefix <- "sepsis_"

# set paths
data_path <- "/Users/lucileneyton/Box Sync/EARLI_VALID/data/"
results_path <-
  "/Users/lucileneyton/Box Sync/EARLI_VALID/de_analysis_sepsis/results/"

# functions
# generates a scatter plot of log2FC gene expression values
# between two datasets
# df1 and df2 are the two datasets of interest
# df1_labels and df2_labels are the corresponding meta data
# df1_group_label and df2_group_label are the column names 
#   for the group variables used to compute log2FC values
# df1_group1 and df2_group1 are the group variable first value
#   used to compute log2FC values
# df1_group2 and df2_group2 are the group variable second value
#   used to compute log2FC values
# plot controls whether a correlation plot is produced
# plot_label is the name of the plot
# df1_name and df2_name are the names of the datasets
log2fc_corr <- function(df1, df2, df1_labels, df2_labels,
                        df1_group_label, df2_group_label,
                        df1_group1, df1_group2,
                        df2_group1 = NULL, df2_group2 = NULL,
                        plot = FALSE, plot_label = NULL,
                        df1_name = NULL, df2_name = NULL,
                        logtr_data=F) {
  # if the group labels for the second dataset are not listed
  #   use the same ones as the first dataset
  if (is.null(df2_group1) & is.null(df2_group2)) {
    df2_group1 <- df1_group1
    df2_group2 <- df1_group2
  }
  
  # split data per group value
  # first group value
  df2_data_1 <- df2[, !is.na(df2_labels[, df2_group_label]) &
                      (df2_labels[, df2_group_label] == df2_group1)]
  df1_data_1 <- df1[, !is.na(df1_labels[, df1_group_label]) &
                      (df1_labels[, df1_group_label] == df1_group1)]
  
  # second group value
  df2_data_2 <- df2[, !is.na(df2_labels[, df2_group_label]) &
                      (df2_labels[, df2_group_label] == df2_group2)]
  df1_data_2 <- df1[, !is.na(df1_labels[, df1_group_label]) &
                      (df1_labels[, df1_group_label] == df1_group2)]
  
  # calculate average expression values per group
  df2_means_1 <- rowMeans(df2_data_1)
  df1_means_1 <- rowMeans(df1_data_1)
  
  df2_means_2 <- rowMeans(df2_data_2)
  df1_means_2 <- rowMeans(df1_data_2)
  
  # compute log2FC
  if (!logtr_data){
    df2_log2fc <- log2(df2_means_1 / df2_means_2)
  }else{
    print("diff")
    df2_log2fc <- df2_means_1 - df2_means_2
  }
  df1_log2fc <- log2(df1_means_1 / df1_means_2)
  
  df_log2fc <- data.frame("df2" = df2_log2fc,
                          "df1" = df1_log2fc[names(df2_log2fc)])
  
  # compute Spearman's correlation value
  cor_val <- cor(df_log2fc$df1, df_log2fc$df2, method = "spearman")
  
  if (plot) {
    if (!is.null(df2_group1) & !is.null(df2_group2)) {
      x_label <- paste(paste(df1_name, paste(df1_group1, df1_group2, sep=" / ")), "Log2 Fold Change")
      y_label <- paste(paste(df2_name, paste(df2_group1, df2_group2, sep=" / ")), "Log2 Fold Change")
    }else{
      x_label <- paste(paste(df1_name, paste(df1_group1, df1_group2, sep=" / ")), "Log2 Fold Change")
      y_label <- paste(paste(df2_name, paste(df1_group1, df1_group2, sep=" / ")), "Log2 Fold Change")
    }
    
    # compare rankings
    pdf(plot_label, width = 0.6*fig_width, height = 0.6*fig_height,
        fonts = "Helvetica", pointsize = 8)
    print(ggscatter(df_log2fc,
                    x = "df1", y = "df2",
                    xlab = x_label,
                    ylab = y_label
    ) +
      stat_cor(method="spearman", aes(label = ..r.label..), size = 3, family="Helvetica") +
      font("title", size = 8, family="Helvetica")+
      font("subtitle", size = 8, family="Helvetica")+
      font("caption", size = 8, family="Helvetica")+
      font("xlab", size = 8, family="Helvetica")+
      font("ylab", size = 8, family="Helvetica")+
      font("xy.text", size = 8, family="Helvetica")+
      font("legend.title", size = 8, family="Helvetica")+
      font("legend.text", size = 8, family="Helvetica")
    )
    dev.off()
  }
  
  return(list(cor_val=cor_val, df_log2fc=df_log2fc))
}

# generate correlation values from permuted datasets
# parameters are the same as for the log2fc_corr function
# nperm is the number of permuted sets to generate
# seed is the random state, specified by the user
perm_log2fc_corr <- function(df1, df2, df1_labels, df2_labels,
                             df1_group_label, df2_group_label,
                             df1_group1, df1_group2,
                             df2_group1 = NULL, df2_group2 = NULL,
                             nperm, seed) {
  cor_vals <- c()
  
  set.seed(seed)
  # for each permutation
  for (i in c(1:nperm)) {
    target_labels <- df2_labels
    # shuffle the labels of the second dataframe
    target_labels[, df2_group_label] <- sample(target_labels[, df2_group_label])
    
    # call the correlation function using the shuffled labels
    cor_val <- log2fc_corr(
      df1, df2, df1_labels, target_labels,
      df1_group_label, df2_group_label, df1_group1, df1_group2,
      df2_group1, df2_group2
    )[["cor_val"]]
    cor_vals <- c(cor_vals, cor_val)
  }
  
  return(cor_vals)
}

# compute empirical p-value
# corr_val is the value obtained from the original datasets
# perm_cor_vals are the permuted correlation values
# two-sided
emp_pval <- function(cor_val, perm_cor_vals) {
  cor_val <- cor_val[1]
  
  greater_equal_cnt <- sum(abs(perm_cor_vals) >= abs(cor_val))
  
  print(binom.confint(greater_equal_cnt, length(perm_cor_vals),
                      methods = "wilson"))
  
  return(greater_equal_cnt / length(perm_cor_vals))
}

#########################
# DATA LOADING
#########################
# list data files
# read EARLI data in
cnt_data_path <- paste(data_path, "raw/earli_counts_kallisto.csv", sep = "")
cnt_data <- read.csv(cnt_data_path, row.names = 1)
earli_labels <- read.csv(paste(data_path, "processed/covs_df_sepsis_cleaned.csv",
                                    sep = ""),
                              row.names = 1)

# external gene expression data
mars_data_path <- paste(data_path, "mars_paper/geneLevel_data.csv", sep = "")
srs_data_path <- paste(data_path, "srs_paper/geneLevel_data.csv", sep = "")

mars_data <- read.csv(mars_data_path, row.names = 1)
srs_data <- read.csv(srs_data_path, row.names = 1)

# external labels and/or metadata
mars_orig_labels_path <-
  paste(data_path, "mars_paper/class_labels.csv", sep = "")
mars2_orig_labels_path <-
  paste(data_path, "raw/MARS_subphenotype_classes_LDB20220301.csv", sep = "")
srs_orig_labels_path <-
  paste(data_path, "srs_paper/class_labels.csv", sep = "")

mars_orig_labels <- read.csv(mars_orig_labels_path, row.names = 1)
mars2_orig_labels <- read.csv(mars2_orig_labels_path, row.names = 1, sep=";")
srs_orig_labels <- read.csv(srs_orig_labels_path, row.names = 1)

#########################
# DATA PREPROCESSING
#########################
# make sure our variable of interest is a factor
earli_labels <- earli_labels[!is.na(earli_labels$LCA_label), ]
earli_labels$LCA_label <- as.factor(earli_labels$LCA_label)

# drop rows of meta data pertaining to samples without sequencing data
earli_labels <- earli_labels[rownames(earli_labels) %in% colnames(cnt_data), ]

# filter count data to exclude samples without a label
cnt_data <- cnt_data[, rownames(earli_labels)]

#########################
# DE ANALYSIS
#########################
# for mars
# keep only genes in common
cnt_data_filt_mars <- cnt_data[rownames(cnt_data) %in% colnames(mars_data), ]

# build DESeq2 object
dds_mars <- DESeqDataSetFromMatrix(countData = cnt_data_filt_mars,
                                   colData = earli_labels,
                                   design = ~ 1)

# choose the reference level for the factor of interest
dds_mars$LCA_label <- relevel(dds_mars$LCA_label, ref = "Hypo")

# run DESeq
dds_mars <- DESeq(dds_mars)

# transform data
vsd_mars <- vst(dds_mars, blind = TRUE)

# for srs
# keep only genes in common
cnt_data_filt_srs <- cnt_data[rownames(cnt_data) %in% colnames(srs_data), ]

# build DESeq2 object
dds_srs <- DESeqDataSetFromMatrix(countData = cnt_data_filt_srs,
                                  colData = earli_labels,
                                  design = ~ 1)

# choose the reference level for the factor of interest
dds_srs$LCA_label <- relevel(dds_srs$LCA_label, ref = "Hypo")

# run DESeq
dds_srs <- DESeq(dds_srs)

# transform data
vsd_srs <- vst(dds_srs, blind = TRUE)

#########################
# COMPARE EARLI FOLD-CHANGE VALUES TO MARS/SRS
#########################
##### MARS
# filter for genes in common
mars_data_filt <- t(mars_data[, colnames(mars_data) %in% rownames(vsd_mars)])

vsd_mars <- assay(vsd_mars)

# compare class labels
mars_orig_labels <- mars_orig_labels[!is.na(mars_orig_labels$subphenotype), ,
                                     drop = F]

# compare MARS4/2 vs Hyper/Hypo
corr_val <- log2fc_corr(vsd_mars, mars_data_filt,
                        earli_labels, mars_orig_labels,
                        "LCA_label", "subphenotype", "Hyper", "Hypo",
                        df2_group1 = "Mars4", df2_group2 = "Mars2",
                        plot = TRUE, plot_label = paste(results_path,
                                                        paste0(file_prefix, "log2fc_cor_mars42_hyperhypo.pdf"),
                                                        sep = ""),
                        df1_name = "EARLI", df2_name = "MARS", logtr_data=T
)
df_log2fc <- corr_val[["df_log2fc"]]
df_mars <- corr_val[["df_log2fc"]][, "df2"]
corr_val <- corr_val[["cor_val"]]

corr_vals <- perm_log2fc_corr(vsd_mars, mars_data_filt,
                              earli_labels, mars_orig_labels,
                              "LCA_label", "subphenotype", "Hyper", "Hypo",
                              df2_group1 = "Mars4", df2_group2 = "Mars2",
                              nperm = 1000, seed = 123
)

# generate p-value
res_p <- emp_pval(corr_val, corr_vals)

# which are DE in same or different direction?
# wilcox DE
# hyper vs hypo
labels_oi_earli_for_mars <- earli_labels[!is.na(earli_labels[, "LCA_label"]) & (earli_labels[, "LCA_label"] %in% c("Hypo", "Hyper")), "LCA_label", drop=F]
wilcox_df_earli_for_mars <- data.frame(t(vsd_mars[, rownames(labels_oi_earli_for_mars)]))

wilcox_tests_earli_for_mars <- sapply(colnames(wilcox_df_earli_for_mars), function(x) wilcox.test(as.numeric(wilcox_df_earli_for_mars[labels_oi_earli_for_mars$LCA_label=="Hypo", x]), as.numeric(wilcox_df_earli_for_mars[labels_oi_earli_for_mars$LCA_label=="Hyper", x]))$p.value)

wilcox_tests_adj_earli_for_mars <- p.adjust(wilcox_tests_earli_for_mars, "BH")

# mars4 vs mars2
labels_oi_mars <- mars_orig_labels[!is.na(mars_orig_labels[, "subphenotype"]) & (mars_orig_labels[, "subphenotype"] %in% c("Mars2", "Mars4")), "subphenotype", drop=F]

# numbers
print(dim(mars_orig_labels[!is.na(mars_orig_labels[, "subphenotype"]) & (mars_orig_labels[, "subphenotype"] %in% c("Mars2", "Mars4", "Mars3", "Mars1")), "subphenotype", drop=F]))
wilcox_df_mars <- data.frame(t(mars_data_filt[, rownames(labels_oi_mars)]))

wilcox_tests_mars <- sapply(colnames(wilcox_df_mars), function(x) wilcox.test(as.numeric(wilcox_df_mars[labels_oi_mars$subphenotype=="Mars2", x]), as.numeric(wilcox_df_mars[labels_oi_mars$subphenotype=="Mars4", x]))$p.value)

wilcox_tests_adj_mars <- p.adjust(wilcox_tests_mars, "BH")

# pie chart
wilcox_tests_adj_earli_for_mars <- wilcox_tests_adj_earli_for_mars[rownames(df_log2fc)]
wilcox_tests_adj_mars <- wilcox_tests_adj_mars[rownames(df_log2fc)]

de_samedir_mars <- sum(sign(df_log2fc$df1)==sign(df_log2fc$df2) & wilcox_tests_adj_earli_for_mars<0.05 & wilcox_tests_adj_mars<0.05)
notdeboth_samedir_mars <- sum(sign(df_log2fc$df1)==sign(df_log2fc$df2) & !(wilcox_tests_adj_earli_for_mars<0.05 & wilcox_tests_adj_mars<0.05))
diffdir_mars <- sum(sign(df_log2fc$df1)!=sign(df_log2fc$df2))

pie_df_mars <- data.frame(value = c(de_samedir_mars, notdeboth_samedir_mars, diffdir_mars),
                group = c("DE Same Direction","Non-DE Same Direction","Different Direction"))

pdf(paste0(results_path, "pie_mars.pdf"), height = 2, width=3)
print(ggplot(pie_df_mars, aes(x = "", y = value, fill = group, alpha=0.8)) +
        scale_alpha(range = c(0, 1)) +
        guides(alpha = FALSE) +
 geom_col(color = "black") +
 geom_text(aes(label = value),
           position = position_stack(vjust = 0.5),
           alpha=1) +
 coord_polar(theta = "y") +
   scale_fill_manual(breaks=c("DE Same Direction","Non-DE Same Direction","Different Direction"), values = c("red", "pink", "white")) +
  theme_void() + 
  ggtitle(paste0("N=", length(wilcox_tests_adj_earli_for_mars))) +
  theme(plot.title = element_text(hjust = 1), legend.title= element_blank()) +
   guides(fill = guide_legend(override.aes = list(alpha = 0.6))))
dev.off()

df_log2fc$cat[sign(df_log2fc$df1)==sign(df_log2fc$df2) & wilcox_tests_adj_earli_for_mars<0.05 & wilcox_tests_adj_mars<0.05] <- "DE Same Direction"
df_log2fc$cat[sign(df_log2fc$df1)==sign(df_log2fc$df2) & !(wilcox_tests_adj_earli_for_mars<0.05 & wilcox_tests_adj_mars<0.05)] <- "Non-DE Same Direction"
df_log2fc$cat[sign(df_log2fc$df1)!=sign(df_log2fc$df2)] <- "Different Direction"
df_log2fc$cat <- factor(df_log2fc$cat)

# compare rankings
pdf("mars24_color.pdf", width = 0.6*fig_width, height = 0.6*fig_height, fonts = "Helvetica", pointsize = 8)
print(ggscatter(df_log2fc, x = "df1", y = "df2", fill="cat", color="black",
                xlab = paste(paste("EARLI", paste("Hyper", "Hypo", sep=" / ")), "Log2 Fold Change"),
                ylab = paste(paste("MARS", paste("Mars4", "Mars2", sep=" / ")), "Log2 Fold Change"),
                title="MARS vs. EARLI",
                shape=21)+
  scale_fill_manual(breaks=c("DE Same Direction","Non-DE Same Direction","Different Direction"), values = c("red", "pink", "white"))+
  stat_cor(method="spearman", aes(label = ..r.label..), size = 3, family="Helvetica") +
  font("title", size = 12, family="Helvetica")+
  font("subtitle", size = 8, family="Helvetica")+
  font("caption", size = 8, family="Helvetica")+
  font("xlab", size = 8, family="Helvetica")+
  font("ylab", size = 8, family="Helvetica")+
  font("xy.text", size = 8, family="Helvetica")+
  font("legend.title", size = 8, family="Helvetica")+
    font("legend.text", size = 8, family="Helvetica") + theme(legend.position = "none"))
dev.off()

##### MARS REACTIVE UNINFLAMED
# filter for genes in common
mars_orig_labels <- read.csv(mars_orig_labels_path, row.names = 1)
mars_orig_labels$gsm_id <- rownames(mars_orig_labels)
mars2_all_orig_labels <- merge(mars_orig_labels, mars2_orig_labels, by.x="orig_label", by.y="ICU_ID_from_datasource", all.x=T)

# compare class labels
mars2_all_orig_labels <- mars2_all_orig_labels[!is.na(mars2_all_orig_labels$predicted_cluster), ,
                                               drop = F]
mars2_all_orig_labels$predicted_cluster[mars2_all_orig_labels$predicted_cluster==1] <- "reactive"
mars2_all_orig_labels$predicted_cluster[mars2_all_orig_labels$predicted_cluster==2] <- "uninflamed"

# compare reactive/uninflamed vs Hyper/Hypo
mars2_data_filt <- mars_data_filt[, colnames(mars_data_filt) %in% mars2_all_orig_labels$gsm_id]
rownames(mars2_all_orig_labels) <- mars2_all_orig_labels$gsm_id
mars2_all_orig_labels <- mars2_all_orig_labels[colnames(mars2_data_filt), ]
corr_val <- log2fc_corr(vsd_mars, mars2_data_filt[, rownames(mars2_all_orig_labels)],
                        earli_labels, mars2_all_orig_labels,
                        "LCA_label", "predicted_cluster", "Hyper", "Hypo",
                        df2_group1 = "reactive", df2_group2 = "uninflamed",
                        plot = TRUE, plot_label = paste(results_path,
                                                        paste0(file_prefix, "log2fc_cor_marsru_hyperhypo.pdf"),
                                                        sep = ""),
                        df1_name = "EARLI", df2_name = "MARS", logtr_data=T
)

df_log2fc <- corr_val[["df_log2fc"]]
df_mars2 <- corr_val[["df_log2fc"]][, "df2"]
corr_val <- corr_val[["cor_val"]]

corr_vals <- perm_log2fc_corr(vsd_mars, mars2_data_filt[, rownames(mars2_all_orig_labels)],
                              earli_labels, mars2_all_orig_labels,
                              "LCA_label", "predicted_cluster", "Hyper", "Hypo",
                              df2_group1 = "reactive", df2_group2 = "uninflamed",
                              nperm = 1000, seed = 123
)

# generate p-value
res_p <- emp_pval(corr_val, corr_vals)

# which are DE in same or different direction?
# wilcox DE
# reactive vs uninflamed
labels_oi_mars2 <- mars2_all_orig_labels[!is.na(mars2_all_orig_labels[, "predicted_cluster"]) & (mars2_all_orig_labels[, "predicted_cluster"] %in% c("uninflamed", "reactive")), "predicted_cluster", drop=F]
print(dim(labels_oi_mars2))
wilcox_df_mars2 <- data.frame(t(mars_data_filt[, rownames(labels_oi_mars2)]))

wilcox_tests_mars2 <- sapply(colnames(wilcox_df_mars2), function(x) wilcox.test(as.numeric(wilcox_df_mars2[labels_oi_mars2$predicted_cluster=="uninflamed", x]), as.numeric(wilcox_df_mars2[labels_oi_mars2$predicted_cluster=="reactive", x]))$p.value)

wilcox_tests_adj_mars2 <- p.adjust(wilcox_tests_mars2, "BH")

# pie chart
wilcox_tests_adj_earli_for_mars2 <- wilcox_tests_adj_earli_for_mars[rownames(df_log2fc)]
wilcox_tests_adj_mars2 <- wilcox_tests_adj_mars2[rownames(df_log2fc)]

de_samedir_mars2 <- sum(sign(df_log2fc$df1)==sign(df_log2fc$df2) & wilcox_tests_adj_earli_for_mars2<0.05 & wilcox_tests_adj_mars2<0.05)
notdeboth_samedir_mars2 <- sum(sign(df_log2fc$df1)==sign(df_log2fc$df2) & !(wilcox_tests_adj_earli_for_mars2<0.05 & wilcox_tests_adj_mars2<0.05))
diffdir_mars2 <- sum(sign(df_log2fc$df1)!=sign(df_log2fc$df2))

pie_df_mars2 <- data.frame(value = c(de_samedir_mars2, notdeboth_samedir_mars2, diffdir_mars2),
                          group = c("DE Same Direction","Non-DE Same Direction","Different Direction"))

pdf(paste0(results_path, "pie_mars2.pdf"), height = 2, width=3)
print(ggplot(pie_df_mars2, aes(x = "", y = value, fill = group, alpha=0.8)) +
        scale_alpha(range = c(0, 1)) +
        guides(alpha = FALSE) +
        geom_col(color = "black") +
        geom_text(aes(label = value),
                  position = position_stack(vjust = 0.5),
                  alpha=1) +
        coord_polar(theta = "y") +
        scale_fill_manual(breaks=c("DE Same Direction","Non-DE Same Direction","Different Direction"), values = c("red", "pink", "white")) +
        theme_void() + 
        ggtitle(paste0("N=", length(wilcox_tests_adj_earli_for_mars2))) +
        theme(plot.title = element_text(hjust = 1), legend.title= element_blank()) +
        guides(fill = guide_legend(override.aes = list(alpha = 0.6))))
dev.off()

df_log2fc$cat[sign(df_log2fc$df1)==sign(df_log2fc$df2) & wilcox_tests_adj_earli_for_mars2<0.05 & wilcox_tests_adj_mars2<0.05] <- "DE Same Direction"
df_log2fc$cat[sign(df_log2fc$df1)==sign(df_log2fc$df2) & !(wilcox_tests_adj_earli_for_mars2<0.05 & wilcox_tests_adj_mars2<0.05)] <- "Non-DE Same Direction"
df_log2fc$cat[sign(df_log2fc$df1)!=sign(df_log2fc$df2)] <- "Different Direction"
df_log2fc$cat <- factor(df_log2fc$cat)

pdf("marsru_color.pdf", width = 0.6*fig_width, height = 0.6*fig_height, fonts = "Helvetica", pointsize = 8)
print(ggscatter(df_log2fc, x = "df1", y = "df2", fill="cat", color="black",
                xlab = paste(paste("EARLI", paste("Hyper", "Hypo", sep=" / ")), "Log2 Fold Change"),
                ylab = paste(paste("MARS", paste("reactive", "uninflamed", sep=" / ")), "Log2 Fold Change"),
                title="Reactive/Uninflamed vs. EARLI",
                shape=21)+
        scale_fill_manual(breaks=c("DE Same Direction","Non-DE Same Direction","Different Direction"), values = c("red", "pink", "white"))+
        stat_cor(method="spearman", aes(label = ..r.label..), size = 3, family="Helvetica") +
        font("title", size = 12, family="Helvetica")+
        font("subtitle", size = 8, family="Helvetica")+
        font("caption", size = 8, family="Helvetica")+
        font("xlab", size = 8, family="Helvetica")+
        font("ylab", size = 8, family="Helvetica")+
        font("xy.text", size = 8, family="Helvetica")+
        font("legend.title", size = 8, family="Helvetica")+
        font("legend.text", size = 8, family="Helvetica") + theme(legend.position = "none"))
dev.off()

##### SRS
# filter for genes in common
srs_data_filt <- t(srs_data[, colnames(srs_data) %in% rownames(vsd_srs)])

vsd_srs <- assay(vsd_srs)

# compare class labels
srs_orig_labels <- srs_orig_labels[!is.na(srs_orig_labels$subphenotype), ,
                                   drop = F]

# compare SRS1/2 vs Hyper/Hypo
corr_val <- log2fc_corr(vsd_srs, srs_data_filt, earli_labels, srs_orig_labels,
                        "LCA_label", "subphenotype", "Hyper", "Hypo",
                        df2_group1 = 1, df2_group2 = 2,
                        plot = TRUE, plot_label =
                          paste(results_path, paste0(file_prefix, "log2fc_cor_srs12_hyperhypo.pdf"), sep = ""),
                        df1_name = "EARLI", df2_name = "SRS"
)

df_log2fc <- corr_val[["df_log2fc"]]
df_srs <- corr_val[["df_log2fc"]][, "df2"]
corr_val <- corr_val[["cor_val"]]

corr_vals <- perm_log2fc_corr(vsd_srs, srs_data_filt,
                              earli_labels, srs_orig_labels,
                              "LCA_label", "subphenotype", "Hyper", "Hypo",
                              df2_group1 = 1, df2_group2 = 2,
                              nperm = 1000, seed = 123
)

# generate p-value
res_p <- emp_pval(corr_val, corr_vals)

# which are DE in same or different direction?
# wilcox DE
# hyper vs hypo
labels_oi_earli_for_srs <- earli_labels[!is.na(earli_labels[, "LCA_label"]) & (earli_labels[, "LCA_label"] %in% c("Hypo", "Hyper")), "LCA_label", drop=F]
wilcox_df_earli_for_srs <- data.frame(t(vsd_srs[, rownames(labels_oi_earli_for_srs)]))

wilcox_tests_earli_for_srs <- sapply(colnames(wilcox_df_earli_for_srs), function(x) wilcox.test(as.numeric(wilcox_df_earli_for_srs[labels_oi_earli_for_srs$LCA_label=="Hypo", x]), as.numeric(wilcox_df_earli_for_srs[labels_oi_earli_for_srs$LCA_label=="Hyper", x]))$p.value)

wilcox_tests_adj_earli_for_srs <- p.adjust(wilcox_tests_earli_for_srs, "BH")

# srs 1 vs 2
labels_oi_srs <- srs_orig_labels[!is.na(srs_orig_labels[, "subphenotype"]) & (srs_orig_labels[, "subphenotype"] %in% c(2, 1)), "subphenotype", drop=F]
print(dim(labels_oi_srs))

wilcox_df_srs <- data.frame(t(srs_data_filt[, rownames(labels_oi_srs)]))

wilcox_tests_srs <- sapply(colnames(wilcox_df_srs), function(x) wilcox.test(as.numeric(wilcox_df_srs[labels_oi_srs$subphenotype==2, x]), as.numeric(wilcox_df_srs[labels_oi_srs$subphenotype==1, x]))$p.value)

wilcox_tests_adj_srs <- p.adjust(wilcox_tests_srs, "BH")

# pie chart
wilcox_tests_adj_earli_for_srs <- wilcox_tests_adj_earli_for_srs[rownames(df_log2fc)]
wilcox_tests_adj_srs <- wilcox_tests_adj_srs[rownames(df_log2fc)]

de_samedir_srs <- sum(sign(df_log2fc$df1)==sign(df_log2fc$df2) & wilcox_tests_adj_earli_for_srs<0.05 & wilcox_tests_adj_srs<0.05)
notdeboth_samedir_srs <- sum(sign(df_log2fc$df1)==sign(df_log2fc$df2) & !(wilcox_tests_adj_earli_for_srs<0.05 & wilcox_tests_adj_srs<0.05))
diffdir_srs <- sum(sign(df_log2fc$df1)!=sign(df_log2fc$df2))

pie_df_srs <- data.frame(value = c(de_samedir_srs, notdeboth_samedir_srs, diffdir_srs),
                          group = c("DE Same Direction","Non-DE Same Direction","Different Direction"))

pie_df_srs$group <- factor(pie_df_srs$group, levels=c("DE Same Direction","Non-DE Same Direction","Different Direction"))

pdf(paste0(results_path, "pie_srs.pdf"), height = 2, width=3)
print(ggplot(pie_df_srs, aes(x = "", y = value, fill = group, alpha=0.8)) +
        scale_alpha(range = c(0, 1)) +
        guides(alpha = FALSE) +
        geom_col(color = "black") +
        geom_text(aes(label = value),
                  position = position_stack(vjust = 0.5),
                  alpha=1) +
        coord_polar(theta = "y") +
        scale_fill_manual(breaks=c("DE Same Direction","Non-DE Same Direction","Different Direction"), values = c("red", "pink", "white")) +
        theme_void() + 
        ggtitle(paste0("N=", length(wilcox_tests_adj_earli_for_srs))) +
        theme(plot.title = element_text(hjust = 1), legend.title= element_blank()) +
        guides(fill = guide_legend(override.aes = list(alpha = 0.6))))
dev.off()

df_log2fc$cat[sign(df_log2fc$df1)==sign(df_log2fc$df2) & wilcox_tests_adj_earli_for_srs<0.05 & wilcox_tests_adj_srs<0.05] <- "DE Same Direction"
df_log2fc$cat[sign(df_log2fc$df1)==sign(df_log2fc$df2) & !(wilcox_tests_adj_earli_for_srs<0.05 & wilcox_tests_adj_srs<0.05)] <- "Non-DE Same Direction"
df_log2fc$cat[sign(df_log2fc$df1)!=sign(df_log2fc$df2)] <- "Different Direction"
df_log2fc$cat <- factor(df_log2fc$cat)

pdf("srs12_color.pdf", width = 0.6*fig_width, height = 0.6*fig_height, fonts = "Helvetica", pointsize = 8)
print(ggscatter(df_log2fc, x = "df1", y = "df2", fill="cat", color="black",
                xlab = paste(paste("EARLI", paste("Hyper", "Hypo", sep=" / ")), "Log2 Fold Change"),
                ylab = paste(paste("SRS", paste("1", "2", sep=" / ")), "Log2 Fold Change"),
                title="SRS vs. EARLI",
                shape=21)+
        scale_fill_manual(breaks=c("DE Same Direction","Non-DE Same Direction","Different Direction"), values = c("red", "pink", "white"))+
        stat_cor(method="spearman", aes(label = ..r.label..), size = 3, family="Helvetica") +
        font("title", size = 12, family="Helvetica")+
        font("subtitle", size = 8, family="Helvetica")+
        font("caption", size = 8, family="Helvetica")+
        font("xlab", size = 8, family="Helvetica")+
        font("ylab", size = 8, family="Helvetica")+
        font("xy.text", size = 8, family="Helvetica")+
        font("legend.title", size = 8, family="Helvetica")+
        font("legend.text", size = 8, family="Helvetica") + theme(legend.position = "none"))
dev.off()

#########################
# ONLY CONSIDER DIFFERENTIALLY EXPRESSED GENES
#########################
dgea_data <- read.csv(paste0(results_path, "DGEA_results.csv"), row.names = 1)

vsd_mars_de <- vsd_mars[rownames(vsd_mars) %in% rownames(dgea_data), ]
vsd_srs_de <-  vsd_srs[rownames(vsd_srs) %in% rownames(dgea_data), ]

mars_data_filt_de <- mars_data_filt[rownames(mars_data_filt) %in% rownames(dgea_data), ]
srs_data_filt_de <- srs_data_filt[rownames(srs_data_filt) %in% rownames(dgea_data), ]

# compare MARS4/2 vs Hyper/Hypo
corr_val <- log2fc_corr(vsd_mars_de, mars_data_filt_de,
                        earli_labels, mars_orig_labels,
                        "LCA_label", "subphenotype", "Hyper", "Hypo",
                        df2_group1 = "Mars2", df2_group2 = "Mars4",
                        plot = TRUE, plot_label = paste(results_path,
                                                        paste0(file_prefix, "log2fc_cor_mars24_hyperhypo_de.pdf"),
                                                        sep = ""),
                        df1_name = "EARLI", df2_name = "MARS", logtr_data=T
)
df_log2fc <- corr_val[["df_log2fc"]]
corr_val <- corr_val[["cor_val"]]

corr_vals <- perm_log2fc_corr(vsd_mars_de, mars_data_filt_de,
                              earli_labels, mars_orig_labels,
                              "LCA_label", "subphenotype", "Hyper", "Hypo",
                              df2_group1 = "Mars2", df2_group2 = "Mars4",
                              nperm = 1000, seed = 123
)

# generate p-value
print(emp_pval(corr_val, corr_vals))

# compare reactive/uninflamed vs Hyper/Hypo
mars2_data_filt <- mars_data_filt_de[, colnames(mars_data_filt_de) %in% mars2_all_orig_labels$gsm_id]
rownames(mars2_all_orig_labels) <- mars2_all_orig_labels$gsm_id
mars2_all_orig_labels <- mars2_all_orig_labels[colnames(mars2_data_filt), ]
corr_val <- log2fc_corr(vsd_mars_de, mars2_data_filt[, rownames(mars2_all_orig_labels)],
                        earli_labels, mars2_all_orig_labels,
                        "LCA_label", "predicted_cluster", "Hyper", "Hypo",
                        df2_group1 = "reactive", df2_group2 = "uninflamed",
                        plot = TRUE, plot_label = paste(results_path,
                                                        paste0(file_prefix, "log2fc_cor_marsru_hyperhypo_de.pdf"),
                                                        sep = ""),
                        df1_name = "EARLI", df2_name = "MARS", logtr_data=T
)

df_log2fc <- corr_val[["df_log2fc"]]
corr_val <- corr_val[["cor_val"]]

corr_vals <- perm_log2fc_corr(vsd_mars_de, mars2_data_filt[, rownames(mars2_all_orig_labels)],
                              earli_labels, mars2_all_orig_labels,
                              "LCA_label", "predicted_cluster", "Hyper", "Hypo",
                              df2_group1 = "reactive", df2_group2 = "uninflamed",
                              nperm = 1000, seed = 123
)
print(emp_pval(corr_val, corr_vals))

# compare SRS1/2 vs Hyper/Hypo
corr_val <- log2fc_corr(vsd_srs_de, srs_data_filt_de, earli_labels, srs_orig_labels,
                        "LCA_label", "subphenotype", "Hyper", "Hypo",
                        df2_group1 = 1, df2_group2 = 2,
                        plot = TRUE, plot_label =
                          paste(results_path, paste0(file_prefix, "log2fc_cor_srs12_hyperhypo_de.pdf"), sep = ""),
                        df1_name = "EARLI", df2_name = "SRS"
)

df_log2fc <- corr_val[["df_log2fc"]]
corr_val <- corr_val[["cor_val"]]

corr_vals <- perm_log2fc_corr(vsd_srs_de, srs_data_filt_de,
                              earli_labels, srs_orig_labels,
                              "LCA_label", "subphenotype", "Hyper", "Hypo",
                              df2_group1 = 1, df2_group2 = 2,
                              nperm = 1000, seed = 123
)
print(emp_pval(corr_val, corr_vals))

#########################
# GSEA ON WILCOXON-DERIVED P-VALUES
#########################
# IPA res
hypervshypo_cp_sig <- read.csv(paste(results_path, "hypervshypo_cp_sig.csv", sep = ""))

ipa_res_list <- list("mars"=read.csv(paste0(results_path, "ipa_mars.txt"), sep = "\t", skip = 2),
                     "mars2"=read.csv(paste0(results_path, "ipa_mars2.txt"), sep = "\t", skip = 2),
                     "srs"=read.csv(paste0(results_path, "ipa_srs.txt"), sep = "\t", skip = 2))

for (ext_ds_name in names(ipa_res_list)){
  ext_ds <- ipa_res_list[[ext_ds_name]]
  
  ext_ds_filt <- ext_ds[ext_ds$Ingenuity.Canonical.Pathways %in% hypervshypo_cp_sig$Ingenuity.Canonical.Pathways, ]
  
  for (path_ in hypervshypo_cp_sig$Ingenuity.Canonical.Pathways){
    if (!(path_ %in% ext_ds_filt$Ingenuity.Canonical.Pathways)){
      ext_ds_filt <- rbind(ext_ds_filt, c("Ingenuity.Canonical.Pathways"=path_, "X.log.p.value."=NA, "Ratio"=NA,  "z.score"=0, "Molecules"=NA, "X"=NA))
    }
  }
  
  ext_ds_filt$sig <- (abs(as.numeric(ext_ds_filt$z.score))>2)
  
  ext_ds_filt$cp <- factor(ext_ds_filt$Ingenuity.Canonical.Pathways, levels=hypervshypo_cp_sig[order(hypervshypo_cp_sig$zScore, decreasing = T), "Ingenuity.Canonical.Pathways"])
  
  ext_ds_filt$zscore_sign[ext_ds_filt$z.score > 0] <- "+"
  ext_ds_filt$zscore_sign[ext_ds_filt$z.score < 0] <- "-"
  
  ext_ds_filt$zscore_sign <- as.factor(ext_ds_filt$zscore_sign)
  ext_ds_filt$zscore_sign <- relevel(ext_ds_filt$zscore_sign, "+")
  
  ext_ds_filt$alpha[ext_ds_filt$sig==T] <- 1
  ext_ds_filt$alpha[ext_ds_filt$sig==F] <- 0.5
  
  ext_ds_filt$z.score <- as.numeric(ext_ds_filt$z.score)
  
  # dot plot
  pdf(paste(results_path, paste0(ext_ds_name, "_ipa_cp_dotplot.pdf"), sep = ""), width = fig_width*0.8, height = fig_height)
  print(ggplot(ext_ds_filt, aes(x=z.score, y=cp)) +
          geom_point(aes(color=zscore_sign, alpha=alpha), size=3, show.legend=F) +
          scale_alpha(range = c(0.5, 1)) +
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
          scale_x_continuous(breaks = c(-5, -2.5, 0, 2.5, 5)) +
          scale_y_discrete(limits = rev(levels(ext_ds_filt$cp))))
  dev.off()
  
}










