######################################################
# metag_script.R
# created on Sept 1 2022
# lucile.neyton@ucsf.edu
######################################################

rm(list = ls())
setwd("/Users/lucileneyton/Box Sync/EARLI_VALID/de_analysis_sepsis/")

library(ggplot2)
library(dplyr)
library(tidyr)
library(data.table)
library(phyloseq)
library(ggpubr)
library(vegan)

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

data_path <- "/Users/lucileneyton/Box Sync/EARLI_VALID/data/"
results_path <- "/Users/lucileneyton/Box Sync/EARLI_VALID/de_analysis_sepsis/results/"

#########################
# DATA LOADING
#########################
# background-corrected data for alpha and beta diversity calculations
bgc_data <- read.csv(paste0(data_path, "/processed/all_bgc_taxa.csv"))

# meta data
meta_data <- read.csv(paste0(data_path, "/processed/covs_df_sepsis_cleaned.csv"), row.names = 1)

# microbial mass data
microbial_mass_data <- read.csv("/Users/lucileneyton/Box Sync/EARLI_clustering/data/processed/EARLI_data_merged_updated_6.2.22.csv")
rownames(microbial_mass_data) <- paste0("EARLI_", microbial_mass_data$Barcode)
microbial_mass_data <- microbial_mass_data[rownames(meta_data), ]

# taxonomic data from NCBI
# Sept 1, 2022
#tax_data <- classification(db = "ncbi", sci_id = bgc_data$tax_id)
#saveRDS(tax_data, paste0(data_path, "/processed/tax_data.RDS"))
tax_data <- readRDS(paste0(data_path, "/processed/tax_data.RDS"))

# format for phyloseq
ranklist_ <- c("superkingdom", "kingdom", "phylum", "class", "order", "family", "genus", "species", "taxon")
tax_data <- rbindlist(tax_data, idcol = T) %>%
  group_by(.id) %>%
  distinct(rank, .keep_all = T) %>%
  dplyr::select(-id) %>%
  pivot_wider(values_from = "name", names_from = "rank") %>%
  dplyr::rename(taxon = .id) %>%
  dplyr::select(all_of(ranklist_)) %>%
  ungroup

# rbm data for nt rpm and differential abundance analysis
rbm_data <- read.csv(paste0(data_path, "processed/all_rbm_results_w_gram.csv"))

rbm_data$EARLI_Barcode <- paste0("EARLI_", sapply(rbm_data$sample_name, function(x) strsplit(x, "_")[[1]][2]))
rbm_data <- rbm_data[rbm_data$sepsis_list, ]

rbm_data_cast <- dcast(rbm_data, genus~EARLI_Barcode, value.var="nt_rpm")
rbm_data_cast[is.na(rbm_data_cast)] <- 0
rownames(rbm_data_cast) <- rbm_data_cast$genus
rbm_data_cast <- rbm_data_cast[, colnames(rbm_data_cast)!="genus"]

# fill in with zeros for non-detected species
for (row_ in rownames(meta_data)){
  if (!(row_ %in% colnames(rbm_data_cast))){
    rbm_data_cast <- cbind(rbm_data_cast, rep(0, nrow(rbm_data_cast)))
    colnames(rbm_data_cast)[ncol(rbm_data_cast)] <- row_
  }
}

rbm_data_cast <- rbm_data_cast[, rownames(meta_data)]

#########################
# DATA PRE-PROCESSING
#########################
# format and melt
bgc_data$EARLI_Barcode <- sapply(bgc_data$sample_name, function(x) paste0("EARLI_", strsplit(x, "_", fixed=T)[[1]][2]))
bgc_data <- bgc_data[bgc_data$EARLI_Barcode %in% rownames(meta_data), ]

bgc_data_cast <- dcast(bgc_data, tax_id~EARLI_Barcode, value.var="nt_rpm")
bgc_data_cast[is.na(bgc_data_cast)] <- 0
rownames(bgc_data_cast) <- bgc_data_cast$tax_id
bgc_data_cast <- bgc_data_cast[, colnames(bgc_data_cast)!="tax_id"]

# samples with data
samples_with_data <- intersect(rownames(meta_data), colnames(bgc_data_cast))
bgc_data_cast <- bgc_data_cast[, samples_with_data]

# taxonomic data
tax_data <- as.data.frame(tax_data)
rownames(tax_data) <- tax_data$taxon
tax_data <- tax_data[rownames(bgc_data_cast), ]

# filter meta data
meta_data <- meta_data[colnames(bgc_data_cast), ]

# microbial mass
microbial_mass_data <- microbial_mass_data[rownames(meta_data), ]

# create phyloseq object
# need OTU, taxonomy, and metadata
bgc_data_cast <- otu_table(as.matrix(bgc_data_cast), taxa_are_rows = TRUE)
tax_data <- tax_table(as.matrix(tax_data))
meta_data_phyloseq <- sample_data(meta_data)
microb_data <- phyloseq(bgc_data_cast, tax_data, meta_data_phyloseq)

# aggregate at genus level
microb_data_genus <- tax_glom(microb_data, taxrank = "genus")
microb_data_glom <- list(genus=microb_data_genus)

rbm_data_cast <- rbm_data_cast[, rownames(meta_data)]

#########################
# DATA ANALYSIS
#########################
# compute in pg
meta_data$microbial_mass <- 1000*microbial_mass_data$microbial_mass

meta_data$LCA_label <- factor(meta_data$LCA_label, levels=c("Hypo", "Hyper"))

# microbial mass
meta_data$log10_microbial_mass <- log10(meta_data$microbial_mass+0.1)

pdf(paste0(results_path, "microb_mass_boxplot.pdf"), height = 3.5, width = 3)
print(ggplot(meta_data, aes(x = LCA_label, y = log10_microbial_mass,
                            fill = LCA_label)) +
        geom_boxplot() +
        geom_point() +
        ylab("Log10 Microbial Mass") +
        stat_compare_means(method = "wilcox.test") +
        scale_fill_manual(breaks=c("Hyper", "Hypo"), values = c("#EE7D31","#4472C4")) +
        theme_bw() + theme(legend.position="none",
                           text = element_text(family = "Helvetica", size=8),
                           axis.text.x = element_text(family = "Helvetica", size=8),
                           axis.text.y = element_text(family = "Helvetica", size=8),
                           axis.title.x = element_text(family = "Helvetica", size=8),
                           axis.title.y = element_text(family = "Helvetica", size=8),
                           legend.text = element_text(family = "Helvetica", size=8),
                           legend.title = element_text(family = "Helvetica", size=8),
                           strip.text.x = element_text(family = "Helvetica", size=8),
                           strip.text.y = element_text(family = "Helvetica", size=8),
                           plot.title = element_text(family = "Helvetica", size=8)))
dev.off()

# rank_names(microb_data)
tax_data %>%
  as.data.frame %>%
  filter(superkingdom == "Bacteria") %>% 
  rownames() -> bacterialist

tax_data %>%
  as.data.frame %>%
  filter(superkingdom == "Viruses") %>% 
  rownames() -> viruslist

microb_data <- list(bacteria=bacterialist, virus=viruslist)

# for bacteria at the genus level
for (microb_type in names(microb_data)){
  print(microb_type)
  microb_list <- microb_data[[microb_type]]
  
  for (rank_ in names(microb_data_glom)){
    print(rank_)

    microb_data_rank <- microb_data_glom[[rank_]]
    
    # data viz
    bact_data_rank <- prune_taxa(microb_list, microb_data_rank)
    
    # alpha diversity -> species diversity within samples
    bact_data_gen_alpha <- estimate_richness(bact_data_rank, measures=c("Shannon"))
    
    print(wilcox.test(bact_data_gen_alpha$Shannon ~ sample_data(bact_data_rank)$LCA_label))
    
    bact_data_gen_alpha$LCA_label <- factor(sample_data(bact_data_rank)$LCA_label, levels=c("Hypo", "Hyper"))
    
    print(median(bact_data_gen_alpha[bact_data_gen_alpha$LCA_label=="Hyper", "Shannon"]))
    print(median(bact_data_gen_alpha[bact_data_gen_alpha$LCA_label=="Hypo", "Shannon"]))
    
    print(median(apply(otu_table(bact_data_rank)[, sample_data(bact_data_rank)$LCA_label=="Hypo"], 2, max)/colSums(otu_table(bact_data_rank)[, sample_data(bact_data_rank)$LCA_label=="Hypo"])))
    print(median(apply(otu_table(bact_data_rank)[, sample_data(bact_data_rank)$LCA_label=="Hyper"], 2, max)/colSums(otu_table(bact_data_rank)[, sample_data(bact_data_rank)$LCA_label=="Hyper"])))
    wilcox.test(apply(otu_table(bact_data_rank)[, sample_data(bact_data_rank)$LCA_label=="Hypo"], 2, max)/colSums(otu_table(bact_data_rank)[, sample_data(bact_data_rank)$LCA_label=="Hypo"]), apply(otu_table(bact_data_rank)[, sample_data(bact_data_rank)$LCA_label=="Hyper"], 2, max)/colSums(otu_table(bact_data_rank)[, sample_data(bact_data_rank)$LCA_label=="Hyper"]))            
    
    pdf(paste0(results_path, paste0(paste(microb_type, rank_, sep="_"), "_shannon_boxplot.pdf")), height = 3.5, width = 3)
    print(ggplot(bact_data_gen_alpha, aes(x = LCA_label, y = Shannon,
                                          fill = LCA_label)) +
            geom_boxplot() +
            geom_jitter(width=0.25) +
            stat_compare_means(method = "wilcox.test") +
            scale_fill_manual(breaks=c("Hyper", "Hypo"), values = c("#EE7D31","#4472C4")) +
            theme_bw() + theme(legend.position="none",
                               text = element_text(family = "Helvetica", size=8),
                               axis.text.x = element_text(family = "Helvetica", size=8),
                               axis.text.y = element_text(family = "Helvetica", size=8),
                               axis.title.x = element_text(family = "Helvetica", size=8),
                               axis.title.y = element_text(family = "Helvetica", size=8),
                               legend.text = element_text(family = "Helvetica", size=8),
                               legend.title = element_text(family = "Helvetica", size=8),
                               strip.text.x = element_text(family = "Helvetica", size=8),
                               strip.text.y = element_text(family = "Helvetica", size=8),
                               plot.title = element_text(family = "Helvetica", size=8)))
    dev.off()
    
    # beta diversity -> compare diversities bw samples
    # add centroids to plot
    ord_df <- ordinate(prune_samples(sample_sums(bact_data_rank)>0, bact_data_rank), method = "NMDS", distance = "bray")
    scrs <- scores(ord_df)
    centroids <- setNames(aggregate(scrs$sites, by=list(meta_data[["LCA_label"]]), FUN=mean), c("LCA_label", "NMDS1", "NMDS2"))
    
    rownames(centroids) <- centroids$LCA_label
    centroids$LCA_label <- factor(centroids$LCA_label, levels=c("Hypo", "Hyper"))
    
    line_df <- data.frame(scores(ord_df)$sites)
    line_df$LCA_label <- meta_data$LCA_label
    
    line_df$NMDS1c <- centroids[line_df$LCA_label, "NMDS1"]
    line_df$NMDS2c <- centroids[line_df$LCA_label, "NMDS2"]
    
    p_nmds <- plot_ordination(prune_samples(sample_sums(bact_data_rank)>0, bact_data_rank), ord_df, color="LCA_label") +
      scale_color_manual(breaks=c("Hyper", "Hypo"), values = c("#EE7D31","#4472C4")) +
      theme_bw() + theme(text = element_text(family = "Helvetica", size=8),
                         axis.text.x = element_text(family = "Helvetica", size=8),
                         axis.text.y = element_text(family = "Helvetica", size=8),
                         axis.title.x = element_text(family = "Helvetica", size=8),
                         axis.title.y = element_text(family = "Helvetica", size=8),
                         legend.text = element_text(family = "Helvetica", size=8),
                         legend.title = element_text(family = "Helvetica", size=8),
                         strip.text.x = element_text(family = "Helvetica", size=8),
                         strip.text.y = element_text(family = "Helvetica", size=8),
                         plot.title = element_text(family = "Helvetica", size=8))
    
    p_nmds <- p_nmds + xlim(min(scrs$sites[,1]),max(scrs$sites[,1])) + ylim(min(scrs$sites[,2]),max(scrs$sites[,2]))
    p_nmds <- p_nmds + geom_point(data=centroids, size=4)
    
    p_nmds <- p_nmds + geom_segment(data=line_df, 
                                    aes(x=NMDS1c, y=NMDS2c, xend=NMDS1, yend=NMDS2, colour=LCA_label), 
                                    size=0.25, show.legend=FALSE) +
      
    pdf(paste0(results_path, paste0(paste(microb_type, rank_, sep="_"), "_nmds_boxplot.pdf")), height = 3, width = 4)
    print(p_nmds)
    dev.off()
    
    # permanova
    bact_data_gen_beta <- phyloseq::distance(prune_samples(sample_sums(bact_data_rank)>0, bact_data_rank), 
                                             method="bray")
    
    # sig
    set.seed(123)
    adonis2(formula=bact_data_gen_beta ~ as(sample_data(prune_samples(sample_sums(bact_data_rank)>0, bact_data_rank)), "data.frame")$LCA_label)
    
    # wilcox on sepsis rbm species only
    #nt_rpm
    otu_df <- rbm_data_cast
    otu_pvals <- c()
    otu_lfc <- c()
    for (otu_row in rownames(otu_df)){
      otu_row_df <- data.frame(t(otu_df[otu_row, , drop=F]))
      colnames(otu_row_df) <- c("otu_row")
      otu_row_df$LCA_label <- sample_data(bact_data_rank)[colnames(otu_df), "LCA_label"]$LCA_label
      
      otu_row_df$otu_row <- log10(otu_row_df$otu_row + 0.1)
      
      wilcox_otu_row <- wilcox.test(otu_row~LCA_label, otu_row_df)
      otu_pvals <- c(otu_pvals, wilcox_otu_row$p.value)
      
      otu_lfc <- c(otu_lfc, mean(otu_row_df$otu_row[otu_row_df$LCA_label=="Hyper"])-mean(otu_row_df$otu_row[otu_row_df$LCA_label=="Hypo"]))
    }
    otu_padj <- p.adjust(otu_pvals, "BH")
    
    wilcox_otu_df <- data.frame(cbind(rownames(otu_df), 
                                      otu_lfc,
                                      otu_pvals,
                                      otu_padj))
    colnames(wilcox_otu_df) <- c("genus", "lfc", "pval", "padj")
    wilcox_otu_df$lfc <- as.numeric(wilcox_otu_df$lfc)
    wilcox_otu_df$padj <- as.numeric(wilcox_otu_df$padj)
    
    wilcox_otu_df_sig <- wilcox_otu_df[wilcox_otu_df$padj<0.05, ]
    
    wilcox_otu_df$sig <- wilcox_otu_df$padj < 0.05
    
    # barplot log fold-change
    wilcox_otu_df <- wilcox_otu_df[order(abs(wilcox_otu_df$lfc), decreasing = T), ]
    wilcox_otu_df <- wilcox_otu_df[order(wilcox_otu_df$padj, decreasing = F), ]
    wilcox_otu_hypo_df <- wilcox_otu_df[wilcox_otu_df$lfc<0, ]
    wilcox_otu_hyper_df <- wilcox_otu_df[wilcox_otu_df$lfc>0, ]
    
    wilcox_otu_filt_df <- wilcox_otu_df[wilcox_otu_df$genus %in% c(wilcox_otu_hypo_df$genus[1:5], wilcox_otu_hyper_df$genus[1:5]), ]
    wilcox_otu_filt_df <- wilcox_otu_filt_df[order(wilcox_otu_filt_df$lfc, decreasing = T), ]
    wilcox_otu_filt_df$genus <- factor(wilcox_otu_filt_df$genus, levels=rev(wilcox_otu_filt_df$genus))
    
    wilcox_otu_filt_df$LCA_label <- wilcox_otu_filt_df$lfc
    wilcox_otu_filt_df$LCA_label[wilcox_otu_filt_df$LCA_label<0] <- 'Hypo'
    wilcox_otu_filt_df$LCA_label[wilcox_otu_filt_df$LCA_label!="Hypo"] <- 'Hyper'
    
    wilcox_otu_filt_df$alpha[wilcox_otu_filt_df$padj<0.05] <- 1
    wilcox_otu_filt_df$alpha[wilcox_otu_filt_df$padj>=0.05] <- 0.5
    
    pdf(paste0(results_path, paste0(paste(microb_type, rank_, sep="_"), "_top5_barplot.pdf")), height = 3, width = 4.5)
    print(ggplot(wilcox_otu_filt_df, aes(x=lfc, y=genus, fill=LCA_label, alpha=alpha)) + 
            geom_col() + 
            scale_fill_manual(breaks=c("Hyper", "Hypo"), values = c("#EE7D31","#4472C4")) +
            scale_alpha(range = c(0.5, 1)) +
            guides(alpha = FALSE) +
            xlab("Log2 Fold-Change")  +
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
                  plot.title = element_text(family = "Helvetica", size=8)))
    dev.off()
    
    # dominant species per sample considering sepsis-causing pathogens only
    dominant_species <- rownames(data.frame(rbm_data_cast))[apply(data.frame(rbm_data_cast), 2, which.max)]
    
    dominant_species[colnames(rbm_data_cast) %in% colnames(rbm_data_cast)[colSums(rbm_data_cast)==0]] <- NA
    
    table(dominant_species, sample_data(bact_data_rank)[colnames(otu_df), "LCA_label"]$LCA_label)
    
    otu_df <- data.frame(rbm_data_cast)
    
    # only keep maximum value per sample and replace rest with zeros
    for (otu_col in colnames(otu_df)){
      otu_df[, otu_col][otu_df[, otu_col]!=max(otu_df[, otu_col])] <- 0
    }
    
    otu_df[is.na(otu_df)] <- 0
    
    otu_df <- otu_df[rowSums(otu_df)>0, ]
    otu_df <- (otu_df > 0)
    
    otu_pvals <- c()
    otu_names <- c()
    for (otu_row in rownames(otu_df)){
      
      otu_row_df <- data.frame(t(otu_df[otu_row, , drop=F]))
      colnames(otu_row_df) <- c("otu_row")
      otu_row_df[, "LCA_label"] <- meta_data[colnames(otu_df), "LCA_label"]
      
      if (sum(otu_row_df$otu_row)>1){
        otu_row_test <- chisq.test(table(otu_row_df$otu_row, otu_row_df$LCA_label))
        
        otu_pvals <- c(otu_pvals, otu_row_test$p.value)
        otu_names <- c(otu_names, otu_row)
        
        if (otu_row_test$p.value<0.05){
          print(otu_row)
          print(otu_row_test)
        }
      }
    }
    
    otu_padj <- p.adjust(otu_pvals, "BH")
    names(otu_padj) <- otu_names
    otu_padj[otu_padj<0.05]

    otu_df_melt <- melt(otu_df)
    colnames(otu_df_melt) <- c(rank_, "Patient", "Dominant")
    
    otu_df_melt$LCA_label <- meta_data[otu_df_melt$Patient, "LCA_label"]
    
    pdf(paste0(results_path, paste0(paste(microb_type, rank_, sep="_"), "_dominant_matrix.pdf")))
    print(ggplot(otu_df_melt, aes_string(x="Patient", y=rank_, fill = "Dominant")) + 
            geom_tile(colour = "black") + 
            scale_fill_manual(breaks=c("TRUE", "FALSE"), values = c("red","white")) +
            facet_grid(~LCA_label, scales = "free_x", space = "free_x") +
            theme(strip.background = element_blank(), 
                  axis.text.x.bottom = element_blank(), 
                  axis.ticks.x = element_blank()))
    dev.off()
  }
  
  # nt_rpm sum per sample
  print(wilcox.test(colSums(rbm_data_cast) ~ meta_data[colnames(rbm_data_cast), ]$LCA_label))
  
  sum_rpm <- data.frame(sum_rpm=log10(colSums(rbm_data_cast) + 0.1), 
                        LCA_label=factor(meta_data[colnames(rbm_data_cast), ]$LCA_label, levels=c("Hypo", "Hyper")))
  
  pdf(paste0(results_path, paste0(microb_type, "_ntrpm_boxplot.pdf")), height = 3.5, width = 3)
  print(ggplot(sum_rpm, aes(x = LCA_label, y = sum_rpm,
                            fill = LCA_label)) +
          geom_boxplot() +
          geom_jitter(width=0.25) +
          stat_compare_means(method = "wilcox.test") +
          scale_fill_manual(breaks=c("Hyper", "Hypo"), values = c("#EE7D31","#4472C4")) +
          ylab("NT RPM") + 
          theme_bw() + theme(legend.position="none",
                             text = element_text(family = "Helvetica", size=8),
                             axis.text.x = element_text(family = "Helvetica", size=8),
                             axis.text.y = element_text(family = "Helvetica", size=8),
                             axis.title.x = element_text(family = "Helvetica", size=8),
                             axis.title.y = element_text(family = "Helvetica", size=8),
                             legend.text = element_text(family = "Helvetica", size=8),
                             legend.title = element_text(family = "Helvetica", size=8),
                             strip.text.x = element_text(family = "Helvetica", size=8),
                             strip.text.y = element_text(family = "Helvetica", size=8),
                             plot.title = element_text(family = "Helvetica", size=8)))
  dev.off()
}








