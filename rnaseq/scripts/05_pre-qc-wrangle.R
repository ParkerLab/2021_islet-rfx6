#! /usr/bin/env Rscript

suppressPackageStartupMessages({
  library(glue)
  library(tidyverse)
  library(cowplot)
  library(rtracklayer)
  library(gghalves)
  library(ggbeeswarm)
})

theme_set(theme_custom())
setwd("~/analyses/2020-01_vanderbilt_rna/")

# Get sample.IDs
samples <- list.dirs("work/2020-01-26_remap/qorts/pre-prune/", full.names = F, recursive = F)

# Read Donor.ID and Sample.ID mapping --> convert to long format
rna_files <- read.table("sample_info/rna_samples.txt", header = T)
rna_files <- rna_files %>% pivot_longer(names_to = "Type", cols = -Donor.ID, values_to = "Sample.ID") %>%
  filter(!is.na(Sample.ID))

# Read co-variate file
covs <- read.table("sample_info/covariates.txt", header = T, sep = '\t')

# Read RIN values
rin <- read.table("sample_info/rin.txt", sep = '\t', header = T)
rin$RIN <- as.numeric(rin$RIN)

# Read TIN values
tin <- do.call(rbind, lapply(samples, function(x) {
  read.table(paste0("work/2020-01-26_remap/tin/", x, "/Aligned.sortedByCoord.out.summary.txt"), header = T,
             col.names = c("Sample.ID", "TIN_mean", "TIN_median", "TIN_sdev")) %>%
    mutate(Sample.ID = x)
}))


## PCR batches
pcr_batches <- read.table("sample_info/pcr_batches.txt", header = T)

# Read QoRTs
qorts <- read.table("work/2020-02-14_qc_summary/qorts.txt", header = T)

# Read PCA
pca_df <- read.table("work/2020-02-14_qc_summary/pca_protein-coding_tmm.txt", header = T)

# Read
counts_df <- read.table("work/2020-02-14_qc_summary/counts_summary.txt", header = T)

# Read FACS cell coutn
facs_count <- read.table("sample_info/sample_cell-counts.txt", header = T)

# Check what data is missing in above
## dplyr::select(all_df, Type, Sample.ID, Num.cells) %>% filter(Type != 'Whole.Islet.RNA') %>% filter(is.na(Num.cells))

## Combine all
all_df <- left_join(rna_files, covs, by = c("Donor.ID" = "Donor")) %>%
  left_join(., qorts) %>%
  left_join(., rin) %>%
  left_join(., tin) %>%
  left_join(., counts_df) %>%
  left_join(., pcr_batches) %>%
  left_join(., facs_count) %>%
  left_join(., dplyr::select(pca_df, Sample.ID, PC1, PC2, PC3, PC4)) %>%
  mutate(
    Age = ifelse(Age.Units == "years", Age, Age/12),
    Age.Units = 'years',
    Sequencing.Batch = PCR.group,
    C.peptide.ng.mL = as.numeric(C.peptide.ng.mL),
    Date.Islets.Received = ifelse(
      is.na(parse_date(Date.Islets.Received, "%m/%d/%Y")),
      parse_date(Date.Islets.Received, "%m/%d/%y"),
      parse_date(Date.Islets.Received, "%m/%d/%Y"))
  ) %>%
  mutate(
    Islet.Days.Since.First = as.numeric(Date.Islets.Received - min(Date.Islets.Received))
  ) %>%
  dplyr::select(-gbcbias_5p, -gbcbias_3p, -insertionEventCt, -TIN_median, -TIN_sdev, -Date.Islets.Received, -Time.on.Ventilator..Days.)


write.table(all_df, file = 'work/2020-02-14_qc_summary/all_df.txt', sep = '\t', quote = F, row.names = F)

## Plots

ggplot(all_df, aes(' ', gbcbias_5p_3p)) +
  geom_quasirandom(width = .25, aes(col = Disease.Status), size = 1) + geom_boxplot(width = .1, notch = T, alpha = .2, col = 'black', outlier.color = 'NA') +
  ggrepel::geom_text_repel(aes(label = ifelse(gbcbias_5p_3p < .95, Sample.ID, '')), nudge_x = .2) +
  labs(x = '', y = "5'/3' bias (high expression genes)",
       title = "**5'-3' coverage bias** <br> <span style='font-size:11pt'>Lower is poor</span>") +
  theme(plot.title = element_markdown(lineheight = 1.1)) + facet_wrap(~Type)
save_plot("figures/qc/gbc-5p-3p-bias", base_height = 6, base_width = 8)

# Check Xist/Sry expression with sex
ggplot(all_df, aes(xist_exp, chrY_mean_exp)) + geom_point(aes(col = Sex)) +
  labs(x = "XIST", y = "Mean expression of ChrY", title = "Sex-linked gene expression")
save_plot("figures/qc/xist_by_meanChrY_sex", base_height = 4, base_width = 5)


# Check CGD K-S Statistic Z-scores
ggplot(all_df, aes(x = '', y=cgd_zscore)) +
  geom_quasirandom(width = .25, aes(col = Disease.Status), size = 1) + geom_boxplot(width = .1, notch = T, alpha = .2, col = 'black', outlier.color = 'NA') +
  ggrepel::geom_text_repel(aes(label = ifelse(cgd_zscore > 2, Sample.ID, '')), nudge_x = 2) +
  labs(x = '', y = "Z-score of K-S Statistic rel. to median",
       title = "**Cumulative Gene Diveristy** <br> <span style='font-size:11pt'>Higher is poor</span>") +
  theme(plot.title = element_markdown(lineheight = 1.1)) + facet_wrap(~Type)
save_plot("figures/qc/cgd_z_score", base_height = 6, base_width = 8)

## Read pairs (unique genes)
ggplot(all_df, aes(x = 'Samples', y=ReadPairs_UniqueGene)) +
  geom_quasirandom(width = .1, col = gray(.1), size = 1) + geom_boxplot(width = .1, notch = T, alpha = .2, col = 'red', outlier.color = 'NA') +
  ggrepel::geom_text_repel(aes(label = ifelse(ReadPairs_UniqueGene < .5e7 | ReadPairs_UniqueGene > 2.9e7, Sample.ID, '')), nudge_x = 2) +
  scale_y_continuous(labels = scales::unit_format(scale = 1e-6, unit = "M")) +
  labs(x='', y = 'Read Pairs (mapping to unique genes)')
save_plot("figures/qc/readpairs_unique-gene", base_height = 6, base_width = 4)

# Insert Size
ggplot(all_df, aes(x = 'Samples', y=InsertSize_Median)) +
  geom_quasirandom(width = .1, col = gray(.1), size = 1) + geom_boxplot(width = .1, notch = T, alpha = .2, col = 'red', outlier.color = 'NA') +
  #ggrepel::geom_text_repel(aes(label = ifelse(InsertSize_Median < 200, Donor.ID, '')), size = 3) +
  #scale_y_continuous(labels = scales::unit_format(scale = 1e-6, unit = "M"))
  labs(x='', y = 'Insert Size (median)') +
  facet_wrap(~Type)
save_plot("figures/qc/insert_size_by_type", base_height = 5, base_width = 7)

# Insertion event ct
#ggplot(all_df, aes(x = 'Samples', insertionEventCt)) +
#  geom_quasirandom(width = .1, col = gray(.1), size = 1) + geom_boxplot(width = .1, notch = T, alpha = .2, col = 'red', outlier.color = 'NA') +
#  ggrepel::geom_text_repel(aes(label = ifelse(insertionEventCt > 3e5, Sample.ID, '')), nudge_x = 2) +
#  labs(x='', y = 'Insertion Event Count')
#save_plot("figures/qc/insertion-evt-count", base_height = 6, base_width = 4)

## RIN
ggplot(all_df, aes(x = 'Samples', RIN)) +
  geom_quasirandom(width = .1, col = gray(.1), size = 1) + geom_boxplot(width = .1, notch = T, alpha = .2, col = 'red', outlier.color = 'NA') +
  ggrepel::geom_text_repel(aes(label = ifelse(RIN < 5, Sample.ID, '')), nudge_x = 2) +
  labs(x='', y = 'RIN')
save_plot("figures/qc/rin", base_height = 6, base_width = 4)

## TIN
ggplot(all_df, aes(x = 'Samples', TIN_mean)) +
  geom_quasirandom(width = .1, col = gray(.1), size = 1) + geom_boxplot(width = .1, notch = T, alpha = .2, col = 'red', outlier.color = 'NA') +
  ggrepel::geom_text_repel(aes(label = ifelse(TIN_mean < 65, Sample.ID, '')), nudge_x = 2) +
  labs(x='', y = 'TIN')
save_plot("figures/qc/tin", base_height = 6, base_width = 4)

## RIN vs TIN
ggplot(data = NULL, aes(all_df$RIN, all_df$TIN_mean)) + geom_point() + labs(x = "RIN", y = "TIN (mean)") +
  ggrepel::geom_text_repel(aes(label = ifelse(all_df$TIN_mean < 65 | all_df$RIN < 5, all_df$Sample.ID, '')))
save_plot("figures/qc/rin-vs-tin", base_height = 5, base_width = 5)

pca_plot <- ggplot(all_df, aes(PC1, PC2)) + geom_point(aes(col = Disease.Status))
save_plot("figures/qc/pca_disease-status", base_height = 5, base_width = 6.5)
