#! /usr/bin/env Rscript
library(optparse)

option_list <- list(
  make_option(c("--qorts-dir"),
              action = "store", dest = "qorts_dir", type = "character",
              help = "[Required] Path to the QoRTs output directory"
  ),
  make_option(c("--out-dir"),
              action = "store", dest = "out_dir",
              help = "[Required] Output directory"
  ),
  make_option(c("--covariates-file"),
              action = "store", type = "character", dest = "covariates_file",
              help = "[Required] Covariate file containing donor information"
  ),
  make_option(c("--donor-sample-file"),
              action = "store", type = "character", dest = "donor_sample_file",
              help = "[Required] File containing donor and library information pair"
  )
)

option_parser <- OptionParser(
  usage = "usage: Rscript %prog [options]",
  option_list = option_list, add_help_option = T
)

opts <- parse_args(option_parser)

if(length(opts) != 5) {
  parse_args(option_parser, args = c("--help"))
  stop()
}

suppressPackageStartupMessages({
  library(glue)
  library(tidyverse)
  library(cowplot)
  library(QoRTs)
  library(rtracklayer)
})

theme_set(cowplot::theme_cowplot())

# Get sample.IDs
samples <- list.dirs(opts$qorts_dir, full.names = F, recursive = F)

# Read Donor.ID and Sample.ID mapping --> convert to long format
rna_files <- read.table(opts$donor_sample_file, header = T)

# Read co-variate file
covs <- read.table(opts$covariates_file, header = T, sep = '\t')

# Make a Venn diagram
function() {
  tmp <- rna_files[,]
  tmp[!is.na(tmp)] <- 1L
  tmp[is.na(tmp)] <- 0L
  tmp$Donor.ID <- rna_files$Donor.ID

  dir.create(paste0(opts$out_dir, "figures"), recursive = T, showWarnings = F)

  pdf(paste0(opts$out_dir, "figures", "libraries_venn.pdf", sep = "/"), width = 6, height = 6)
  limma::vennDiagram(limma::vennCounts(data.matrix(tmp[, 2:3])), circle.col = 1:4)
  dev.off()
}()

# Pivot to longer dataframe
rna_files <- rna_files %>%
  pivot_longer(names_to = "Type", cols = -Donor.ID, values_to = "Sample.ID") %>%
  filter(!is.na(Sample.ID)) %>%
  left_join(., covs, by = c("Donor.ID" = "Donor"))

# Create QoRTs metadata/decoder information
decoder_obj <- rna_files %>%
  dplyr::select(Donor.ID, Sample.ID, Type, Disease.Status) %>%
  dplyr::mutate(Type = gsub("\\.RNA", "", Type)) %>%
  dplyr::rename(
    unique.ID = Sample.ID,
    sample.ID = Donor.ID,
    lane.ID = Type,
    group.ID = Disease.Status
  ) %>%
  as.data.frame()

## Load all the QC information
message("info: reading QoRTs data..")
pre_prune <- read.qc.results.data(paste0(opts$qorts_dir, "/"), decoder = decoder_obj)

saveRDS(pre_prune, file = paste0(opts$out_dir, "pre-prune.rds", sep = "/"))

## Generate QoRTs plot for pre-pruned libraries
figure_dir <- paste0(opts$out_dir, "figures", "qorts-pre-hg19", sep = "/")
dir.create(figure_dir, recursive = T, showWarnings = F)

message("info: generating QoRTs basic figures..")
#makeMultiPlot.all(pre_prune, outfile.dir = figure_dir, plot.device.name = "pdf", rasterize.large.plots = T)

## Make certain diagnostic figures separately figures
plotter.colorByGroup <- build.plotter.colorByGroup(pre_prune)
plotter.colorBySample <- build.plotter.colorBySample(pre_prune)
plotter.colorByType <- build.plotter.colorByLane(pre_prune)

makeFig <- function(fun, ...) {
  par(mfrow=c(1, 3))
  fun(plotter.colorBySample, ...)
  fun(plotter.colorByGroup, ...)
  makePlot.legend.over("topright", plotter.colorByGroup)
  fun(plotter.colorByType, ...)
  makePlot.legend.over("topright", plotter.colorByType)
}

message("info: creating select QoRTs figures..")

select_figure_dir <- paste0(opts$out_dir, "figures", "qorts-select-figures", sep = "/")
dir.create(select_figure_dir, recursive = T, showWarnings = F)

pdf(paste0(select_figure_dir, 'qorts-insert-size.pdf', sep = '/'), width = 8, height = 4)
makeFig(makePlot.insert.size, xlim = c(0, 1000))
dev.off()

pdf(paste0(select_figure_dir, 'gc-bias.pdf', sep = '/'), width = 8, height = 4)
makeFig(makePlot.gc)
dev.off()

pdf(paste0(select_figure_dir, 'clipping-rates.pdf', sep = '/'), width = 8, height = 4)
makeFig(makePlot.clipping)
dev.off()

pdf(paste0(select_figure_dir, 'cigar-ops-ins.pdf', sep = '/'), width = 8, height = 4)
makeFig(makePlot.cigarOp.byCycle, "Ins")
dev.off()

pdf(paste0(select_figure_dir, 'cigar-ops-del.pdf', sep = '/'), width = 8, height = 4)
makeFig(makePlot.cigarOp.byCycle, "Del")
dev.off()

pdf(paste0(select_figure_dir, 'gene-cdf', sep = '/'), width = 8, height = 4)
makeFig(makePlot.gene.cdf)
dev.off()

pdf(paste0(select_figure_dir, 'mapping-rates.pdf', sep = '/'), width = 8, height = 4)
makeFig(makePlot.gene.assignment.rates)
dev.off()

pdf(paste0(select_figure_dir, 'phred-qual.pdf', sep = '/'), width = 8, height = 4)
makeFig(makePlot.qual.pair, "min")
dev.off()

pdf(paste0(select_figure_dir, 'gene-body-cdf.pdf', sep = '/'), width = 10, height = 4)
makeFig(makePlot.genebody)
dev.off()

pdf(paste0(select_figure_dir, 'gene-body-um-quartile.pdf', sep = '/'), width = 10, height = 4)
makeFig(makePlot.genebody, geneset="50-75")
dev.off()

## Get alignment soft clipping rate
get_soft_clipping_rate <- function(qorts_obj, read = "r1") {
  if (read == "r1") {
    tag = "cigarOpDistribution.byReadCycle.R1"
  } else {
    tag = "cigarOpDistribution.byReadCycle.R2"
  }
  if(!"qc.data" %in% slotNames(qorts_obj) || !tag %in% names(qorts_obj@qc.data)) {
    print("qc.data not found in input object. Is it a QoRTs object?")
    return(NULL)
  }

  do.call(cbind, lapply(qorts_obj@qc.data[[tag]], function(df) {
    rowSums(df[, c("S_S", "S_M", "S_E", "S_B")])/rowSums(df)
  }))
}

## Get 5' bias, 3' bias, and 5'-3' bias; using 10-%ile and 90-%ile values
get_gbc_bias <- function(qorts_obj, gene_category = "X4.high") {
  if(!"qc.data" %in% slotNames(qorts_obj) || !"geneBodyCoverage.pct" %in% names(qorts_obj@qc.data)) {
    print("qc.data not found in input object. Is it QoRTs object?")
    return(NULL)
  }
  do.call(rbind, lapply(qorts_obj@qc.data[["geneBodyCoverage.pct"]], function(x) {
    mean_cov <- mean(x[[gene_category]])
    data.frame("gbcbias_5p" = x[4, gene_category]/mean_cov, "gbcbias_3p" = x[36, gene_category]/mean_cov) %>%
      mutate("gbcbias_5p_3p" = gbcbias_5p / gbcbias_3p)
  })) %>% rownames_to_column(var = "Sample.ID")
}

gbc_plot <- ggplot(get_gbc_bias(pre_prune), aes('', gbcbias_5p_3p)) +
  geom_boxplot(width = .1, notch = T, alpha = .2, col = 'red', outlier.color = 'NA')  +
  ggbeeswarm::geom_quasirandom() +
  ggrepel::geom_text_repel(aes(label = ifelse(abs(gbcbias_5p_3p) > 2.5, Sample.ID, ''))) +
  labs(x = '', y = "5' - 3' bias of gene body coverage")

pdf(paste0(select_figure_dir, "gbc_bias.pdf", sep = "/"), height = 6, width = 3)
print(gbc_plot)
dev.off()

## Compute K-S statistic for Gene Body Coverage deviation from Median Counts (all genes)
get_cdg <- function(qorts_obj) {
  if(!"calc.data" %in% slotNames(qorts_obj) || !"LANEBAM_GENE_CDF" %in% names(qorts_obj@calc.data)) {
    print("calc.data not found in input object. Is it QoRTs object?")
    return(NULL)
  }
  apply(do.call(cbind, qorts_obj@calc.data[["LANEBAM_GENE_CDF"]]), 2, function(x) x/tail(x, 1))
}

# Get cumulative gene diversity
cgd.df <- get_cdg(pre_prune)

### Compute CGD per group
cgd_zscore <- do.call(rbind, lapply(unique(rna_files$Type), function(x) {
  cgd.df.subset <- cgd.df[, rna_files %>% filter(Type == x) %>% dplyr::select(Sample.ID) %>% .[[1]] ]
  median_group <- apply(cgd.df.subset, 1, median)

  data.frame(cgd_zscore = apply(cgd.df.subset, 2, function(x) ks.test(x, median_group)$statistic)) %>%
      rownames_to_column(var = "Sample.ID") %>%
      mutate(cgd_zscore = scale(cgd_zscore))
}))

cgd_plot <- ggplot(cgd_zscore, aes('', cgd_zscore)) +
  geom_boxplot(width = .1, notch = T, alpha = .2, col = 'red', outlier.color = 'NA')  +
  ggbeeswarm::geom_quasirandom() +
  ggrepel::geom_text_repel(aes(label = ifelse(abs(cgd_zscore) > 2.5, Sample.ID, ''))) +
  labs(x = '', y = "Z-score of K-S statistic (wrt. median of each type)")

pdf(paste0(select_figure_dir, "cgd-samples.pdf", sep = "/"), height = 6, width = 3)
print(cgd_plot)
dev.off()


# Extract other relevant metrics from the QoRTs summary table
metrics_of_interest <- c(
  "AVG_GC",
  "InsertSize_Median",
  "ReadPairs_UniqueGene",
  "ReadPairs_NoGene_Intron",
  "READ_PAIR_OK",
  "insertionEventCt"
)

summary_pre <- get.summary.table(pre_prune)

message("info: writing summary data..")
df <- data.frame(t(summary_pre[which(rownames(summary_pre) %in% metrics_of_interest), ])) %>%
  rownames_to_column(var = "Sample.ID") %>%
  mutate(Sample.ID = str_replace_all(Sample.ID, '\\.', '-')) %>%
  left_join(bias) %>%
  left_join(ks_stats_zscore)

# Write the final summary table
write.table(df, file = paste0(opts$out_dir, "qorts_summary.txt", sep = "/"), quote = F, sep = "\t", row.names = F)
