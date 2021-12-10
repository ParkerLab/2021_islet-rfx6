#! /usr/bin/env Rscript

'Prepare counts.

Usage:
  03_counts.R <root> <outdir>
  03_counts.R (-h | --help)

Options:
  -h --help     Show this screen.

' -> doc

library(docopt)
args <- docopt(doc)

writeLines(paste0("info: using output directory: ", outdir))
dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  library(glue)
  library(tidyverse)
  library(cowplot)
  library(rtracklayer)
})

theme_set(cowplot::theme_cowplot())

setwd(args$root)
gtf <- readRDS('../data/gtf/gencode.v19.annotation.gtf.rds')

## Read all counts
counts <- data.table::fread("../work/2020-01-26_remap/collated_counts/gene_counts_all.txt", header = T)

################### Get a few summary statistics from Counts matrix #######################
xist_exp = tpm_counts[which(counts$Geneid == "ENSG00000229807.5"), ]  # XIST gene
sry_exp <- tpm_counts[which(counts$Geneid == "ENSG00000184895.6"), ]  # SRY gene

# mean Y-chromosome exp
y_chr_genes <- genes %>% filter(seqnames == 'chrY') %>% select(gene_id) %>% .[[1]]
mean_y_exp <- apply(tpm_counts[which(y_chr_genes %in% counts$Geneid), ], 2, mean)

chrM_genes <- genes %>% filter(seqnames == "chrM") %>% select(gene_id) %>% .[[1]]
mean_M_exp <- apply(tpm_counts[which(chrM_genes %in% counts$Geneid), ], 2, mean)

out_df <- data.frame(Sample.ID = colnames(tpm_counts),
                     sry_exp = sry_exp,
                     xist_exp = xist_exp,
                     chrY_mean_exp = mean_y_exp,
                     chrM_mean_exp = mean_M_exp)

#write.table(out_df, "work/2020-02-14_qc_summary/counts_summary.txt", row.names = F, quote = F, sep = "\t")
write.table(out_df, paste0(args$outdir, "/", "counts_summary.txt"), row.names = F, quote = F, sep = "\t")
saveRDS(counts, file = paste0(args$outdir, "/", "counts.rds"))

protein_coding_genes <- as.data.frame(
  subset(gtf, type == 'gene' & gene_type == 'protein_coding' & seqnames != "chrM" & seqnames != "chrY")
)[, c(1:5, 10, 14)]

pc_df <- as.data.frame(counts[counts$Geneid %in% protein_coding_genes$gene_id, ])
pc_df$Geneid <- gsub("\\.[0-9]*$", "", pc_df$Geneid)

rownames(pc_df) <- pc_df$Geneid
pc_df <- subset(pc_df, select = -c(Geneid, Length))

saveRDS(pc_df, file = paste0(args$outdir, "/", "counts_protein-coding.rds"))
