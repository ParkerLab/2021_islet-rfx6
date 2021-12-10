#! /usr/bin/env Rscript

root <- "/home/vivekrai/analyses/2020-01_vanderbilt_rna/"

suppressPackageStartupMessages({
  library(tidyverse)
  library(cowplot)
  library(glue)
})

cov_df <- read.table(glue("{root}/work/2020-02-14_qc_summary/all_df.txt"), sep = "\t", header = T)

## Include post QC decisions
exclude_samples <- c(
  "MG-0066_i-ND-15m", "MG-0124_b-ND-3", "MG-0054_a-ND-11", "MG-0055_b-ND-11",
  "MG-0045_a-T2D-51", "MG-0104_a-ND-19", "MG-0008_b-ND-10m"
)

## Include information about functional data availability
cov_df <- cov_df %>% filter(!Sample.ID %in% exclude_samples)

## Make venn diagrams indicating how many libraries met criteria
{
  pdf(glue("{root}figures/qc/03-24_vennDiagram-QC-pass.pdf"), width = 10, height = 3)
  par(mfrow = c(1, 3))
  limma::vennDiagram(limma::vennCounts(df[df$Type == "Alpha.cell.RNA", 4:6]), circle.col = 1:4)
  title("Alpha")
  limma::vennDiagram(limma::vennCounts(df[df$Type == "Beta.cell.RNA", 4:6]), circle.col = 1:4)
  title("Beta")
  limma::vennDiagram(limma:vennCounts(df[df$Type == "Whole.Islet.RNA", 4:6]), circle.col = 1:4)
  title("Whole Islet")
  dev.off()
}

if (sys.nframe() == 0) {
  dir.create(glue("{root}/work/2020-03-25_post-qc"))

  saveRDS(cov_df, file = glue("{root}/work/2020-03-25_post-qc/sample_info.rds"))
  cov_df <- cov_df %>% filter(-starts_with("PC")) # remove bulk PCS

  pc_counts <- readRDS(glue("{root}/work/2020-02-14_qc_summary/counts_protein-coding.rds"))
  pc_counts <- as.matrix(pc_counts[, !(names(pc_counts) %in% exclude_samples)])
  saveRDS(pc_counts, file = glue("{root}/work/2020-03-25_post-qc/counts_protein-coding.rds"))

  # saveRDS(cov_df %>% filter(Type == "Beta.cell.RNA"), file = "work/2020-03-25_post-qc-covariates/beta_sample_info.rds")
  # saveRDS(cov_df %>% filter(Type == "Alpha.cell.RNA"), file = "work/2020-03-25_post-qc-covariates/alpha_sample_info.rds")
  # saveRDS(cov_df %>% filter(Type == "Whole.Islet.RNA"), file = "work/2020-03-25_post-qc-covariates/islet_sample_info.rds")
  # saveRDS(cov_df %>% filter(grepl("Juvenile", Disease.Status)), file = "work/2020-03-25_post-qc-covariates/juvenile_sample_info.rds")
  # saveRDS(cov_df %>% filter(!grepl("Juvenile", Disease.Status)), file = "work/2020-03-25_post-qc-covariates/adult_sample_info.rds")
}
