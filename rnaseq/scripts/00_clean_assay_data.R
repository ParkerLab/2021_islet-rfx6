#! /usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1) {
  stop("invalid argument: provide output directory.")
}

# create dir anyway without warning
dir.create(args[1], recursive = TRUE, showWarnings = FALSE)
outdir <- args[1]

suppressPackageStartupMessages({
  library(tidyverse)
  library(readxl)
})

setwd(Sys.getenv("ROOT"))

df <- readxl::read_xlsx("./data/perifusion/20210507-T2D-consolidatedPerifusionParameters-RNA-donors-updatedAUC-calc.xlsx", col_names = T)
df[df == "N/A"] <- NA

df <- df %>%
  mutate(across(-Donor, as.double))

if (sys.nframe() == 0) {
  saveRDS(file = paste0(outdir, "/perifusion-20210507.rds"), df)
}
