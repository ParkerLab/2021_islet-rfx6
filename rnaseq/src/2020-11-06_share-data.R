#! /usr/bin/env Rscript

analysis_dir <- "/lab/work/vivekrai/2020-01_vanderbilt_rna/work/dge-analysis/2021-02-03_14-52"
out_dir <- "/home/vivekrai/analyses/2020-01_vanderbilt_rna/work/dge-analysis/2021-02-05_normalized-data-share"

beta_k <- 6
alpha_k <- 4
islet_k <- 4

k_vals <- c(
    "Beta.cell.RNA" = 6,
    "Alpha.cell.RNA" = 4,
    "Whole.Islet.RNA" = 4
)

for (celltype in names(k_vals)) {
    path <- glue::glue("{analysis_dir}/{celltype}-Disease-10-0.25/ruv-seq_de.rds")
    print(glue::glue("info: reading {path}..."))

    # sanity check
    stopifnot(file.exists(path))

    data <- readRDS(path)

    saveRDS(
        data[[k_vals[celltype]]],
        file = glue::glue("{out_dir}/{celltype}-{k_vals[celltype]}.rds"),
    )
}