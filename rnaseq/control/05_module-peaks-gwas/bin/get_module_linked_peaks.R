#! /usr/bin/env Rscript

suppressPackageStartupMessages({
  library("optparse")
  library("tidyr")
  library("dplyr")
  library("glue")
  #library("GenomicRanges")
})

option_list <- list(
  make_option("--cicero",
    type = "character",
    help = "Path to annotated Cicero file"
  ),
  make_option("--module-res", dest = "module_res",
    help = "Path to module assignment results file (*.rds)"
  ),
  make_option("--out-dir-prefix", dest = "out_dir_prefix",
    help = "Prefix to the output directory where results are written"
  ),
  make_option("--celltype-label", dest = "celltype_label", default = ".",
    help = "Label to be used in the Celltype column of the link_file"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

if (any(
  is.na(opt$cicero),
  is.na(opt$module_res),
  is.na(opt$out_dir_prefix)
)) {
  message("ERROR: Please supply all arguments")
  stop(0)
}

message("info: reading data..")
# peak 1, peak 2, score, peak1_name, peak2_name
cicero_peaks <- read.table(opt$cicero, header = T)
module_res <- readRDS(opt$module_res)

#
# `%rin%` <- function(pattern, list) {
#   vapply(pattern, function(p) any(grepl(p, list)), logical(1L), USE.NAMES = FALSE)
# }

`%str_in_vec%` <- function(str, list, sep = ",") {
  vapply(
    str, function(s) {
      any(unlist(strsplit(s, split = sep, fixed = TRUE)) %in% list)
    },
    logical(1L),
    USE.NAMES = FALSE
  )
}

create_bed <- function(x, y) {
  gene_symbols <- x$symbol

  df <- subset(
    cicero_peaks,
    peak_1_name %str_in_vec% gene_symbols | peak_2_name %str_in_vec% gene_symbols
  )

  df$module <- y$module
  df
}

message("info: fetching module linked peaks..")
df <- module_res %>%
  dplyr::group_by(module) %>%
  dplyr::group_map(~ create_bed(.x, .y)) %>%
  do.call(rbind, .)

dir.create(opt$out_dir_prefix, recursive = T, showWarnings = F)

message("info: writing peaks..")
df %>%
  group_by(module) %>%
  group_walk(
    ~ readr::write_tsv(
      unique(as.data.frame(stringr::str_split(unique(.x$peak_1, .x$peak_2), "_", simplify = T))),
      paste0(opt$out_dir_prefix, "/", "Module-", .y$module, ".bed"),
      col_names = F
    )
  )

readr::write_tsv(
  unique(as.data.frame(stringr::str_split(
    unique(cicero_peaks$peak_1, cicero_peaks$peak_2), "_", simplify = T))
  ),
  paste0(opt$out_dir_prefix, "/", "Module-ALL", ".bed"),
  col_names = F
)

# module_ids <- c(unique(module_res$module), "ALL")

# annot_df <- data.frame(
#   Index = 1:length(module_ids),
#   Annotation = glue("Module-{module_ids}"),
#   Celltype = opt$celltype_label,
#   Tissue = "Pancreas",
#   Type = "Peaks",
#   Path = glue("{opt$out_dir_prefix}/Module-{module_ids}.bed")
# )

# readr::write_tsv(annot_df, glue("{opt$out_dir_prefix}/annotation_link_file.txt"))
message("info: done.")
