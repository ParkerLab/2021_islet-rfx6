#! /usr/bin/env Rscript

# // Copyright (c) 2021 Vivek Rai
# //
# // This software is released under the MIT License.
# // https://opensource.org/licenses/MIT

# FIXME: IN Progress / Script not complete

suppressPackageStartupMessages({
  library("optparse")
  library("tidyr")
  library("dplyr")
  library("glue")
})

option_list <- list(
  make_option("--input",
    type = "character",
    help = "Path to directory containing bed files"
  ),
  make_option("--outfile", dest = "outfile",
    help = "Path to the output file where annotation link is written"
  ),
  make_option("--celltype", dest = "celltype_label", default = ".",
    help = "Celltype column value"
  ),
  make_option("--tissue", dest = "tissue", default = "Pancreas",
    help = "Tissue column value"
  ),
  make_option("--type", dest = "col_type", default = "Peaks",
    help = "Type column value"
  )
)

opt <- parse_args(OptionParser(option_list = option_list))

if (any(
  is.na(opt$input),
  is.na(opt$outfile)
)) {
  message("ERROR: Please supply all arguments")
  stop(0)
}

opt$input = file.path(normalizePath(dirname(opt$input)), basename(opt$input))

file_list <- list.files(opt$input, "*.bed", full.names = T)

annot_df <- data.frame(
  Index = 1:length(file_list),
  Annotation = glue("{basename(tools::file_path_sans_ext(file_list))}"),
  Celltype = opt$celltype_label,
  Tissue = opt$tissue,
  Type = opt$col_type,
  Path = glue("{file_list}")
)

readr::write_tsv(annot_df, opt$outfile)
message("info: done.")