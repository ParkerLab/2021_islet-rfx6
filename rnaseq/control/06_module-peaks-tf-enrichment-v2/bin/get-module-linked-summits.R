#! /usr/bin/env Rscript

suppressPackageStartupMessages({
  library("optparse")
  library("tidyr")
  library("dplyr")
  library("glue")
  library("magrittr")
  library("GenomicRanges")
})

option_list <- list(
  make_option("--bed",
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
  is.null(opt$bed),
  is.null(opt$module_res),
  is.null(opt$out_dir_prefix)
)) {
  message("ERROR: Please supply all arguments")
  stop(0)
}

message("info: reading data..")
# peak 1, peak 2, score, peak1_name, peak2_name
annotated_peaks <- read.table(opt$bed, header = F)

module_res <- readRDS(opt$module_res)

create_bed <- function(x, y) {
  gene_symbols <- x$symbol
  df <- subset(annotated_peaks, V9 %in% gene_symbols)
  GenomicRanges::reduce(
    GenomicRanges::makeGRangesFromDataFrame(
      df[, c("V1", "V2", "V3")],
      seqnames.field = "V1",
      start.field = "V2",
      end.field = "V3",
      starts.in.df.are.0based = T
    )
  ) %>% as.data.frame() %>%
    mutate(module = y$module)
}

message("info: fetching module linked peaks..")
df <- module_res %>%
  dplyr::group_by(module) %>%
  dplyr::group_map(~ create_bed(.x, .y)) %>%
  do.call(rbind, .)

dir.create(opt$out_dir_prefix, recursive = T, showWarnings = F)

message("info: writing module peaks..")
df %>%
  group_by(module) %>%
  group_walk(
    ~ rtracklayer::export.bed(
        GenomicRanges::makeGRangesFromDataFrame(.x),
        con = paste0(opt$out_dir_prefix, "/", "Module-", .y$module, ".bed")
      )
  )

# Write peaks from ALL modules
rtracklayer::export.bed(
  GenomicRanges::reduce(GenomicRanges::makeGRangesFromDataFrame(df)),
  con = paste0(opt$out_dir_prefix, "/", "Module-ALL", ".bed")
)

message("info: writing peaks not linked to modules..")
# Write peaks that are NOT linked to modules
# Process: This will be the peaks in annotated_peaks that are not in df
peaks_gr <- GenomicRanges::makeGRangesFromDataFrame(
  annotated_peaks[, 1:3], start.field = "V2", end.field = "V3", seqnames.field = "V1"
)

df_gr <- GenomicRanges::makeGRangesFromDataFrame(df)

overlaps <- GenomicRanges::findOverlaps(peaks_gr, df_gr)
peaks_not_in_modules <- peaks_gr[-S4Vectors::queryHits(overlaps)]

rtracklayer::export.bed(
  GenomicRanges::reduce(peaks_not_in_modules),
  con = paste0(opt$out_dir_prefix, "/", "Module-NO-MODULE", ".bed")
)


# Write "control" peaks / peaks shuffled across modules but binned in same
# frequency
# n_peaks <- nrow(df)
# module_peak_counts <- df %>%
#   group_by(module) %>%
#   summarise(peaks = n())
# 
# lapply(1:nrow(module_peak_counts), function(row_idx) {
#   module_id <- module_peak_counts[row_idx, 1]
#   peak_count <- module_peak_counts[row_idx, 2]
# 
#   diff_module_idx <- which(df$module != module_id)
#   sample_idx <- sample(diff_module_idx, peak_count)
# 
#   rtracklayer::export.bed(
#     GenomicRanges::makeGRangesFromDataFrame(sample_df[sample_idx, ]),
#     con = paste0(opt$out_dir_prefix, "/", glue("Module-{sdf}-PERM"), ".bed")
#   )
# 
#   df <- df[-sample_idx, ]
# })

#median_peak_size <- as.integer(median(module_peak_counts$peaks))
#
#df %>%
#  group_by(module) %>%
#  group_walk(
#    ~ rtracklayer::export.bed(
#      GenomicRanges::makeGRangesFromDataFrame(
#        na.omit(.x[sample(median_peak_size), ]),
#      ),
#      con = paste0(opt$out_dir_prefix, "/", "Module-", .y$module, glue("-SMPL.{median_peak_size}"), ".bed")
#    )
#  )

message("info: done.")
